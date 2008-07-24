module tao_lm_optimizer_mod

use tao_mod
use tao_dmerit_mod
use tao_top10_mod
use tao_var_mod
use input_mod

use nrtype; use nrutil, only : assert_eq,diagmult
use nr, only : covsrt,gaussj

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_lm_optimizer (abort)
!
! Subrutine to minimize the merit function by varying variables until
! the "data" as calculated from the model matches the measured data.
! 
! This subroutine is a wrapper for the mrqmin routine of Numerical Recipes.
! See the Numerical Recipes writeup for more details.
! 'lm' stands for Levenburg - Marquardt. Otherwise known as LMDIF. 
!
! Output:
!   abort -- Logical: Set True if an user stop signal detected.
!-

subroutine tao_lm_optimizer (abort)

implicit none

type (tao_universe_struct), pointer :: u

real(rp), allocatable, save :: y(:), weight(:), a(:)
real(rp), allocatable, save :: y_fit(:)
real(rp), allocatable, save :: dy_da(:, :)
real(rp), allocatable, save :: var_value(:), var_weight(:), var_meas_value(:)
real(rp) a_lambda, chi_sq, merit0

integer i, j, k
integer n_data, n_var

logical :: finished, init_needed = .true.
logical limited, limited2, abort

character(20) :: r_name = 'tao_lm_optimizer'
character(80) line
character(1) char

! setup

a_lambda = -1
abort = .false.

merit0 = tao_merit()

call tao_get_vars (var_value, var_weight = var_weight, var_meas_value = var_meas_value)
n_var = size(var_value)

n_data = n_var
do i = lbound(s%u, 1), ubound(s%u, 1)
  if (.not. s%u(i)%is_on) cycle
  n_data = n_data + count(s%u(i)%data(:)%useit_opt .and. s%u(i)%data(:)%weight /= 0)
enddo

if (allocated(y)) deallocate(y, weight, a, tao_com%covar, tao_com%alpha, y_fit, dy_da)
allocate (y(n_data), weight(n_data), y_fit(n_data))
allocate (a(n_var), tao_com%covar(n_var,n_var), tao_com%alpha(n_var,n_var))
allocate (dy_da(n_data, n_var))

! init a and y arrays

a = var_value
y = 0
weight(1:n_var) = var_weight

k = n_var
do j = lbound(s%u, 1), ubound(s%u, 1)
  u => s%u(j)
  if (.not. u%is_on) cycle
  do i = 1, size(u%data)
    if (.not. u%data(i)%useit_opt) cycle
    if (u%data(i)%weight == 0) cycle
    k = k + 1
    weight(k) = u%data(i)%weight
  enddo
enddo

call tao_dModel_dVar_calc (s%global%derivative_recalc)

! run optimizer mrqmin from Numerical Recipes.

finished = .false.
call out_io (s_blank$, r_name, '   Loop      Merit   A_lambda')

do i = 1, s%global%n_opti_cycles+1

  if (a_lambda > 1e10) then
    call out_io (s_blank$, r_name, &
                    'Optimizer at minimum or derivatives need to be recalculated.')
    finished = .true.
  endif

  if (finished .or. i == s%global%n_opti_cycles+1) then
    a_lambda = 0  ! tell mrqmin we are finished
    call tao_var_write (s%global%var_out_file)
  endif

  call tao_mrqmin (y, weight, a, tao_com%covar, tao_com%alpha, chi_sq, a_lambda, limited) 
  call tao_mrq_func (a, y_fit, dy_da, limited2)  ! put a -> model
  write (line, '(i5, es14.4, es10.2)') i, tao_merit(), a_lambda
  call out_io (s_blank$, r_name, line)

  if (finished .or. limited) return

! look for keyboard input to end optimization

#ifndef CESR_WINCVF
  do
    call get_tty_char (char, .false., .false.) 
    if (char == '.') then
      tao_com%optimizer_running = .false.
      call out_io (s_blank$, r_name, 'Optimizer stop signal detected.', 'Stopping now.')
      abort = .true.
      finished = .true.
      exit
    endif
    if (char == achar(0)) exit   ! only exit if there is no more input
  enddo
#endif

! reinit the derivative matrix 

  if (s%global%lm_opt_deriv_reinit > 0 .and. a_lambda > s%global%lm_opt_deriv_reinit) then
    call tao_dmodel_dvar_calc (.true.)
    a_lambda = 1
  endif

enddo

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_mrq_func (a, y_fit, dy_da, limited)
! 
! Subroutine to be called by the Numerical Recipes routine mrqmin

subroutine tao_mrq_func (a, y_fit, dy_da, limited)

implicit none

type (tao_universe_struct), pointer :: u

real(rp), intent(in) :: a(:)
real(rp), intent(out) :: y_fit(:)
real(rp), intent(out) :: dy_da(:, :)
real(rp) merit0
real(rp), allocatable, save :: var_delta(:)

integer i, j, k, n, nn, im, iv, n_var

logical limited
character(80) line

! transfer "a" array to model

call tao_set_vars (a)

! if limited then set y_fit to something large so merit calc gives a large number.

call tao_limit_calc (limited)

if (limited) then
  y_fit = 1e10 
  return
endif

! calculate derivatives

merit0 = tao_merit()

dy_da = 0
n_var = size(a)

call tao_get_vars (var_delta = var_delta)
y_fit(1:n_var) = var_delta

forall (k = 1:n_var) dy_da(k,k) = 1

k = n_var

do j = lbound(s%u, 1), ubound(s%u, 1)
  u => s%u(j)
  if (.not. u%is_on) cycle
  do i = 1, size(u%data)
    if (.not. u%data(i)%useit_opt) cycle
    if (u%data(i)%weight == 0) cycle
    k = k + 1
    y_fit(k) = u%data(i)%delta_merit
    im = u%data(i)%ix_dModel
    nn = 0
    do n = 1, size(s%var)
      if (.not. s%var(n)%useit_opt) cycle
      nn = nn + 1
      iv = s%var(n)%ix_dVar
      dy_da(k, nn) = u%dModel_dVar(im, iv)
    enddo
  enddo
enddo

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

subroutine tao_mrqmin(y,weight,a,covar,alpha,chisq,alamda, limited)

implicit none

real(rp) :: y(:), weight(:)
real(rp) :: a(:)
real(rp) :: covar(:,:), alpha(:,:)
real(rp) :: chisq
real(rp) :: alamda
integer(i4b) :: ma,ndata
integer(i4b), save :: mfit
logical, allocatable, save :: mask(:)
logical limited

real(rp), save :: ochisq
real(rp), dimension(:), allocatable, save :: atry,beta
real(rp), dimension(:,:), allocatable, save :: da

!

ndata=assert_eq(size(y),size(weight),'mrqmin: ndata')
ma=assert_eq((/size(a),size(covar,1),size(covar,2), &
              size(alpha,1),size(alpha,2)/),'mrqmin: ma')
mfit=size(a)

if (allocated (mask)) then
  if (size(mask) /= size(a)) deallocate (mask)
endif
if (.not. allocated(mask)) then
  allocate (mask(size(a)))
  mask = .true.
endif

if (alamda < 0.0) then
  allocate(atry(ma),beta(ma),da(ma,1))
  alamda=0.001_rp
  call tao_mrqcof(a, y, alpha, beta, weight, chisq, limited)
  if (limited) then
    deallocate(atry,beta,da)
    return
  endif
  ochisq=chisq
  atry=a
end if

covar(1:mfit,1:mfit)=alpha(1:mfit,1:mfit)
call diagmult(covar(1:mfit,1:mfit),1.0_rp+alamda)
da(1:mfit,1)=beta(1:mfit)
call gaussj(covar(1:mfit,1:mfit),da(1:mfit,1:1))

if (alamda == 0.0) then
  call covsrt(covar,mask)
  call covsrt(alpha,mask)
  deallocate(atry,beta,da)
  return
end if

atry=a+unpack(da(1:mfit,1),mask,0.0_rp)
call tao_mrqcof(atry, y, covar, da(1:mfit,1), weight, chisq, limited)
if (limited) return

if (chisq < ochisq) then
  alamda=0.1_rp*alamda
  ochisq=chisq
  alpha(1:mfit,1:mfit)=covar(1:mfit,1:mfit)
  beta(1:mfit)=da(1:mfit,1)
  a=atry
else
  alamda=10.0_rp*alamda
  chisq=ochisq
end if

end subroutine tao_mrqmin

!-----------------------------------------------------------

subroutine tao_mrqcof(a, y, alpha, beta, weight, chisq, limited)

implicit none

real(rp) :: y(:), a(:), weight(:)
real(rp) :: beta(:)
real(rp) :: alpha(:,:)
real(rp) chisq
integer(i4b) :: j,k,l,m, nv, nd
real(rp), allocatable, save :: dyda(:,:)
real(rp), allocatable, save :: dy(:),wt(:),ymod(:)
logical limited

!

nd = size(weight)
nv = size(a)

if (allocated(dyda)) then
  if (size(dyda, 1) /= nd .or. size(dyda, 2) /= nv) &
                                        deallocate (dyda, dy, wt, ymod)
endif
if (.not. allocated(dyda)) then
  allocate (dyda(nd,nv), dy(nd),wt(nd),ymod(nd))
endif

!

call tao_mrq_func(a, ymod, dyda, limited)
if (limited) return

dy=y-ymod
j=0

do l=1,nv
  j=j+1
  wt=dyda(:,l) * weight
  k=0
  do m=1,l
    k=k+1
    alpha(j,k)=dot_product(wt,dyda(:,m))
    alpha(k,j)=alpha(j,k)
  end do
  beta(j)=dot_product(dy,wt)
end do

chisq = dot_product(dy**2,weight)

end subroutine tao_mrqcof

end module
