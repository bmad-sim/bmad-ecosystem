module tao_geodesic_lm_optimizer_mod

use tao_interface
use tao_dmerit_mod
use tao_top10_mod
use input_mod
use geodesic_lm

real(rp), allocatable :: dy_da(:, :)
logical, save :: geolevmar_limit_flag = .false.

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_geodesic_lm_optimizer (abort)
!
! Routine to minimize the merit function by varying variables until
! the "data" as calculated from the model matches the measured data.
! 
! This subroutine is a wrapper for the "geodesic"
! Levenburg - Marquardt method.
!
! Output:
!   abort -- Logical: Set True if an user stop signal detected.
!-

subroutine tao_geodesic_lm_optimizer (abort)

implicit none

type (tao_universe_struct), pointer :: u

real(rp), allocatable :: y(:), a(:)
real(rp), allocatable :: y_fit(:)
real(rp), allocatable :: var_value(:), var_weight(:)
real(rp) a_lambda, chi_sq, merit0

real(rp), allocatable :: fjac(:,:), dtd(:,:)

integer i, j, k
integer n_data, n_var, info, niters, nfev, njev, naev, converged
integer, allocatable :: var_ix(:)

logical :: finished, init_needed = .true.
logical abort

character(*), parameter :: r_name = 'tao_geodesic_lm_optimizer'
character(80) line

!

call tao_get_opt_vars (var_value, var_ix = var_ix)
n_var = size(var_value)

call tao_get_opt_vars (var_weight = var_weight, ignore_if_weight_is_zero = .true.)

n_data = size(var_weight)
do i = lbound(s%u, 1), ubound(s%u, 1)
  if (.not. s%u(i)%is_on) cycle
  n_data = n_data + count(s%u(i)%data(:)%useit_opt .and. s%u(i)%data(:)%weight /= 0)
enddo

if (allocated(dy_da)) deallocate(dy_da)
if (allocated(s%com%covar)) deallocate (s%com%covar, s%com%alpha)
allocate (y(n_data), y_fit(n_data))
allocate (a(n_var), s%com%covar(n_var,n_var), s%com%alpha(n_var,n_var))
allocate (dy_da(n_data, n_var))

allocate (fjac(n_data, n_var), dtd(n_var, n_var))

! init a and y arrays

a = var_value
y = 0

! run geodesic_lm

info = 0
niters = 0
geolevmar_limit_flag = .false.

call run_geodesic_lm (tao_geo_lm_func, jacobian, Avv, a, y_fit, fjac, callback, info, &
                      dtd,  niters, nfev, njev, naev, converged)
call tao_geo_lm_func (n_data, n_var, a, y_fit)  ! put a -> model

! look for keyboard input to end optimization

if (converged /= 0) abort = .true.

if (tao_user_is_terminating_optimization()) then
  abort = .true.
  finished = .true.
endif

! Cleanup: Reinstate vars vetoed with zero dmerit

s%var(:)%good_var = .true.  
call tao_set_var_useit_opt()

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_geo_lm_func (mm, nn, a, y_fit)
! 
! Subroutine to be called by the geodesic_lm code.
!-

subroutine tao_geo_lm_func (mm, nn, a, y_fit)

implicit none

type (tao_universe_struct), pointer :: u

integer mm, nn

real(rp) :: a(nn)
real(rp) :: y_fit(mm)
real(rp) merit0
real(rp), allocatable :: var_delta(:)

integer i, j, k, n, nv, im, iv, n_var, nd

logical limited

character(80) line

! transfer "a" array to model

call tao_set_opt_vars (a, s%global%optimizer_var_limit_warn)

! if limited then set y_fit to something large so merit calc gives a large number.

call tao_limit_calc (limited)

if (limited) then
  y_fit = 1e10 
  return
endif

! Calculate derivatives

merit0 = tao_merit()  ! Calculate %delta_merit values

n_var = size(a)

call tao_get_opt_vars (var_delta = var_delta, ignore_if_weight_is_zero = .true.)
nd = size(var_delta)
y_fit(1:nd) = var_delta

do j = lbound(s%u, 1), ubound(s%u, 1)
  u => s%u(j)
  if (.not. u%is_on) cycle
  do i = 1, size(u%data)
    if (.not. u%data(i)%useit_opt) cycle
    if (u%data(i)%weight == 0) cycle
    nd = nd + 1
    y_fit(nd) = u%data(i)%delta_merit * sqrt(u%data(i)%weight)
  enddo
enddo

end subroutine tao_geo_lm_func

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

subroutine callback (m, n, x, fvec, fjac, accepted, info)

implicit none

integer m, n, info, accepted
real(8) x(n), fvec(m), fjac(m, n)

if (geolevmar_limit_flag .and. accepted > 0) then
  info = 1      ! ends geolevmar, returns control
  write(*, *) "limit hit"
endif

end subroutine callback

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

subroutine jacobian(mm, nn, x, fjac)
        
implicit none

type (tao_universe_struct), pointer :: u

integer :: mm, nn, nd, k, i, j, im, iv, nv, n
real(8) :: x(nn), fjac(mm, nn)
real(rp), allocatable :: var_delta(:)
logical err_flag

!

call tao_dModel_dVar_calc (.true., err_flag)

call tao_get_opt_vars (var_delta = var_delta, ignore_if_weight_is_zero = .true.)
nd = size(var_delta)

dy_da = 0
forall (k = 1:nd) dy_da(k,k) = 1

do j = lbound(s%u, 1), ubound(s%u, 1)
  u => s%u(j)
  if (.not. u%is_on) cycle
  do i = 1, size(u%data)
    if (.not. u%data(i)%useit_opt) cycle
    if (u%data(i)%weight == 0) cycle
    nd = nd + 1

    im = u%data(i)%ix_dModel
    nv = 0
    do n = 1, s%n_var_used
      if (.not. s%var(n)%useit_opt) cycle
      nv = nv + 1
      iv = s%var(n)%ix_dVar
      dy_da(nd, nv) = u%dModel_dVar(im, iv) * sqrt(u%data(i)%weight)
    enddo
  enddo
enddo

fjac = dy_da

end subroutine jacobian

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

subroutine avv(m, n, x, v, acc)
        
implicit none

integer :: m, n
real(8) :: x(n), v(n), acc(n)

end subroutine

end module

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!
! Signal for early termination
!
subroutine user_signal(flag, is_set)

use tao_interface

implicit none

integer :: flag
logical :: is_set

if (tao_user_is_terminating_optimization()) then
  flag = -10
  is_set = .true.
else
  is_set = .false.
endif

end subroutine
