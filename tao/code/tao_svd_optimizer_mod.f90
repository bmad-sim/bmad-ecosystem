module tao_svd_optimizer_mod

use tao_interface
use tao_dmerit_mod

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_svd_optimizer (abort)
!
! Routine to minimize the merit function using svd.
!
! Output:
!   abort -- Logical: Set True if svd step increases the merit function.
!-

subroutine tao_svd_optimizer (abort)

implicit none

type (tao_universe_struct), pointer :: u

real(rp), allocatable :: weight(:), a(:), a_try(:), da(:), b(:)
real(rp), allocatable :: y_fit(:)
real(rp), allocatable :: dy_da(:, :), v(:, :)
real(rp), allocatable :: var_value(:), var_weight(:)
real(rp) merit0, merit

integer i, j, k, status, status2
integer n_data, n_var

logical abort

character(20) :: r_name = 'tao_svd_optimizer'
character(80) line
character(1) char

! Calc derivative matrix

call tao_dModel_dVar_calc (s%global%derivative_recalc, abort)
if (abort) return
call tao_veto_vars_with_zero_dmodel ()

! Setup

merit0 = tao_merit()

call tao_get_opt_vars (var_value)
n_var = size(var_value)

call tao_get_opt_vars (var_weight = var_weight, ignore_if_weight_is_zero = .true., &
                                                             ignore_if_not_limited = .true.)

n_data = size(var_weight)
do i = lbound(s%u, 1), ubound(s%u, 1)
  if (.not. s%u(i)%is_on) cycle
  n_data = n_data + count(s%u(i)%data(:)%useit_opt .and. s%u(i)%data(:)%weight /= 0 .and. &
                                                                     s%u(i)%data(:)%good_model)
enddo

allocate (weight(n_data), y_fit(n_data), b(n_data))
allocate (a(n_var), a_try(n_var), da(n_var))
allocate (dy_da(n_data, n_var), v(n_var, n_var))

! init a and y arrays

a = var_value
k = size(var_weight)
weight(1:k) = var_weight

do j = lbound(s%u, 1), ubound(s%u, 1)
  u => s%u(j)
  if (.not. u%is_on) cycle
  do i = 1, size(u%data)
    if (.not. u%data(i)%useit_opt) cycle
    if (.not. u%data(i)%good_model) cycle
    if (u%data(i)%weight == 0) cycle
    k = k + 1
    weight(k) = u%data(i)%weight
  enddo
enddo

!--------------------------
! SVD 

merit0 = tao_merit()
call out_io (s_blank$, r_name, 'Initial Merit: \es14.4\ ', r_array = [merit0])

call tao_svd_func (a, y_fit, dy_da, status)  ! put a -> model
if (status /= 0) return

do i = 1, n_data
  dy_da(i, :) = dy_da(i, :) * sqrt(weight(i))
enddo

b = y_fit * sqrt(weight)

call svd_fit (dy_da, b, s%global%svd_cutoff, da)
a_try = a - da

!

call tao_set_opt_vars (a_try, s%global%optimizer_var_limit_warn)
merit = tao_merit()
call out_io (s_blank$, r_name, 'Merit after svd: \es14.4\ ', r_array = [merit])

if (status /= 0 .or. (merit > merit0 .and. s%global%svd_retreat_on_merit_increase)) then
  abort = .true.   ! So no more loops
  call tao_set_opt_vars (a, s%global%optimizer_var_limit_warn)
  call out_io (s_blank$, r_name, 'Retreating to initial state.')
endif

if (.not. abort) abort = tao_user_is_terminating_optimization()

! Cleanup: Reinstate vars vetoed with zero dmerit

s%var(:)%good_var = .true.  
call tao_set_var_useit_opt()

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_svd_func (a, y_fit, dy_da, status)
!-

subroutine tao_svd_func (a, y_fit, dy_da, status)

implicit none

type (tao_universe_struct), pointer :: u

real(rp), intent(in) :: a(:)
real(rp), intent(out) :: y_fit(:)
real(rp), intent(out) :: dy_da(:, :)
real(rp) merit0
real(rp), allocatable :: var_delta(:)

integer i, j, k, n, nn, im, iv, n_var, nd
integer status

logical limited

character(80) line

! transfer "a" array to model

call tao_set_opt_vars (a, s%global%optimizer_var_limit_warn)

! if limited then set y_fit to something large so merit calc gives a large number.

call tao_limit_calc (limited)

if (limited) then
  status = not_set$
  y_fit = 1e10 
  return
endif

! Calculate derivatives

merit0 = tao_merit()  ! Calculate %delta_merit values

dy_da = 0
n_var = size(a)

call tao_get_opt_vars (var_delta = var_delta, ignore_if_weight_is_zero = .true., &
                                              ignore_if_not_limited = .true.)
nd = size(var_delta)
y_fit(1:nd) = var_delta

forall (k = 1:nd) dy_da(k,k) = 1

do j = lbound(s%u, 1), ubound(s%u, 1)
  u => s%u(j)
  if (.not. u%is_on) cycle
  do i = 1, size(u%data)
    if (.not. u%data(i)%useit_opt) cycle
    if (.not. u%data(i)%good_model) cycle
    if (u%data(i)%weight == 0) cycle
    nd = nd + 1
    y_fit(nd) = u%data(i)%delta_merit
    im = u%data(i)%ix_dModel
    nn = 0
    do n = 1, s%n_var_used
      if (.not. s%var(n)%useit_opt) cycle
      nn = nn + 1
      iv = s%var(n)%ix_dVar
      dy_da(nd, nn) = u%dModel_dVar(im, iv)
    enddo
  enddo
enddo

status = 0

end subroutine

end module
