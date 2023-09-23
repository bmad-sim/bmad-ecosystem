module tao_lm_optimizer_mod

use tao_interface
use tao_dmerit_mod
use tao_top10_mod
use input_mod
use super_recipes_mod

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_lm_optimizer (abort)
!
! Routine to minimize the merit function by varying variables until
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
type (super_mrqmin_storage_struct) storage

real(rp), allocatable :: y(:), weight(:), a(:)
real(rp), allocatable :: y_fit(:)
real(rp), allocatable :: dy_da(:, :)
real(rp), allocatable :: var_value(:), var_weight(:)
real(rp) a_lambda, chi_sq, merit0

integer i, j, k, status, status2
integer n_data, n_var
integer, allocatable :: var_ix(:)

logical :: finished, init_needed = .true.
logical abort

character(20) :: r_name = 'tao_lm_optimizer'
character(80) line

! Calc derivative matrix

call tao_dModel_dVar_calc (s%global%derivative_recalc, abort)
if (abort) return
call tao_veto_vars_with_zero_dmodel ()

! setup

a_lambda = -1

merit0 = tao_merit()

call tao_get_opt_vars (var_value, var_ix = var_ix)
n_var = size(var_value)

call tao_get_opt_vars (var_weight = var_weight, ignore_if_weight_is_zero = .true.)

n_data = size(var_weight)
do i = lbound(s%u, 1), ubound(s%u, 1)
  if (.not. s%u(i)%is_on) cycle
  n_data = n_data + count(s%u(i)%data(:)%useit_opt .and. s%u(i)%data(:)%weight /= 0)
enddo

if (allocated(y)) deallocate(y, weight, a, y_fit, dy_da)
if (allocated(s%com%covar)) deallocate (s%com%covar, s%com%alpha)
allocate (y(n_data), weight(n_data), y_fit(n_data))
allocate (a(n_var), s%com%covar(n_var,n_var), s%com%alpha(n_var,n_var))
allocate (dy_da(n_data, n_var))

! init a and y arrays

a = var_value
y = 0

k = size(var_weight)
weight(1:k) = var_weight
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

! run optimizer mrqmin from Numerical Recipes.

finished = .false.
call out_io (s_blank$, r_name, '  Cycle      Merit    A_lambda')

do i = 1, s%global%n_opti_cycles+1

  if (a_lambda > 1e10) then
    call out_io (s_blank$, r_name, 'Optimizer at minimum or derivatives need to be recalculated.')
    finished = .true.
  endif

  ! Underflow can cause lm to bomb so limit min of a_lambda
  if (a_lambda < 1d-98 .and. a_lambda > 0) a_lambda = 1d-98

  if (finished .or. i == s%global%n_opti_cycles+1) then
    a_lambda = 0  ! tell mrqmin we are finished
    if (s%opti_write_var_file) call tao_var_write (s%global%var_out_file)
  endif

  call super_mrqmin (y, weight, a, chi_sq, tao_mrq_func, storage, a_lambda, status)
  if (status == 0) then  ! If everything OK
    s%com%covar = storage%covar
    s%com%alpha = storage%alpha
  endif

  call tao_mrq_func (a, y_fit, dy_da, status2)  ! put a -> model
  write (line, '(i5, es14.4, es11.2)') i, tao_merit(), a_lambda
  call out_io (s_blank$, r_name, line)

  ! Try to find problem variable!

  if (status == -1 .or. status == -2) then  ! gaussj singular error
    do k = 1, n_var
      if (all (s%com%covar(k,:) == 0) .or. all(s%com%covar(:,k) == 0)) then
        call out_io (s_error$, r_name, 'Problem variable: ' // tao_var1_name(s%var(var_ix(k))))
      endif
    enddo
    abort = .true.
  endif

  ! look for keyboard input to end optimization

  if (tao_user_is_terminating_optimization()) then
    abort = .true.
    finished = .true.
  endif

  if (finished .or. status /= 0) exit

  ! reinit the derivative matrix 

  if (s%global%lm_opt_deriv_reinit > 0 .and. a_lambda > s%global%lm_opt_deriv_reinit) then
    call tao_dmodel_dvar_calc (.true., abort)
    a_lambda = 1
  endif

enddo

! Cleanup: Reinstate vars vetoed with zero dmerit

s%var(:)%good_var = .true.  
call tao_set_var_useit_opt()

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_mrq_func (a, y_fit, dy_da, status)
! 
! Subroutine to be called by the Numerical Recipes routine mrqmin
!-

subroutine tao_mrq_func (a, y_fit, dy_da, status)

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

call tao_get_opt_vars (var_delta = var_delta, ignore_if_weight_is_zero = .true.)
nd = size(var_delta)
y_fit(1:nd) = var_delta

forall (k = 1:nd) dy_da(k,k) = 1

do j = lbound(s%u, 1), ubound(s%u, 1)
  u => s%u(j)
  if (.not. u%is_on) cycle
  do i = 1, size(u%data)
    if (.not. u%data(i)%useit_opt) cycle
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
