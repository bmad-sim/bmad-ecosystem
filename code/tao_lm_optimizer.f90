!+
! Subroutine tao_lm_optimizer ()
!
! Subrutine to minimize the merit function by varying variables until
! the "data" as calculated from the model matches the measured data.
! 
! This subroutine is a wrapper for the mrqmin routine of Numerical Recipes.
! See the Numerical Recipes writeup for more details.
! 'lm' stands for Levenburg - Marquardt. Otherwise known as LMDIF. 
!
! Note: LM assumes 
!
! Input:
!
! Output:
!-

subroutine tao_lm_optimizer ()

use tao_mod
use nr
use tao_dmerit_mod
use single_char_input_mod

implicit none

type (tao_universe_struct), pointer :: u

real(rp), allocatable, save :: x(:), y(:), sig(:), a(:)
real(rp), allocatable, save :: covar(:,:), alpha(:,:)
real(rp), allocatable, save :: y_fit(:)
real(rp), allocatable, save :: dy_da(:, :)
real(rp), allocatable, save :: var_value(:), weight(:), var_meas_value(:)
real(rp) a_lambda, chi_sq, merit0

integer i, j, k
integer n_data, n_var

logical, allocatable, save :: mask_a(:)
logical :: finished, init_needed = .true.

character(20) :: r_name = 'tao_lm_optimizer'
character(80) line
character(1) char

! setup

a_lambda = -1

call tao_get_vars (var_value, var_weight = weight, var_meas_value = var_meas_value)
n_var = size(var_value)

n_data = n_var
do i = 1, size(s%u)
  n_data = n_data + count(s%u(i)%data(:)%useit_opt)
enddo

if (allocated(x)) deallocate(x, y, sig, a, covar, alpha, mask_a, y_fit, dy_da)
allocate (x(n_data), y(n_data), sig(n_data), y_fit(n_data))
allocate (a(n_var), covar(n_var,n_var), alpha(n_var,n_var), mask_a(n_var))
allocate (dy_da(n_data, n_var))
mask_a = .true.

! init a and y arrays

a(1:n_var) = var_value
y(1:n_var) = var_meas_value
sig(1:n_var) = 1e10  ! something large
where (weight /= 0) sig(1:n_var) = sqrt(1/weight)

merit0 = tao_merit()

k = n_var
do j = 1, size(s%u)
  u => s%u(j)
  do i = 1, size(u%data)
    if (.not. u%data(i)%useit_opt) cycle
    k = k + 1
    y(k) = 0
    if (u%data(i)%weight == 0) then
      sig(k) = 1e10
    else
      sig(k) = sqrt(1/u%data(i)%weight)
    endif
  enddo
enddo

call tao_dModel_dVar_calc (s%global%derivative_recalc)

! run optimizer mrqmin from Numerical Recipes.

call out_io (s_blank$, r_name, '   Loop      Merit   A_lambda')

do i = 1, s%global%n_opti_cycles
  if (i == s%global%n_opti_cycles) a_lambda = 0  ! tell mrqmin we are finished
  call mrqmin (x, y, sig, a, mask_a, covar, alpha, chi_sq, tao_mrq_func, a_lambda) 
  call tao_mrq_func (x, a, y_fit, dy_da)  ! put a -> model
  write (line, '(i5, es14.4, es10.2)'), i, tao_merit(), a_lambda
  call out_io (s_blank$, r_name, line)

! look for keyboard input to end optimization

  do
    call get_tty_char (char, .false., .false.) 
    if (char == '.') then
      s%global%optimizer_running = .false.
      call out_io (s_blank$, r_name, 'Optimizer stop signal detected.', &
                                                             'Stopping now.')
      return
    endif
    if (char == achar(0)) exit   ! only exit if there is no more input
  enddo

enddo

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_mrq_func (x, a, y_fit, dy_da)
! 
! Subroutine to be called by the Numerical Recipes routine mrqmin


subroutine tao_mrq_func (x, a, y_fit, dy_da)

use tao_mod

implicit none

type (tao_universe_struct), pointer :: u

real(rp), intent(in) :: x(:), a(:)
real(rp), intent(out) :: y_fit(:)
real(rp), intent(out) :: dy_da(:, :)
real(rp) merit0

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
forall (k = 1:n_var) dy_da(k,k) = 1

y_fit(1:n_var) = a

k = n_var

do j = 1, size(s%u)
  u => s%u(j)
  do i = 1, size(u%data)
    if (.not. u%data(i)%useit_opt) cycle
    k = k + 1
    y_fit(k) = u%data(i)%delta
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

