module tao_svd_optimizer_mod

use tao_mod
use tao_dmerit_mod
use tao_var_mod
use tao_lm_optimizer_mod

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_svd_optimizer 
!
! Routine to minimize the merit function using svd.
!-

subroutine tao_svd_optimizer ()

use nr

implicit none

type (tao_universe_struct), pointer :: u

real(rp), allocatable, save :: weight(:), a(:), a_try(:), da(:), b(:), w(:)
real(rp), allocatable, save :: y_fit(:)
real(rp), allocatable, save :: dy_da(:, :), v(:, :)
real(rp), allocatable, save :: var_value(:), var_weight(:)
real(rp) merit0, merit
real(rp), parameter :: tol = 1d-5

integer i, j, k, status, status2
integer n_data, n_var
integer, allocatable, save :: var_ix(:)

logical :: finished, init_needed = .true.

character(20) :: r_name = 'tao_svd_optimizer'
character(80) line
character(1) char

! Calc derivative matrix

call tao_dModel_dVar_calc (s%global%derivative_recalc, .true.)

! setup

merit0 = tao_merit()

call tao_get_vars (var_value, var_weight = var_weight, var_ix = var_ix)
n_var = size(var_value)

n_data = n_var
do i = lbound(s%u, 1), ubound(s%u, 1)
  if (.not. s%u(i)%is_on) cycle
  n_data = n_data + count(s%u(i)%data(:)%useit_opt .and. s%u(i)%data(:)%weight /= 0)
enddo

if (allocated(weight)) deallocate(weight, a, a_try, da, y_fit, dy_da, b, w, v)
allocate (weight(n_data), y_fit(n_data))
allocate (a(n_var), a_try(n_var), da(n_var), b(n_var), w(n_var))
allocate (dy_da(n_data, n_var), v(n_var, n_var))

! init a and y arrays

a = var_value
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

!--------------------------
! SVD 

merit0 = tao_merit()
call out_io (s_blank$, r_name, 'Initial Merit:   \es14.4\ ', r_array = [merit0])

call tao_mrq_func (a, y_fit, dy_da, status)  ! put a -> model
if (status /= 0) return

b = y_fit * sqrt(weight)
do j = 1, n_data
  dy_da(j, :) = dy_da(j, :) / sqrt(weight(j))
enddo
call svdcmp(dy_da, w, v)
where (w < tol * maxval(w)) w = 0
call svbksb (dy_da, w, v, b, da)
a_try = a + da

call tao_mrq_func (a_try, y_fit, dy_da, status)  ! put a -> model
merit = tao_merit()
call out_io (s_blank$, r_name, 'Merit after svd: \es14.4\ ', r_array = [merit])

if (status /= 0 .or. merit > merit0) then
  call tao_mrq_func (a, y_fit, dy_da, status2)  ! put a -> model
  call out_io (s_blank$, r_name, 'Retreating to initial state.')
endif

end subroutine

end module
