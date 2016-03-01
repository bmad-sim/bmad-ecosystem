!+
! Program wake_fit
!
! Program to fit a wakefield to a set of Bmad pseudo wake modes.
!-

program wake_fit

use sim_utils
use super_recipes_mod

implicit none

type wake_mode_struct
  real(rp) amp    ! Amplitude
  real(rp) damp   ! Damping factor
  real(rp) k      ! wave number
end type

type (wake_mode_struct) mode(100)

real(rp) amp_weight, chisq, a_lambda, merit_min, merit_now, negligible_merit
real(rp), allocatable :: y(:),  weight(:), a(:)
real(rp), allocatable :: y_fit(:)
real(rp), allocatable :: dy_da(:,:), alpha(:,:), covar(:,:)

integer i, j, n_mode, n_var, n_data, n_cycles, n_pts, ix_data_start, ix_min
integer status, n_y, id, iv, iv0

character(100) param_file

namelist / params / mode, n_mode, n_data, ix_data_start, amp_weight, n_cycles, negligible_merit

! Read in parameters

param_file = 'wake_fit.init'

open (1, file = param_file)
read (1, nml = params)
close (1)

! Var_array = amp_i, damp_i, k_i,  i = 1, ..., n_mode
! y_array = amp_i, data

! Fit

n_var = 3 * n_mode
n_y = n_mode + n_data

allocate (weight(n_y), y(n_y), y_fit(n_y), a(n_var))
allocate (covar(n_var, n_var), alpha(n_var, n_var))
allocate (dy_da(n_y, n_var))

weight = 1
weight(1:n_mode) = amp_weight

y = 0
do id = 1, n_data
  y(id+n_mode) = 1 / sqrt(real(ix_data_start + id - 1, rp))
enddo

do iv = 1, n_mode
  iv0 = 3 * (iv - 1)
  a(iv0+1) = mode(iv)%amp
  a(iv0+2) = mode(iv)%damp
  a(iv0+3) = mode(iv)%k
enddo

!

merit_min = merit()
ix_min = 0

do i = 1, 1000000
  call super_mrqmin (y, weight, a, covar, alpha, chisq, mrq_func, a_lambda, status)
  call mrq_func (a, y_fit, dy_da, status)

  merit_now = merit()
  print '(i5, es14.4, es10.2, 10x, 100f10.4)', i, merit_now, a_lambda, (y(j)-y_fit(j), j = 1, n_y)

  if (merit_now < 0.99 * merit_min) then
    merit_min = merit_now
    ix_min = i
  endif

  if (ix_min - i > n_cycles .or. a_lambda > 1d10) then
    print *, 'At local minimum. Stopping.'
    exit
  endif

  if (merit_now < negligible_merit) then
    print *, 'Success! Merit below negligible_merit.'
    exit
  endif

enddo

! Write parameters

print *
do i = 1, n_mode
  print '(i5, 3f10.5)', i, mode(i)%amp, mode(i)%damp, mode(i)%k
enddo

!----------------------------------------------------------------------------------
contains

subroutine mrq_func (a, y_fit, dy_da, status)

real(rp), intent(in) :: a(:)
real(rp), intent(out) :: y_fit(:)
real(rp), intent(out) :: dy_da(:, :)
real(rp) damp, cosine

integer status
integer j, iv, id, iid, iy, iv0

!

status = 0
dy_da = 0
y_fit = 0

do iv = 1, n_mode
  iv0 = 3*(iv-1)
  mode(iv) = wake_mode_struct(a(iv0+1), a(iv0+2), a(iv0+3))
  y_fit(iv) = a(iv0+1)
  dy_da(iv,iv0+1) = 1
enddo

do id = 1, n_data
  iy = id + n_mode
  iid = id + ix_data_start - 1
  do iv = 1, n_mode
    iv0 = 3*(iv-1)
    damp = exp(-iid * a(iv0+2)) 
    cosine = cos(iid * a(iv0+3))
    y_fit(iy) = y_fit(iy) + a(iv0+1) * damp * cosine
    dy_da(iy, iv0+1) = damp * cosine
    dy_da(iy, iv0+2) = -iid * a(iv0+1) * damp * cosine
    dy_da(iy, iv0+3) = -a(iv0+1) * damp * sin(iid * a(iv0+3))
  enddo
enddo

end subroutine mrq_func

!----------------------------------------------------------------------------------
! contains

function merit() result (this_merit)

real(rp) this_merit
integer iv, id, iid
real(rp) y_fit, damp, c

!

this_merit = 0
do iv = 1, n_mode
  this_merit = this_merit + amp_weight * mode(iv)%amp**2
enddo

do id = 1, n_data
  iid = id + ix_data_start - 1
  y_fit = 0
  do iv = 1, n_mode
    damp = exp(-iid * mode(iv)%damp) 
    c = cos(iid * mode(iv)%k)
    y_fit = y_fit + mode(iv)%amp * damp * c
  enddo
  this_merit = this_merit + (1 / sqrt(real(iid, rp)) - y_fit)**2
enddo

end function merit

end program
