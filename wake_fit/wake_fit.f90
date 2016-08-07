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
  real(rp) :: amp = 0    ! Amplitude
  real(rp) :: damp = 0   ! Damping factor
  real(rp) :: k = 0      ! wave number
  real(rp) :: phi = 0    ! phase
end type

type (wake_mode_struct) mode(100)
type (super_mrqmin_storage_struct) storage

real(rp) amp_weight, chisq, a_lambda, merit_min, merit_now, negligible_merit
real(rp), allocatable :: y(:),  weight(:), a(:)
real(rp), allocatable :: y_fit(:)
real(rp), allocatable :: dy_da(:,:)

integer i, j, k, n, n_var, n_data, n_cycles, n_pts, ix_data_start, ix_min
integer status, n_y, id, iv, iv0, n_mode, n_var_in_mode

logical use_phase

character(100) param_file

namelist / params / mode, n_data, ix_data_start, amp_weight, n_cycles, negligible_merit, use_phase

! Read in parameters

param_file = 'wake_fit.init'

open (1, file = param_file)
read (1, nml = params)
close (1)

! Var_array = amp_i, damp_i, k_i,  i = 1, ..., n_mode
! y_array = amp_i, data

! Fit

do iv = 1, size(mode)
  if (mode(iv)%amp /= 0 .or. mode(iv)%damp /= 0 .or. mode(iv)%k /= 0 .or. mode(iv)%phi /= 0) cycle
  n_mode = iv - 1
  exit
end do

print *, 'Number of modes:', n_mode

if (use_phase) then
  n_var_in_mode = 4
else
  n_var_in_mode = 3
  mode%phi = 0
endif

n_var = n_var_in_mode * n_mode
n_y = n_mode + n_data

allocate (weight(n_y), y(n_y), y_fit(n_y), a(n_var))
allocate (dy_da(n_y, n_var))

weight = 1
weight(1:n_mode) = amp_weight

y = 0
do id = 1, n_data
  y(id+n_mode) = 1 / sqrt(real(ix_data_start + id - 1, rp))
enddo

do iv = 1, n_mode
  iv0 = n_var_in_mode * (iv - 1)
  a(iv0+1) = mode(iv)%amp
  a(iv0+2) = mode(iv)%damp
  a(iv0+3) = mode(iv)%k
  if (use_phase) a(iv0+4) = mode(iv)%phi
enddo

!

merit_min = merit()
ix_min = 0
a_lambda = -1

do i = 1, 1000000

  if (modulo(i, 100) == 0) then
    open (1, file = 'wake.fit')
    print *

    do iv = 1, n_mode
      print '(a, i0, a, 4f18.11)', '  mode(', iv, ') =', mode(iv)%amp, mode(iv)%damp, mode(iv)%k, mode(iv)%phi
      write (1, '(a, i0, a, 4f18.11)') '  mode(', iv, ') =', mode(iv)%amp, mode(iv)%damp, mode(iv)%k, mode(iv)%phi
    enddo

    do k = 1, nint(log10(n_data-n_mode-1.0_rp)) - 1
      j = nint(10.0**k)
      n = n_mode + 1
      print '(2i10, es10.3)', k, j, sum(abs(y(n+j:n+10*j)-y_fit(n+j:n+10*j)))/sum(y(n+j:n+10*j))
      write (1, '(2i10, es10.3)') k, j, sum(abs(y(n+j:n+10*j)-y_fit(n+j:n+10*j)))/sum(y(n+j:n+10*j))
    enddo

    close (1)
  endif

  call super_mrqmin (y, weight, a, chisq, mrq_func, storage, a_lambda, status)
  if (status /= 0) then
    a_lambda = 2 * a_lambda
    cycle
  endif
  call mrq_func (a, y_fit, dy_da, status)

  ! To look nice, restrict phase to be in the range [-pi/2, pi/2]
  if (use_phase) then
    do iv = 1, n_mode
      iv0 = 4 * (iv - 1)
      a(iv0+4) = modulo2 (a(iv0+4), pi)  ! phi
      if (abs(a(iv0+4)) > pi/2) then
        a(iv0+4) = modulo2 (a(iv0+4), pi/2)  ! phi
        a(iv0+1) = -a(iv0+1)
      endif
    enddo
  endif

  merit_now = merit()
  n = n_mode
  print '(i5, 2es14.4, es10.2, 10x, 10f10.4, 2x, a, 3f10.4)', i, sqrt(merit_now)/n_data, sum(abs(y(n+1:)-y_fit(n+1:)))/sum(y(n+1:)), a_lambda, &
              (y(j)-y_fit(j), j = n+1, n+min(n_data, 10)), '|', &
              sum(abs(y(n+1:n+10)-y_fit(n+1:n+10)))/10, sum(abs(y(n+11:n+100)-y_fit(n+11:n+100)))/90, sum(abs(y(n+101:n+1000)-y_fit(n+101:n+1000)))/900

  if (merit_now < 0.99 * merit_min) then
    merit_min = merit_now
    ix_min = i
  endif

  if (i - ix_min > n_cycles .or. a_lambda > 1d10) then
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
  print '(a, i0, a, 4f18.11)', '  mode(', i, ') =', mode(i)%amp, mode(i)%damp, mode(i)%k, mode(i)%phi
enddo

!----------------------------------------------------------------------------------
contains

subroutine mrq_func (a, y_fit, dy_da, status)

real(rp), intent(in) :: a(:)
real(rp), intent(out) :: y_fit(:)
real(rp), intent(out) :: dy_da(:, :)
real(rp) damp, c, s, phi

integer status
integer j, iv, id, iid, iy, iv0

!

status = 0
dy_da = 0
y_fit = 0

do iv = 1, n_mode
  iv0 = n_var_in_mode * (iv-1)
  if (a(iv0+2) < 0) then   ! negative damp factor not allowed.
    status = 1
    return
  endif
enddo

do iv = 1, n_mode
  iv0 = n_var_in_mode * (iv-1)
  if (use_phase) then
    mode(iv) = wake_mode_struct(a(iv0+1), a(iv0+2), a(iv0+3), a(iv0+4))
  else
    mode(iv) = wake_mode_struct(a(iv0+1), a(iv0+2), a(iv0+3), 0.0_rp)
  endif
  y_fit(iv) = a(iv0+1)
  dy_da(iv,iv0+1) = 1
enddo


do id = 1, n_data
  iy = id + n_mode
  iid = id + ix_data_start - 1
  do iv = 1, n_mode
    if (use_phase) then
      phi = a(iv0+4)
    else
      phi = 0
    endif
    iv0 = n_var_in_mode * (iv-1)
    damp = exp(-iid * a(iv0+2)) 
    s = sin(iid * a(iv0+3) + phi)
    c = cos(iid * a(iv0+3) + phi)
    y_fit(iy) = y_fit(iy) + a(iv0+1)   * damp * s
    dy_da(iy, iv0+1) =                   damp * s
    dy_da(iy, iv0+2) = -iid * a(iv0+1) * damp * s
    dy_da(iy, iv0+3) =  iid * a(iv0+1) * damp * c
    if (use_phase) then
      dy_da(iy, iv0+4) =        a(iv0+1) * damp * c
    endif
  enddo
enddo

end subroutine mrq_func

!----------------------------------------------------------------------------------
! contains

function merit() result (this_merit)

real(rp) this_merit
integer iv, id, iid
real(rp) y_fit, damp, s

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
    s = sin(iid * mode(iv)%k + mode(iv)%phi)
    y_fit = y_fit + mode(iv)%amp * damp * s
  enddo
  this_merit = this_merit + (1 / sqrt(real(iid, rp)) - y_fit)**2
enddo

end function merit

end program
