!+
! Module spline_mod
!
! Module of cubic spline routines.
!-

module spline_mod

use sim_utils

implicit none

! Given a spline s and some point x_eval, y_spline is:
!   y_spline = Sum: s%coef(i) * dx**i, i = [0:3]
! where dx = x_eval - s%x0 

type spline_struct
  real(rp) :: x0 = 0, y0 = 0     ! Point at start of spline
  real(rp) :: x1 = 0             ! Point at end of spline
  real(rp) :: coef(0:3) = 0      ! coefficients for cubic spline
end type

private akima_spline_coef23_calc, akima_spline_slope_calc, bracket_index_for_spline

contains

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine reallocate_spline (spline, n, n_min, exact)
!
! Subroutine to allocate an allocatable spline_struct array.
! The data of the array is preserved but data at the end of the
! array will be lost if n is less than the original size of the array
! 
!
! Input:
!   spline(:)   -- spline_struct, allocatable: Spline to reallocate.
!   n           -- integer: Upper bound needed for 1-dimensional arrays.
!   n_min       -- integer, optional: Lower bound of spline array. Default is 1.
!   exact       -- logical, optional: If present and False then the size of 
!                    the output array is permitted to be larger than n. Default is True.
!
! Output:
!   spline(:)   -- spline_struct, allocatable: Allocated spline.
!-

subroutine reallocate_spline (spline, n, n_min, exact)

type (spline_struct), allocatable :: spline(:)
type (spline_struct), allocatable :: temp_spline(:)

integer, optional :: n_min
integer :: n, n1, n2
integer i, n1_old, n2_old, n1_save, n2_save

logical, optional :: exact

character(*), parameter :: r_name = 'reallocate_spline'

!

n1 = integer_option(1, n_min); n2 = n1 + n - 1

if (allocated (spline)) then
  n1_old = lbound(spline,1); n2_old = ubound(spline,1)
  if (.not. logic_option(.true., exact) .and. n1_old <= n1 .and. n2 <= n2_old .and. n1 == n1_old) return
  call move_alloc(spline, temp_spline)
  allocate (spline(n1:n2))
  n1_save = max(n1, n1_old); n2_save = min(n2, n2_old)
  spline(n1_save:n2_save) = temp_spline(n1_save:n2_save)
  deallocate (temp_spline)

else
  allocate (spline(n1:n2))
endif

end subroutine reallocate_spline

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!+
! Function create_a_spline (r0, r1, slope0, slope1) result (spline)
!
! Routine to create a single spline given end point positions and slopes.
! The spline will pass through the data points and have the given slopes
! at these points.
!
! Modules used:
!   use spline_mod
!
! Input:
!   r0(2)   -- real(rp): Start (x, y) point.
!   r1(2)   -- real(rp): End (x, y) point.
!   slope0  -- real(rp): Starting slope.
!   slope1  -- real(rp): End slope.
!
! Output:
!   spline  -- spline_struct: Spline.
!-

function create_a_spline (r0, r1, slope0, slope1) result (spline)

type (spline_struct) spline
real(rp) r0(:), r1(:), slope0, slope1, dx, dy

character(*), parameter :: r_name = 'create_a_spline'

!

spline%x0 = r0(1)
spline%y0 = r0(2)
spline%x1 = r1(1)

dx = r1(1) - r0(1)
dy = r1(2) - r0(2)

if (dx == 0) then
  call out_io (s_fatal$, r_name, 'X DISTANCE BETWEEN POINTS IS ZERO.')
  if (global_com%exit_on_error) call err_exit
  return
endif

spline%coef(0) = r0(2)
spline%coef(1) = slope0
spline%coef(2) = (3*dy / dx - 2*slope0 - slope1) / dx
spline%coef(3) = (slope0 + slope1 - 2*dy / dx) / (dx**2)

end function create_a_spline

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!+
! Subroutine spline_akima_interpolate (x_knot, y_knot, x, ok, y, dy)
!
! Routine to interpolate using an akima spline. 
!
! When evaluating at enough points, this routine is slower than calling spline_akima to
! first evaluate the spline coefficients and then repeatedly calling spline_evaluate.
!
! The advantage of this routine is that only the (x, y) knot points need to be stored
! and it will be faster if the number of evaluations is small.
!
! This routine will extrapolate past the range of x_knot(:) up to a distance equal to the 
! length between an end point and the point just inside the end point. 
!
! Input:
!   x_knot(:)   -- real(rp): Array of x values for the knot points. 
!                    Must have more than 2 points and be in asending order.
!   y_knot(:)   -- real(rp): Array of y values for the knot points. Must be same size as x_knot(:).
!   x           -- real(rp): Point to evaluate at.
!
! Output:
!   ok        -- Logical: Set .true. if everything ok, That is, x is within the spline range.
!   y         -- Real(rp), optional: Spline interpolation.
!   dy        -- Real(rp), optional: Spline derivative interpolation.
!-

subroutine spline_akima_interpolate (x_knot, y_knot, x, ok, y, dy)

type (spline_struct) spline(3), spline0

real(rp) x_knot(:), y_knot(:)
real(rp) x
real(rp), optional :: y, dy
real(rp) slope(5), slope0, slope1

integer ix0, ix_min, ix_max, n
logical ok

! Check if x value out of bounds.
          
ok = bracket_index_for_spline(x_knot, x, ix0)
if (.not. ok) return

! And evaluate...
! If up to three points

n = size(x_knot)
if (n < 4) then
  spline(1:n)%x0 = x_knot
  spline(1:n)%y0 = y_knot
  call spline_akima (spline(1:n), ok)
  if (.not. ok) return
  call spline_evaluate (spline(1:n), x, ok, y, dy)
  return
endif

! Regular case

if (ix0 == 1) then
  slope(3:5) = (y_knot(2:4) - y_knot(1:3)) / (x_knot(2:4) - x_knot(1:3))
  slope(2) = 2*slope(3) - slope(4)
  slope(1) = 2*slope(2) - slope(3)

elseif (ix0 == 2) then
  if (n == 4) then
    slope(2:4) = (y_knot(2:4) - y_knot(1:3)) / (x_knot(2:4) - x_knot(1:3))
    slope(5) = 2*slope(4) - slope(3)
  else
    slope(2:5) = (y_knot(2:5) - y_knot(1:4)) / (x_knot(2:5) - x_knot(1:4))
  endif
  slope(1) = 2*slope(2) - slope(3)

elseif (ix0 == n - 2) then
  slope(1:4) = (y_knot(n-3:n) - y_knot(n-4:n-1)) / (x_knot(n-3:n) - x_knot(n-4:n-1))
  slope(5) = 2*slope(4) - slope(3)

elseif (ix0 == n - 1) then
  slope(1:3) = (y_knot(n-2:n) - y_knot(n-3:n-1)) / (x_knot(n-2:n) - x_knot(n-3:n-1))
  slope(4) = 2*slope(3) - slope(2)
  slope(5) = 2*slope(4) - slope(3)

else
  slope(1:5) = (y_knot(ix0-1:ix0+3) - y_knot(ix0-2:ix0+2)) / (x_knot(ix0-1:ix0+3) - x_knot(ix0-2:ix0+2))
endif

!

slope0 = this_slope_calc(slope(1:4))
slope1 = this_slope_calc(slope(2:5))
spline0 = create_a_spline ([x_knot(ix0), y_knot(ix0)], [x_knot(ix0+1), y_knot(ix0+1)], slope0, slope1)

if (present(y))  y  = spline1 (spline0, x)
if (present(dy)) dy = spline1 (spline0, x, 1)

!----------------------------------------------------------
contains

function this_slope_calc(m) result (this_slope)

real(rp) m(4), this_slope
real(rp) m43, m21

!

if (m(1) == m(2) .and. m(3) == m(4)) then  ! special case
  this_slope = (m(2) + m(3)) / 2
else
  m43 = abs(m(4) - m(3))
  m21 = abs(m(2) - m(1)) 
  this_slope = (m43 * m(2) + m21 * m(3)) / (m43 + m21)
endif

end function this_slope_calc

end subroutine spline_akima_interpolate 

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!+
! Subroutine spline_evaluate (spline, x, ok, y, dy)
!
! Subroutine to evalueate a spline at a set of points.
!
! A point outside of the range of knot points is an error.
! Also see: 
!   spline1
!   spline_akima_interpolate
!
! A spline may be generated using, for example, the spline_akima routine.
!
! Modules used:
!   use spline_mod
!
! Input:
!   spline(:) -- Spline_struct: Spline structure.
!   x         -- Real(rp): point for evaluation.
!
! Output:
!   ok        -- Logical: Set .true. if everything ok. That is, x is within the spline range.
!   y         -- Real(rp), optional: Spline interpolation.
!   dy        -- Real(rp), optional: Spline derivative interpolation.
!-

subroutine spline_evaluate (spline, x, ok, y, dy)

type (spline_struct), target :: spline(:)

real(rp) :: x
real(rp), optional :: y, dy

integer ix0, ix_max
logical ok

! Check if x value out of bounds.
          
ok = bracket_index_for_spline(spline%x0, x, ix0)
if (.not. ok) return

! And evaluate

if (present(y))  y  = spline1 (spline(ix0), x)
if (present(dy)) dy = spline1 (spline(ix0), x, 1)

end subroutine spline_evaluate

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!+
! Function bracket_index_for_spline (x_knot, x, ix0) result (ok)
!
! Routine for internal use only.
!-

function bracket_index_for_spline (x_knot, x, ix0) result (ok)

real(rp) x_knot(:), x

integer ix0, ix_max
logical ok

character(*), parameter :: r_name = 'bracket_index_for_spline'

!

ok = .false.
ix_max = size(x_knot)

if (x < x_knot(1) - (x_knot(2) - x_knot(1))) then
  call out_io (s_error$, r_name, 'X EVALUATION POINT (\es12.4\) IS MUCH LESS THAN LOWER BOUND OF SPLINE INTERVAL (\es12.4\)', &
                                 r_array = [x, x_knot(1)])
  return
endif
                              
if (x > x_knot(ix_max) + (x_knot(ix_max) - x_knot(ix_max-1))) then
  call out_io (s_error$, r_name, 'X EVALUATION POINT (\es12.4\) IS MUCH GREATER THAN UPPER BOUND OF SPLINE INTERVAL (\es12.4\) ', &
                                 r_array = [x, x_knot(ix_max)])
  return
endif

ix0 = bracket_index (x, x_knot, 1)
if (ix0 == 0) ix0 = 1
if (ix0 == ix_max) ix0 = ix_max - 1

ok = .true.

end function bracket_index_for_spline

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!+
! Function spline1 (a_spline, x, n) result (y)
!
! Function for spline evaluation using a single spline (instead of a spline array).
! Also see: 
!   spline_evaluate
!   spline_akima_interpolate
!
! Modules used:
!   use spline_mod
!
! Input:
!   a_spline  -- spline_struct: Single spline structure.
!   x         -- real(rp): Point for evaluation.
!   n         -- integer, optional: Output derivative order. May be -1, 0, 1, 2, or 3. Default is 0.
!                   n = -1 => output is integral of y from a_spline%x0 to x.
!                   n = 1 => output is dy/dx, n = 2 => output is d^2y/dx^2, etc.
!
! Output:
!   y         -- real(rp), optional: Interpolated spline value or derivative.
!-

function spline1 (a_spline, x, n) result (y)

type (spline_struct), target :: a_spline

real(rp) :: x, y
real(rp) :: c(0:3)
real(rp) dx       
real(rp), parameter :: a1 = 0.5_rp, a2 = 1.0_rp / 3.0_rp, a3 = 0.25_rp

integer, optional :: n

character(16) :: r_name = 'spline1'

!

dx = x - a_spline%x0
c = a_spline%coef

select case (integer_option(0, n))
case (-1)
  y = ((((a3 * c(3) * dx) + a2 * c(2)) * dx + a1 * c(1)) * dx + c(0)) * dx

case (0)
  y = (((c(3) * dx) + c(2)) * dx + c(1)) * dx + c(0)

case (1)
  y = ((3*c(3) * dx) + 2*c(2)) * dx + c(1)

case (2)
  y = 6*c(3) * dx + 2*c(2)

case (3)
  y = 6*c(3)

end select

end function spline1

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!+
! Subroutine spline_akima (spline, ok)
!
! Given a set of (x,y) points we want to interpolate between the points.
! This subroutine computes the semi-hermite cubic spline developed by 
! Hiroshi Akima. The spline goes thorugh all the points (that is, it is 
! not a smoothing spline). For interpolation use:
!   spline_evaluate
!   spline_akima_interpolate ! You do not need to call spline_akima if you use this routine.
!
! Reference: 
!   H Akima, "A New Method of Interpolation and Smooth Curve Fitting Based 
!   on Local Procedures", J. Assoc. Comp. Mach., Vol 17(4), 589-602 (1970).
!
! Modules used:
!   use spline_mod
!
! Input:
!   spline(:) -- Spline_struct: 
!     %x0  -- X-component of a point. Note: points must be in assending order.
!     %y0  -- Y-component of a point.
!
! Output:
!   spline(:) -- Spline_struct:
!     %coef(0:3)  -- Spline coefficients at a point.
!   ok        -- Logical: Set .false. if something is wrong (like less than 2 points used).
!
!-

subroutine spline_akima (spline, ok)

type (spline_struct) :: spline(:)

real(rp) x21, x31, x32

logical ok

integer i, nmax

character(*), parameter :: r_name = 'spline_akima'

! init
                     
ok = .false.  ! assume the worst
nmax = size(spline)

if (nmax < 2) then
  call out_io (s_error$, r_name, 'LESS THAN 2 DATA POINTS USED!')
  return
endif

do i = 2, nmax
  if (spline(i-1)%x0 .ge. spline(i)%x0) then
    call out_io (s_error$, r_name, 'DATA POINTS NOT IN ASENDING ORDER!', &
           'FOR POINTS: (\i0\, \es10.2\, \es10.2\), and (\i0\, \es10.2\, \es10.2\)', &
           i_array = [i-1, i], r_array = [spline(i-1)%x0, spline(i)%y0, spline(i)%x0, spline(i)%y0])
    return
  endif
enddo

spline(:)%coef(0) = spline(:)%y0  ! Spline must pass through all the data points
spline(1:nmax)%x1 = [spline(2:nmax)%x0, spline(nmax)%x0]

! special case for 2 two points: use a straight line

if (nmax .eq. 2) then
  spline(1)%coef(1) = (spline(2)%y0 - spline(1)%y0) / (spline(2)%x0 - spline(1)%x0)
  spline(1)%coef(2:3) = 0
  ok = .true.
  return
endif

! special case for 3 points: use a quadratic

if (nmax .eq. 3) then
  x21 = spline(2)%x0 - spline(1)%x0
  x31 = spline(3)%x0 - spline(1)%x0
  x32 = spline(3)%x0 - spline(2)%x0

  spline(1:2)%coef(3) = 0

  spline(1)%coef(2) = spline(1)%y0 / (x21 * x31) - spline(2)%y0 / (x21 * x32) + spline(3)%y0 / (x31 * x32)
  spline(2)%coef(2) = spline(1)%coef(2)

  spline(1)%coef(1) = ((spline(3)%y0 - spline(2)%y0) - spline(1)%coef(2) * (x31**2 - x21**2)) / x32
  spline(2)%coef(1) = ((spline(3)%y0 - spline(1)%y0) - spline(2)%coef(2) * (x32**2 - x21**2)) / x31

  ok = .true.
  return
endif

! load coef0 and calc spline at ends

call end_akima_spline_calc (spline, 0)
call end_akima_spline_calc (spline, 1)

! calc spline everywhere else
     
do i = 3, nmax-2
  call akima_spline_slope_calc(spline(i-2:i+2))
enddo

do i = 1, nmax-1
  call akima_spline_coef23_calc(spline(i:i+1))
enddo

ok = .true.

end subroutine spline_akima

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!+
! Subroutine end_akima_spline_calc (spline, which_end)
!
! Routine to calculate the slopes at the ends of a spline array
!
! Input:
!   spline(:)   -- spline_struct: Array of splines.
!   which_end   -- integer: 0 => calculate slopes for the start end of the array.
!                           1 => calculate slopes for the end end of the array.
!
! Output:
!   spline(:)   -- spline_struct: Array with slopes at end calculated.
!-

subroutine end_akima_spline_calc (spline, which_end)

type (spline_struct) :: spline(:)
type (spline_struct), target :: end(0:5)
real(rp), pointer :: x(:), y(:)
real(rp) rk
integer which_end, nmax

!                                                                 

select case (which_end)
case (0)
  end(0:3) = spline(4:1:-1)
case (1)
  nmax = ubound(spline, 1)
  end(0:3) = spline(nmax-3:nmax)
case default
  call err_exit
end select

!

x => end(1:5)%x0
y => end(1:5)%y0

x(4) = x(3) - x(1) + x(2)
x(5) = 2*x(3) - x(1)

rk = (y(3) - y(2)) / (x(3) - x(2)) - (y(2) - y(1)) / (x(2) - x(1))
y(4) = y(3) + (x(4) - x(3)) * ((y(3) - y(2)) / (x(3) - x(2)) + rk)
y(5) = y(4) + (x(5) - x(4)) * ((y(4) - y(3)) / (x(4) - x(3)) + rk)

call akima_spline_slope_calc(end(0:4))
call akima_spline_slope_calc(end(1:5))

!

select case (which_end)
case (0)
  spline(1)%coef(1) = end(3)%coef(1)
  spline(2)%coef(1) = end(2)%coef(1)
case (1)
  spline(nmax)%coef(1) = end(3)%coef(1)
  spline(nmax-1)%coef(1) = end(2)%coef(1)
end select

end subroutine end_akima_spline_calc 

!--------------------------------------------------------------------
! contains

subroutine akima_spline_slope_calc (spl)

type (spline_struct), target :: spl(1:5)
real(rp), pointer :: xx(:), yy(:)
real(rp) m(4), m43, m21

!

xx => spl(:)%x0
yy => spl(:)%y0

m(:) = (yy(2:5) - yy(1:4)) / (xx(2:5) - xx(1:4))

if (m(1) == m(2) .and. m(3) == m(4)) then  ! special case
  spl(3)%coef(1) = (m(2) + m(3)) / 2
else
  m43 = abs(m(4) - m(3))
  m21 = abs(m(2) - m(1)) 
  spl(3)%coef(1) = (m43 * m(2) + m21 * m(3)) / (m43 + m21)
endif

end subroutine akima_spline_slope_calc

!-----------------------------------------------------------------------
! contains

subroutine akima_spline_coef23_calc (s2)

type (spline_struct) :: s2(2)
real(rp) x21, y21, t1, t2
                        
!

x21 = s2(2)%x0 - s2(1)%x0
y21 = s2(2)%y0 - s2(1)%y0
t1 = s2(1)%coef(1)
t2 = s2(2)%coef(1)

s2(1)%coef(2) = (3*y21 / x21 - 2*t1 - t2) / x21
s2(1)%coef(3) = (t1 + t2 - 2*y21 / x21) / (x21**2)

end subroutine akima_spline_coef23_calc

end module
