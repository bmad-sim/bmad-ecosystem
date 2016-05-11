module spline_mod

use sim_utils

! Given a spline s and some point x_eval, y_spline is:
!   y_spline = Sum: s%coef(i) * dx**i, i = [0:3]
! where dx = x_eval - s%x0 

type spline_struct
  real(rp) x0, y0       ! Point at start of spline
  real(rp) x1           ! Point at end of spline
  real(rp) coef(0:3)    ! coefficients for cubic spline
end type

private akima_spline_coef23_calc, akima_spline_slope_calc, end_akima_spline_calc 

contains

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!+
! Subroutine create_a_spline (spline, r0, r1, slope0, slope1)
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

subroutine create_a_spline (spline, r0, r1, slope0, slope1)

implicit none

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



end subroutine create_a_spline

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!+
! Subroutine spline_evaluate (spline, x, ok, y, dy)
!
! Subroutine to evalueate a spline at a set of points. 
! Also see spline1
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
!   ok        -- Logical: Set .true. if everything ok
!   y         -- Real(rp), optional: Spline interpolation.
!   dy        -- Real(rp), optional: Spline derivative interpolation.
!
! Note:
!   The point x must lie between spline(1)%x0 and spline(max)%x0
!-

subroutine spline_evaluate (spline, x, ok, y, dy)

implicit none

type (spline_struct), target :: spline(:)

real(rp) :: x
real(rp), optional :: y, dy

real(rp) eps

integer ix0, ix_max
                  
logical ok       
character(16) :: r_name = 'spline_evaluate'

! Check if x value out of bounds.
          
ok = .false.

ix_max = ubound(spline(:), 1)
eps = 1d-6 * (spline(ix_max)%x0 - spline(1)%x0)   ! something small

if (x < spline(1)%x0 - eps) then
  call out_io (s_error$, r_name, 'X EVALUATION POINT LESS THAN LOWER BOUND OF SPLINE INTERVAL')
  return
endif
                              
if (x > spline(ix_max)%x0 + eps) then
  call out_io (s_error$, r_name, 'X EVALUATION POINT GREATER THAN UPPER BOUND OF SPLINE INTERVAL')
  return
endif

! Find correct interval and evaluate

call bracket_index (spline%x0, 1, ix_max, x, ix0)
if (present(y))  y  = spline1 (spline(ix0), x)
if (present(dy)) dy = spline1 (spline(ix0), x, 1)

ok = .true.

end subroutine spline_evaluate

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!+
! Function spline1 (a_spline, x, n) result (y)
!
! Function for spline evaluation using a single spline (instead of a spline array).
! Also see: spline_evaluate
!
! Modules used:
!   use spline_mod
!
! Input:
!   a_spline  -- spline_struct: Single spline structure.
!   x         -- real(rp): Point for evaluation.
!   n         -- integer, optional: Output derivative order. May be 0, 1, 2, or 3. Default is 0.
!                   n = 1 => output is dy/dx, n = 2 => output is d^2y/dx^2, etc.
!
! Output:
!   y         -- real(rp), optional: Interpolated spline value or derivative.
!-

function spline1 (a_spline, x, n) result (y)

implicit none

type (spline_struct), target :: a_spline

real(rp) :: x, y
real(rp) :: c(0:3)
real(rp) dx       

integer, optional :: n

character(16) :: r_name = 'spline1'

!

dx = x - a_spline%x0
c = a_spline%coef

select case (integer_option(0, n))
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
!           spline_evaluate
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

implicit none

type (spline_struct) :: spline(:)
type (spline_struct) :: end(0:5)

real(rp) y21, y32, x21, x32, x221, x232

logical ok

integer i, nmax

! init
                     
ok = .false.  ! assume the worst
nmax = ubound(spline, 1)

if (nmax < 2) then
  print *, 'ERROR IN SPLINE_AKIMA: LESS THAN 2 DATA POINTS USED!'
  return
endif

do i = 2, nmax
  if (spline(i-1)%x0 .ge. spline(i)%x0) then
    print *, 'ERROR IN SPLINE_AKIMA: DATA POINTS NOT IN ASENDING ORDER!'
    print *, i-1, spline(i-1)%x0, spline(i)%y0
    print *, i, spline(i)%x0, spline(i)%y0
    return
  endif
enddo

spline(:)%coef(0) = spline(:)%y0  ! spline passes through all the data points

! special case for 2 two points: use a straight line

if (nmax .eq. 2) then
  spline(1)%coef(1) = (spline(2)%y0 - spline(1)%y0) / (spline(2)%x0 - spline(1)%x0)
  spline(1)%coef(2:3) = 0
  return
endif

! special case for 3 points: use a quadratic

if (nmax .eq. 3) then
  y21 = spline(2)%y0 - spline(1)%y0
  y32 = spline(3)%y0 - spline(2)%y0
  x21 = spline(2)%x0 - spline(1)%x0
  x32 = spline(3)%x0 - spline(2)%x0
  x221 = spline(2)%x0**2 - spline(1)%x0**2
  x232 = spline(3)%x0**2 - spline(2)%x0**2
  spline(1)%coef(2) = 2 * (y21*x32 - y32*x21) / (x221*x32 - x232*x21)
  spline(2)%coef(1) = (x32*y21/x21 + x21*y32/x32) / (x32 + x21) 
  spline(1)%coef(1) = spline(2)%coef(1) - spline(1)%coef(2) * x21
  spline(2)%coef(2) = spline(1)%coef(2)
  spline(1:2)%coef(3) = 0
  ok = .true.
  return
endif

! load coef0 and calc spline at ends

end(0:3) = spline(4:1:-1)
call end_akima_spline_calc (end)
spline(1)%coef(1) = end(3)%coef(1)
spline(2)%coef(1) = end(2)%coef(1)

end(0:3) = spline(nmax-3:nmax)
call end_akima_spline_calc (end)
spline(nmax)%coef(1) = end(3)%coef(1)
spline(nmax-1)%coef(1) = end(2)%coef(1)

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
! Subroutine end_akima_spline_calc (end)
!
! Private routine.
!-

subroutine end_akima_spline_calc (end)

implicit none

type (spline_struct), target :: end(0:5)
real(rp), pointer :: x(:), y(:)
real(rp) rk

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

end subroutine end_akima_spline_calc 

!--------------------------------------------------------------------
! contains

subroutine akima_spline_slope_calc (spl)

implicit none

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

implicit none

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
