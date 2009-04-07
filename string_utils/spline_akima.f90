!+
! Subroutine spline_akima (spline, stat)
!
! Given a set of (x,y) points we want to interpolate between the points.
! This subroutine computes the semi-hermite cubic spline developed by 
! Hiroshi Akima. The spline goes thorugh all the points (that is, it is 
! not a smoothing spline). For interpolation use:
!           SPLINE_EVALUATE
!
! Reference: 
!   H Akima, "A New Method of Interpolation and Smooth Curve Fitting Based 
!   on Local Procedures", J. Assoc. Comp. Mach., Vol 17(4), 589-602 (1970).
!
! Modules used:
!   use dcslib
!
! Input:
!   spline(:) -- Spline_struct: 
!     %x  -- X-component of a point. Note: points must be in assending order.
!     %y  -- Y-component of a point.
!
! Output:
!   spline(:) -- Spline_struct:
!     %coef(0:3)  -- Spline coefficients at a point.
!   stat      -- Logical: Set .false. if something is wrong (like less than 2
!                    points used).
!
!-

#include "CESR_platform.inc"

subroutine spline_akima (spline, ok)

  use dcslib_struct
  use dcslib_interface, except => spline_akima
  use precision_def

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
    if (spline(i-1)%x .ge. spline(i)%x) then
      print *, 'ERROR IN SPLINE_AKIMA: DATA POINTS NOT IN ASENDING ORDER!'
      print *, i-1, spline(i-1)%x, spline(i)%y
      print *, i, spline(i)%x, spline(i)%y
      return
    endif
  enddo

  spline(:)%coef(0) = spline(:)%y  ! spline passes through all the data points

! special case for 2 two points: use a straight line

  if (nmax .eq. 2) then
    spline(1)%coef(1) = (spline(2)%y - spline(1)%y) / &
                                  (spline(2)%x - spline(1)%x)
    spline(1)%coef(2:3) = 0
    return
  endif

! special case for 3 points: use a quadratic

  if (nmax .eq. 3) then
    y21 = spline(2)%y - spline(1)%y
    y32 = spline(3)%y - spline(2)%y
    x21 = spline(2)%x - spline(1)%x
    x32 = spline(3)%x - spline(2)%x
    x221 = spline(2)%x**2 - spline(1)%x**2
    x232 = spline(3)%x**2 - spline(2)%x**2
    spline(1)%coef(2) = 2 * (y21*x32 - y32*x21) / (x221*x32 - x232*x21)
    spline(2)%coef(1) = (x32*y21/x21 + x21*y32/x32) / (x32 + x21) 
    spline(1)%coef(1) = spline(2)%coef(1) - spline(1)%coef(2) * x21
    spline(2)%coef(2) = spline(1)%coef(2)
    spline(1:2)%coef(3) = 0
    return
  endif

! load coef0 and calc spline at ends

  end(0:3) = spline(4:1:-1)
  call end_spline_calc (end)
  spline(1)%coef(1) = end(3)%coef(1)
  spline(2)%coef(1) = end(2)%coef(1)

  end(0:3) = spline(nmax-3:nmax)
  call end_spline_calc (end)
  spline(nmax)%coef(1) = end(3)%coef(1)
  spline(nmax-1)%coef(1) = end(2)%coef(1)

! calc spline everywhere else
       
  do i = 3, nmax-2
    call slope_calc(spline(i-2:i+2))
  enddo

  do i = 1, nmax-1
    call coef23_calc(spline(i:i+1))
  enddo

  ok = .true.
  return                            

contains

!--------------------------------------------------------------------

  subroutine end_spline_calc (end)
                

  use precision_def

    implicit none

    type (spline_struct), target :: end(0:5)
    real(rp), pointer :: x(:), y(:)
    real(rp) rk

!                                                                 

    x => end(1:5)%x
    y => end(1:5)%y

    x(4) = x(3) - x(1) + x(2)
    x(5) = 2*x(3) - x(1)

    rk = (y(3) - y(2)) / (x(3) - x(2)) - (y(2) - y(1)) / (x(2) - x(1))
    y(4) = y(3) + (x(4) - x(3)) * ((y(3) - y(2)) / (x(3) - x(2)) + rk)
    y(5) = y(4) + (x(5) - x(4)) * ((y(4) - y(3)) / (x(4) - x(3)) + rk)

    call slope_calc(end(0:4))
    call slope_calc(end(1:5))

    return

  end subroutine

!--------------------------------------------------------------------

  subroutine slope_calc (spl)


  use precision_def

    implicit none

    type (spline_struct), target :: spl(1:5)
    real(rp), pointer :: xx(:), yy(:)
    real(rp) m(4), m43, m21

!

    xx => spl(:)%x
    yy => spl(:)%y

    m(:) = (yy(2:5) - yy(1:4)) / (xx(2:5) - xx(1:4))

    if (m(1) == m(2) .and. m(3) == m(4)) then  ! special case
      spl(3)%coef(1) = (m(2) + m(3)) / 2
    else
      m43 = abs(m(4) - m(3))
      m21 = abs(m(2) - m(1)) 
      spl(3)%coef(1) = (m43 * m(2) + m21 * m(3)) / (m43 + m21)
    endif

    return

  end subroutine slope_calc

!-----------------------------------------------------------------------

  subroutine coef23_calc (s2)


  use precision_def

    implicit none

    type (spline_struct) :: s2(2)
    real(rp) x21, y21, t1, t2
                            
!

    x21 = s2(2)%x - s2(1)%x
    y21 = s2(2)%y - s2(1)%y
    t1 = s2(1)%coef(1)
    t2 = s2(2)%coef(1)

    s2(1)%coef(2) = (3*y21 / x21 - 2*t1 - t2) / x21
    s2(1)%coef(3) = (t1 + t2 - 2*y21 / x21) / (x21**2)

    return

  end subroutine coef23_calc

end subroutine spline_akima
