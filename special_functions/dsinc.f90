!+
! function dsinc(x) result (y)
!
! Routine to calculate d(sin(x) / x)
!
! Input:
!   x -- real(rp): Number
!
! Output:
!   y -- real(rp): d(sin(x) / x)
!-

elemental function dsinc(x) result (y)

use precision_def
implicit none

real(rp), intent(in) :: x
real(rp) y, x2
real(rp), parameter :: c1 = -1.0_rp / 3.0_rp, c2 = 1.0_rp / 30.0_rp
real(rp), parameter :: c3 = -1.0_rp / 840.0_rp, c4 = 1.0_rp / 45360.0_rp

!

x2 = x*x

if (abs(x) < 0.1_rp) then
  y = x * (c1 + x2 * (c2 + x2 * (c3 + x2 * c4)))
else
  y = (x*cos(x) - sin(x)) / x2
endif

end function
