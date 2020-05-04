!+
! function dasinc(x) result (y)
!
! Routine to calculate d(asin(x) / x)
!
! Input:
!   x -- real(rp): Number
!
! Output:
!   y -- real(rp): d(asin(x) / x)
!-

elemental function dasinc(x) result (y)

use precision_def
implicit none

real(rp), intent(in) :: x
real(rp) y, x2
real(rp), parameter :: c1 = 1.0_rp / 3.0_rp, c2 = 0.3_rp
real(rp), parameter :: c3 = 15.0_rp / 56.0_rp, c4 = 35.0_rp / 144.0_rp

!

x2 = x*x

if (abs(x) < 0.03_rp) then
  y = x * (c1 + x2 * (c2 + x2 * (c3 + x2 * c4)))
else
  y = 1 / (x * sqrt(1 - x2)) - asin(x) / x2
endif

end function
