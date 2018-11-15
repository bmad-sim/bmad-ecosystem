!+
! function sinhx_x(x) result (y)
!
! Routine to calculate sinh(x) / x
!
! Input:
!   x -- real(rp): Number
!
! Output:
!   y -- real(rp): sinh(x) / x
!-

elemental function sinhx_x(x) result (y)

use precision_def
implicit none

real(rp), intent(in) :: x
real(rp) y

!

if (abs(x) < 1d-8) then
  y = 1
else
  y = sinh(x) / x
endif

end function
