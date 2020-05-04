!+
! function asinc(x) result (y)
!
! Routine to calculate asin(x) / x
!
! Input:
!   x -- real(rp):
!
! Output:
!   y -- real(rp): asin(x) / x
!-

elemental function asinc(x) result (y)

use precision_def
implicit none

real(rp), intent(in) :: x
real(rp) y
real(rp), parameter :: asinc_eps = sqrt(6.0_rp*epsilon(1.0_rp))

!

if (abs(x) < asinc_eps) then
   y = 1.0_rp + x**2/6.0_rp
else
   y = asin(x)/x
end if

end function
