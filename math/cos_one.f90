!+
! function cos_one(angle) result (cos1)
!
! Routine to calculate Cos[angle] - 1 to machine precision.
!
! This is usful if angle can be near zero where the direct evaluation of
! cos[angle] - 1 is inaccurate.
!
! Input:
!   angle -- real(rp): Angle.
!
! Output:
!   cos1  -- real(rp): Result.
!-

elemental function cos_one(angle) result (cos1)

use precision_def
implicit none

real(rp), intent(in) :: angle
real(rp) cos1

!

cos1 = -2 * sin(angle/2)**2

end function
