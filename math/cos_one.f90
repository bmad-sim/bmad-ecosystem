!+
! function cos_one(angle) result (cos1)
!
! Routine to calculate Cos[angle] - 1 to near machine precision.
! Maximum loss of precision is about 1 digit.
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
real(rp) cos1, a2
real(rp), parameter :: c2 = 1.0_rp/2.0_rp, c4 = 1.0_rp/24.0_rp, c6 = 1.0_rp/720.0_rp
real(rp), parameter :: c8 = 1.0_rp/40320.0_rp, c10 = 1.0_rp/3628800.0_rp, c12 = 1.0/479001600.0_rp

!

if (abs(angle) < 0.01) then
  a2 = angle**2
  cos1 = (c2 + (c4 + (c6 * a2)) * a2) * a2

elseif (abs(angle) < 0.3) then
  a2 = angle**2
  cos1 = (c2 + (c4 + (c6 + (c8 + (c10 + (c12 * a2)) * a2) * a2) * a2) * a2) * a2

else
  cos1 = cos(angle) - 1
endif

end function
