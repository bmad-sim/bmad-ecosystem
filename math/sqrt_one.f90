!+
! function sqrt_one(eps) result (ds1)
!
! Routine to calculate Sqrt[1+eps] - 1 to machine precision.
!
! Input:
!   eps -- real(rp):
!
! Output:
!   ds1 -- real(rp): Result.
!-

elemental function sqrt_one(eps) result (ds1)

use precision_def
implicit none

real(rp), intent(in) :: eps
real(rp) ds1

ds1 = eps / (sqrt(1 + eps) + 1)

end function
