!+
! function sqrt1(eps) result (ds1)
!
! Routine to calculate Sqrt[1+eps] - 1 with a precision of about 1e-14.
!
! Input:
!   eps -- real(rp):
!
! Output:
!   ds1 -- real(rp): Result.
!-

function sqrt1(eps) result (ds1)

use precision_def

real(rp) eps, ds1
real(rp), parameter :: c1 = 1.0_rp / 2.0_rp, c2 = -1.0_rp / 8.0_rp 
real(rp), parameter :: c3 = 1.0_rp / 16.0_rp, c4 = -5.0_rp / 128.0_rp
real(rp), parameter :: c5 = 7.0_rp / 256.0_rp

!

if (abs(eps) < 1d-2) then
  ds1 = (c1 + (c2 + (c3 + (c4 + c5 * eps) * eps) * eps) * eps) * eps
else
  ds1 = sqrt(1 + eps) - 1
endif

end function
