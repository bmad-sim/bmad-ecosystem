!+
! function sinc(x) result (y)
!
! Routine to calculate sin(x) / x
!
! Input:
!   x -- real(rp):
!
! Output:
!   y -- real(rp): sin(x) / x
!-

function sinc(x) result (y)

use precision_def
implicit none

real(rp) x, y
real(rp), parameter :: sinc_eps = sqrt(6.0_rp*epsilon(1.0_rp))

!

if (abs(x) < sinc_eps) then
   y = 1.0_rp - x**2/6.0_rp
else
   y = sin(x)/x
end if

end function
