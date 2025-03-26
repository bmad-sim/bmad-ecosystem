!+
! Subroutine multipole1_kt_to_ab (knl, knsl, tn, n, an, bn)
!
! Subroutine to convert kt (MAD standard) multipoles to ab type multipoles.
! Also see: multipole_kt_to_ab.
!
! Input:
!   knl  -- Real(rp): Normal multitude component.
!   knsl -- Real(rp): Skew multitude component.
!   tn   -- Real(rp): Multipole angle.
!   n    -- Integer: Multipole order.
!
! Output:
!   an -- Real(rp): Skew multipole component.
!   bn -- Real(rp): Normal multipole component.
!-

subroutine multipole1_kt_to_ab (knl, knsl, tn, n, an, bn)

use sim_utils

implicit none

real(rp) an, bn
real(rp) knl, knsl, tn
real(rp) angle

integer n

!

if (knl == 0 .and. knsl == 0) then
  bn = 0
  an = 0
else
  angle = -tn * (n + 1)
  bn = (knl * cos(angle) - knsl * sin(angle)) / factorial(n)
  an = (knl * sin(angle) + knsl * cos(angle)) / factorial(n)
endif

end subroutine multipole1_kt_to_ab

