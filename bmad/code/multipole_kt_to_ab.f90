!+
! Subroutine multipole_kt_to_ab (knl, knsl, tn, an, bn)
!
! Subroutine to convert kt (MAD standard) multipoles to ab type multipoles.
! Also see: multipole1_kt_to_ab.
!
! Input:
!   knl(0:)  -- Real(rp): Normal multitude component.
!   knsl(0:) -- Real(rp): Skew multitude component.
!   tn(0:)   -- Real(rp): Multipole angle.
!
! Output:
!   an(0:) -- Real(rp): Skew multipole component.
!   bn(0:) -- Real(rp): Normal multipole component.
!-

subroutine multipole_kt_to_ab (knl, knsl, tn, an, bn)

use multipole_mod, dummy => multipole_kt_to_ab

implicit none

real(rp) an(0:), bn(0:)
real(rp) knl(0:), knsl(0:), tn(0:)

integer n

!

do n = lbound(an, 1), ubound(an, 1)
  call multipole1_kt_to_ab (knl(n), knsl(n), tn(n), n, an(n), bn(n))
enddo

end subroutine multipole_kt_to_ab

