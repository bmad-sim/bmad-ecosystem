!+
! Subroutine multipole_ab_to_kt (an, bn, knl, tn)
!
! Subroutine to convert ab type multipoles to kt (MAD standard) multipoles.
! Also see: multipole1_ab_to_kt.
!
! Input:
!   an(0:n_pole_maxx) -- Real(rp): Skew multipole component.
!   bn(0:n_pole_maxx) -- Real(rp): Normal multipole component.
!
! Output:
!   knl(0:n_pole_maxx) -- Real(rp): Multitude magnatude.
!   tn(0:n_pole_maxx)  -- Real(rp): Multipole angle.
!-

subroutine multipole_ab_to_kt (an, bn, knl, tn)

use multipole_mod, dummy => multipole_ab_to_kt

implicit none

real(rp) an(0:), bn(0:)
real(rp) knl(0:), tn(0:)

integer n

!

do n = 0, n_pole_maxx
  call multipole1_ab_to_kt (an(n), bn(n), n, knl(n), tn(n))
enddo

end subroutine multipole_ab_to_kt

