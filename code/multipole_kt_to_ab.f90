!+
! Subroutine multipole_kt_to_ab (knl, tn, an, bn)
!
! Subroutine to convert kt (MAD standard) multipoles to ab type multipoles.
!
! Modules needed:
!   use bmad
!
! Input:
!   knl(0:) -- Real(rdef): Multitude magnatude.
!   tn(0:)  -- Real(rdef): Multipole angle.
!
! Output:
!   an(0:) -- Real(rdef): Skew multipole component.
!   bn(0:) -- Real(rdef): Normal multipole component.
!-

!$Id$
!$Log$
!Revision 1.7  2002/07/23 17:13:34  dcs
!Changed to (0:) type arrays
!
!Revision 1.6  2002/07/23 17:05:37  palmer
!Recover lost version of multipole_kt_to_ab.f90 from attic.
!
!Revision 1.4  2002/02/23 20:32:21  dcs
!Double/Single Real toggle added
!
!Revision 1.3  2002/01/08 21:44:41  dcs
!Aligned with VMS version  -- DCS
!
!Revision 1.2  2001/09/27 18:31:54  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"

subroutine multipole_kt_to_ab (knl, tn, an, bn)

  use bmad

  implicit none

  real(rdef) an(0:), bn(0:)
  real(rdef) knl(0:), tn(0:)
  real(rdef) n_fact, angle, kl

  integer n

!

  n_fact = 1

  do n = lbound(an, 1), ubound(an, 1)

    if (n /= 0) n_fact = n_fact * n

    kl = knl(n) / n_fact

    if (kl == 0) then
      an(n) = 0
      bn(n) = 0
    else
      angle = -tn(n) * (n + 1)
      an(n) = kl * sin(angle)
      bn(n) = kl * cos(angle)
    endif

  enddo

end subroutine
