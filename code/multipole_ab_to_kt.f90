!+
! Subroutine multipole_ab_to_kt (an, bn, knl, tn)
!
! Subroutine to convert ab type multipoles to kt (MAD standard) multipoles.
!
! Modules needed:
!   use bmad
!
! Input:
!   an(0:n_pole_maxx) -- Real(rdef): Skew multipole component.
!   bn(0:n_pole_maxx) -- Real(rdef): Normal multipole component.
!
! Output:
!   knl(0:n_pole_maxx) -- Real(rdef): Multitude magnatude.
!   tn(0:n_pole_maxx)  -- Real(rdef): Multipole angle.
!-

!$Id$
!$Log$
!Revision 1.5  2002/06/13 14:54:26  dcs
!Interfaced with FPP/PTC
!
!Revision 1.4  2002/02/23 20:32:20  dcs
!Double/Single Real toggle added
!
!Revision 1.3  2002/01/08 21:44:41  dcs
!Aligned with VMS version  -- DCS
!
!Revision 1.2  2001/09/27 18:31:54  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine multipole_ab_to_kt (an, bn, knl, tn)

  use bmad

  implicit none

  real(rdef) an(0:), bn(0:)
  real(rdef) knl(0:), tn(0:)
  real(rdef) n_fact, a, b

  integer n

! use a, b as temp values to avoid problems with a call like:
!   call multipole_ab_to_kt (vec1, vec2, vec1, vec2)

  n_fact = 1

  do n = 0, n_pole_maxx

    if (n /= 0) n_fact = n_fact * n

    a = an(n)
    b = bn(n)

    if (a == 0 .and. b == 0) then
      knl(n) = 0
      tn(n) = 0
    else
      knl(n)  = n_fact * sqrt(a**2 + b**2)
      tn(n) = -atan2(a, b) / (n + 1)
    endif

  enddo

end subroutine


