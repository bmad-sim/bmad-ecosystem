!+
! Subroutine multipole_kt_to_ab (knl, tn, an, bn)
!
! Subroutine to convert kt (MAD standard) multipoles to ab type multipoles.
!
! Modules needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!   knl(0:n_pole_maxx) -- Real: Multitude magnatude.
!   tn(0:n_pole_maxx)  -- Real: Multipole angle.
!
! Output:
!   an(0:n_pole_maxx) -- Real: Skew multipole component.
!   bn(0:n_pole_maxx) -- Real: Normal multipole component.
!-

!$Id$
!$Log$
!Revision 1.2  2001/09/27 18:31:54  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine multipole_kt_to_ab (knl, tn, an, bn)

  use bmad_struct
  use bmad_interface

  implicit none

  real an(0:n_pole_maxx), bn(0:n_pole_maxx)
  real knl(0:n_pole_maxx), tn(0:n_pole_maxx)
  real n_fact, angle

  integer n

!

  n_fact = 1

  do n = 0, n_pole_maxx

    if (n /= 0) n_fact = n_fact * n

    if (knl(n) == 0) then
      an(n) = 0
      bn(n) = 0
    else
      angle = -tn(n) * (n + 1)
      an(n) = knl(n) * sin(angle) / n_fact
      bn(n) = knl(n) * cos(angle) / n_fact
    endif

  enddo

end subroutine


