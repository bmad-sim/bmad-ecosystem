!+
! Subroutine multipole_ab_to_kt (an, bn, knl, tn)
!
! Subroutine to convert ab type multipoles to kt (MAD standard) multipoles.
!
! Modules needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!   an(0:n_pole_maxx) -- Real: Skew multipole component.
!   bn(0:n_pole_maxx) -- Real: Normal multipole component.
!
! Output:
!   knl(0:n_pole_maxx) -- Real: Multitude magnatude.
!   tn(0:n_pole_maxx)  -- Real: Multipole angle.
!-

!$Id$
!$Log$
!Revision 1.2  2001/09/27 18:31:54  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine multipole_ab_to_kt (an, bn, knl, tn)

  use bmad_struct
  use bmad_interface

  implicit none

  real an(0:n_pole_maxx), bn(0:n_pole_maxx)
  real knl(0:n_pole_maxx), tn(0:n_pole_maxx)
  real n_fact

  integer n

!

  n_fact = 1

  do n = 0, n_pole_maxx

    if (n /= 0) n_fact = n_fact * n

    if (an(n) == 0 .and. bn(n) == 0) then
      knl(n) = 0
      tn(n) = 0
    else
      knl(n)  = n_fact * sqrt(an(n)**2 + bn(n)**2)
      tn(n) = -atan2(an(n), bn(n)) / (n + 1)
    endif

  enddo


end subroutine


