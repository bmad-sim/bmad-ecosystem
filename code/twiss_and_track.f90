!+
! SUBROUTINE TWISS_AND_TRACK (ring, orb)
!
! Subroutine to calculate the twiss and orbit parameters in a closed ring.
! This is not necessarily the fastest routine.
!
! Modules Needed:
!   use bmad
!
! Input:
!   RING              -- Ring_struct: Ring 
!   ORB(0:n_ele_maxx) -- Coord_struct: Orbit
!
! Output:
!   RING              -- Ring_struct: Ring with computed twiss parameters
!   ORB(0:n_ele_maxx) -- Coord_struct: Orbit
!-

!$Id$
!$Log$
!Revision 1.4  2002/06/13 14:54:30  dcs
!Interfaced with FPP/PTC
!
!Revision 1.3  2002/02/23 20:32:27  dcs
!Double/Single Real toggle added
!
!Revision 1.2  2001/09/27 18:31:59  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine twiss_and_track (ring, orb)

  use bmad

  implicit none

  type (ring_struct) ring
  type (coord_struct) orb(0:*), orb0(0:n_ele_maxx)

!

  call ring_make_mat6 (ring, -1, orb0)
  call twiss_at_start (ring)
  call closed_orbit_at_start (ring, orb(0), 4, .true.)
  call track_all (ring, orb)

  call ring_make_mat6 (ring, -1, orb)
  call twiss_at_start (ring)

  call twiss_propagate_all (ring)

end subroutine
