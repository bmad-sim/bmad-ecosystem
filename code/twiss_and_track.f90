!+
! Subroutine twiss_and_track (ring, orb)
!
! Subroutine to calculate the twiss parameters and orbit for a lattice.
! Note: This is not necessarily the fastest way to do things.
!
! For a closed ring the computed orbit will be the closed orbit.
! For an open ring the orbit will be computed using orb(0) as starting conditions.
!
! Modules Needed:
!   use bmad
!
! Input:
!   ring              -- Ring_struct: Input ring holding the lattice.
!     %param%lattice_type     -- Used to determine if lattice is open or closed.
!                                = circular_lattice$ implies a closed lattice.
!                                all others imply an open lattice.
!   orb(0:n_ele_maxx)  -- Coord_struct: Orbit to be computed
!     orb(0)             -- Initial conditions to be used for an open ring.
!     orb(0)%z%vel       -- For a closed ring: Energy at which the closed orbit 
!                             is computed.
!
! Output:
!   ring              -- Ring_struct: Ring with computed twiss parameters.
!   orb(0:n_ele_maxx) -- Coord_struct: Computed orbit.
!-

!$Id$
!$Log$
!Revision 1.6  2002/12/13 17:09:40  dcs
!Modified to do non-closed lattices.
!
!Revision 1.5  2002/07/16 20:44:02  dcs
!*** empty log message ***
!
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
  type (coord_struct) orb(0:), orb0(0:n_ele_maxx)

! We need to know the orbit first before we can compute the linear
! transfer matrices for the elements.
! However closed_orbit_at_start needs some crude notion of the 1-turn transfer
! matrix in order for it to do the calculation.

  if (ring%param%lattice_type == circular_lattice$) then
    call ring_make_mat6 (ring, -1, orb0)
    call twiss_at_start (ring)
    call closed_orbit_at_start (ring, orb(0), 4, .true.)
  endif

  call track_all (ring, orb)

! now we can compute the Twiss parameters.

  call ring_make_mat6 (ring, -1, orb)
  if (ring%param%lattice_type == circular_lattice$) call twiss_at_start (ring)
  call twiss_propagate_all (ring)

end subroutine
