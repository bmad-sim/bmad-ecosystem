!+                       
! Subroutine track_all (ring, orbit_)
!
! Subroutine to track through the ring.
!
! Modules Needed:
!   use bmad
!
! Input:
!   ring      -- Ring_struct: Ring to track through.
!     %param%aperture_limit_on -- Logical: Sets whether TRACK_ALL looks to
!                                   see whether a particle is lost or not
!   orbit_(0) -- Coord_struct: Coordinates at beginning of ring.
!
! Output:
!   ring
!     %param%lost    -- Logical: Set when a particle is lost with the 
!                         aperture limit on.
!     %param%ix_lost -- Integer: set to index of element where particle is lost
!   orbit_(0:*) -- Coord_struct: Orbit array
!
! Note: If x_limit (or y_limit) for an element is zero then TRACK_ALL will take
!       x_limit (or y_limit) as infinite.
!-

!$Id$
!$Log$
!Revision 1.7  2003/05/02 15:44:03  dcs
!F90 standard conforming changes.
!
!Revision 1.6  2003/01/27 14:40:45  dcs
!bmad_version = 56
!
!Revision 1.5  2002/07/16 20:44:02  dcs
!*** empty log message ***
!
!Revision 1.4  2002/06/13 14:54:29  dcs
!Interfaced with FPP/PTC
!
!Revision 1.3  2002/02/23 20:32:26  dcs
!Double/Single Real toggle added
!
!Revision 1.2  2001/09/27 18:31:58  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"



subroutine track_all (ring, orbit_)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ring_struct)  ring
  type (coord_struct)  orbit_(0:)

  integer n, i

  logical debug / .false. /

! track through elements.

  ring%param%lost = .false.

  do n = 1, ring%n_ele_use

    call track1 (orbit_(n-1), ring%ele_(n), ring%param, orbit_(n))

! check for lost particles

    if (ring%param%lost) then
      ring%param%ix_lost = n
      return
    endif

    if (debug) then
      print *, ring%ele_(n)%name
      print *, (orbit_(n)%vec(i), i = 1, 6)
    endif

  enddo

end subroutine
