!+
! Subroutine track_many (ring, orbit_, ix_start, ix_end, direction)
!
! Subroutine to track from one point in the ring to another.
!
! Note: A faster way for tracking backward is to use a reversed ring using
!   ring_reverse. See ring_reverse for more details.
!
! Note: Starting and ending points are just after the elements with index
!   IX_START and IX_END. For example, if DIRECTION = +1 then the first element
!   tracked through is element ix_start+1. If DIRECTION = -1 then the first
!   element tracked through is element ix_start.
!
! Note: The coordinates for tracking backward are the same as for tracking 
!   forward: That is, +z always points in the direction of increasing s, not 
!   in the actual direction that the particle is traveling.
!
! Note: If needed the subroutine will track through from the end of the ring
!   to the beginning (or vice versa) to get to the end point. 
!   Also: If IX_START = IX_END then the subroutine will track 1 full turn.
!
! Note: If x_limit (or y_limit) for an element is zero then TRACK_MANY will
!   take x_limit (or y_limit) as infinite (this is standard BMAD).
!
! Modules Needed:
!   use bmad
!
! Input:
!   ring             -- Ring_struct: Ring to track through.
!     %param%aperture_limit_on -- Logical: Sets whether TRACK_MANY looks to
!                                 see whether a particle is lost or not
!   orbit_(ix_start) -- Coord_struct: Coordinates at start of tracking.
!   ix_start         -- Integer: Start index (See Note).
!   ix_end           -- Integer: End index (See Note).
!   direction        -- Integer: Direction to track.
!                            = +1  -> Track forward
!                            = -1  -> Track backward
!
! Output:
!   ring          -- Ring_struct:
!     %param%lost    -- Logical: Set when a particle is lost with the 
!                         aperture limit on.
!     %param%ix_lost -- Integer: set to index of element where particle is lost
!   orbit_(ix_end) -- Coord_struct: Coordinates at end of tracking 
!                       Also: the track between IX_START and IX_END
!                       are filled in.
!-

!$Id$
!$Log$
!Revision 1.12  2003/07/09 01:38:22  dcs
!new bmad with allocatable ring%ele_(:)
!
!Revision 1.11  2003/06/04 17:55:55  dcs
!Eliminated x%pos, x%vel, etc. from coord_struct.
!
!Revision 1.10  2003/05/02 15:44:03  dcs
!F90 standard conforming changes.
!
!Revision 1.9  2003/03/06 17:47:40  dcs
!Hiding PTC stuff.
!
!Revision 1.8  2003/01/27 14:40:45  dcs
!bmad_version = 56
!
!Revision 1.7  2003/01/02 16:19:40  dcs
!*** empty log message ***
!
!Revision 1.4  2002/06/13 14:54:29  dcs
!Interfaced with FPP/PTC
!
!Revision 1.3  2002/02/23 20:32:26  dcs
!Double/Single Real toggle added
!
!Revision 1.2  2001/09/27 18:31:59  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"

subroutine track_many (ring, orbit_, ix_start, ix_end, direction)

  use bmad_struct
  use bmad_interface
  use reverse_mod

  implicit none

  record / ring_struct / ring
  record / coord_struct / orbit_(0:)
  record / ele_struct / ele            

  integer ix_start, ix_end, direction
  integer n, i

  real(rp) x_lim, y_lim

  logical debug / .false. /
                            
! track through elements.

  ring%param%lost = .false.

  if (direction == +1) then
    if (ix_start < ix_end) then
      call track_fwd (ix_start+1, ix_end)
      return
    else
      call track_fwd (ix_start+1, ring%n_ele_ring)
      if (ring%param%lost) return
      orbit_(0) = orbit_(ring%n_ele_ring) 
      call track_fwd (1, ix_end)
    endif

  elseif (direction == -1) then
    if (ix_start > ix_end) then
      call track_back (ix_start, ix_end+1)
      return
    else
      call track_back (ix_start, 1)
      if (ring%param%lost) return
      orbit_(ring%n_ele_ring) = orbit_(0)
      call track_back (ring%n_ele_ring, ix_end+1)
    endif

  else
    print *, 'ERROR IN TRACK_MANY: BAD DIRECTION:', direction
    call err_exit
  endif

contains

!--------------------------------------------------------------------------
! tracking forward is simple

subroutine track_fwd (ix1, ix2)

  integer ix1, ix2

  do n = ix1, ix2

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

!--------------------------------------------------------------------------
! reverse_ele is used to reverse an element for tracking backwards.
! However, a reversed element has a different coordinate system so
! we need to transform to the flipped coordinate system, then track, then
! flip back to the standard coord system.

subroutine track_back (ix1, ix2)

  type (ele_struct) ele

  integer ix1, ix2, ix_last 

! flip to reversed coords

  orbit_(ix1)%vec(2) = -orbit_(ix1)%vec(2)
  orbit_(ix1)%vec(4) = -orbit_(ix1)%vec(4)
  orbit_(ix1)%vec(5) = -orbit_(ix1)%vec(5)

! track

  ix_last = ix2-1  ! last index we expect to track.

  do n = ix1, ix2, -1

    ele = ring%ele_(n)
    call reverse_ele (ele)
    call track1 (orbit_(n), ele, ring%param, orbit_(n-1))

! check for lost particles

    if (ring%param%lost) then
      ring%param%ix_lost = n - 1
      ix_last = n-1
      exit
    endif

    if (debug) then
      print *, ring%ele_(n)%name
      print *, (orbit_(n)%vec(i), i = 1, 6)
    endif

  enddo

! flip back to normal coords

  orbit_(ix_last:ix1)%vec(2) = -orbit_(ix_last:ix1)%vec(2)
  orbit_(ix_last:ix1)%vec(4) = -orbit_(ix_last:ix1)%vec(4)
  orbit_(ix_last:ix1)%vec(5) = -orbit_(ix_last:ix1)%vec(5)

end subroutine

end subroutine
