!+
! Subroutine track_long (ring, orbit_, ix_start, direction, mats627)
!
! Subroutine to track for 1-turn. This subroutine is ment for long term
! tracking and uses 2nd order transport matrices for some of the tracking.
!
! Moudules Needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!   ring             -- Ring_struct: Ring.
!     %param%aperture_limit_on -- Logical: Sets whether TRACK_LONG looks to
!                                   see whether a particle is lost or not
!   orbit_(ix_start) -- Coord_struct: Coordinates at start of tracking.
!   ix_start         -- Integer: Start index (See Note).
!   direction        -- Integer: Direction to track.
!                            = +1  -> Track forward
!                            = -1  -> Track backward
!   mats627(0:n_ele_maxx) -- Mat627_struct: array of 6x27 2nd order matrices.
!                       Use RING_MAKE_MAT627 to calculate these matrices 
!                       before calling TRACK_LONG.
!
! Output:
!   ring
!     %param%lost    -- Logical: Set when a particle is lost with the aperture
!                        limit on.
!     %param%ix_lost -- Integer: set to index of element where particle is lost
!   orbit_(0:n_ele_maxx) -- Coord_struct: Coordinates of tracked particle. 
!   orbit_(ix_start) -- Coord_struct: Coordinates at end of tracking.
!
! Note: Starting and ending points are just after the elements with index
!   IX_START. For example, if DIRECTION = +1 then the first element
!   tracked through is element ix_start+1. If DIRECTION = -1 then the first
!   element tracked through is element ix_start.
!
! Note: Remember that if you are tracking backward and you have calculated
!   the orbit by tracking forward (with, for example, CLOSED_ORBIT_AT_START)
!   then you need to reverse velocity signs:
!           orbit_(ix_start)%x_vel = -orbit_(ix_start)%x_vel 
!           orbit_(ix_start)%y_vel = -orbit_(ix_start)%y_vel 
!
! Note: If x_limit (or y_limit) for an element is zero then TRACK_LONG will
!   take x_limit (or y_limit) as infinite (this is standard BMAD).
!-

!$Id$
!$Log$
!Revision 1.2  2001/09/27 18:31:59  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine track_long (ring, orbit_, ix_start, direction, mats627)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ring_struct) ring
  type (coord_struct) orbit_(0:*)
  type (ele_struct) ele            
  type (mat627_struct) mats627(*)

  integer ix_start, direction
  integer n, i

  real x_lim, y_lim

  logical debug / .false. /
                            
! track through elements.

  ring%param%lost = .false.

  if (direction == +1) then
    call track_fwd (ix_start+1, ring%n_ele_ring)
    if (ring%param%lost) return
    orbit_(0) = orbit_(ring%n_ele_ring) 
    call track_fwd (1, ix_start)
  elseif (direction == -1) then
    call track_back (ix_start, 1)
    if (ring%param%lost) return
    orbit_(ring%n_ele_ring) = orbit_(0)
    call track_back (ring%n_ele_ring, ix_start+1)
  else
    type *, 'ERROR IN TRACK_LONG: BAD DIRECTION:', direction
    call err_exit
  endif

contains

!--------------------------------------------------------------------------

subroutine track_fwd (ix1, ix2)

  integer ix1, ix2

  do n = ix1, ix2

    call track1_627 (orbit_(n-1), ring%ele_(n), ring%param, &
                                                 mats627(n)%m, orbit_(n))

! check for lost particles

    if (ring%param%lost) then
      ring%param%ix_lost = n
      return
    endif

    if (debug) then
      type *, ring%ele_(n)%name
      type *, (orbit_(n)%vec(i), i = 1, 6)
    endif

  enddo

end subroutine

!--------------------------------------------------------------------------

subroutine track_back (ix1, ix2)

  integer ix1, ix2
  logical reverse
  real mat_save(6,6)

!

  do n = ix1, ix2, -1

    call track1_627 (orbit_(n), ring%ele_(n), ring%param, mats627(n)%m, &
                                                                orbit_(n-1))

! check for lost particles

    if (ring%param%lost) then
      ring%param%ix_lost = n - 1
      return
    endif

    if (debug) then
      type *, ring%ele_(n)%name
      type *, (orbit_(n)%vec(i), i = 1, 6)
    endif

  enddo

end subroutine

end subroutine
