!+
! Subroutine track_many (ring, orbit, ix_start, ix_end, direction)
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
!   orbit(ix_start)  -- Coord_struct: Coordinates at start of tracking.
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
!     %param%ix_lost -- Integer: Index of element where particle is lost.
!     %param%end_lost_at -- Either entrance_end$ or exit_end$.
!   orbit(:)      -- Coord_struct: Orbit. In particular orbit(ix_end) is
!                       the coordinates at the end of tracking. 
!-

#include "CESR_platform.inc"

subroutine track_many (ring, orbit, ix_start, ix_end, direction)

  use bmad_struct
  use bmad_interface, except => track_many
  use bookkeeper_mod, only: control_lord_bookkeeper
  use reverse_mod, only: reverse_ele

  implicit none

  type (ring_struct) ring
  type (coord_struct) orbit(0:)

  integer ix_start, ix_end, direction
  integer n, i

  logical :: debug = .false.

  character(16) :: r_name = 'track_many'

! init

  call control_lord_bookkeeper (ring)
                            
  ring%param%lost = .false.
  ring%param%ix_lost = -1

! track through elements.

  if (direction == +1) then
    if (ix_start < ix_end) then
      call track_fwd (ix_start+1, ix_end)
      return
    else
      call track_fwd (ix_start+1, ring%n_ele_use)
      if (ring%param%lost) then
        call zero_this_track (0, ix_end)
        return
      endif
      orbit(0) = orbit(ring%n_ele_use) 
      call track_fwd (1, ix_end)
    endif

  elseif (direction == -1) then
    if (ix_start > ix_end) then
      call track_back (ix_start, ix_end+1)
      return
    else
      call track_back (ix_start, 1)
      if (ring%param%lost) then
        call zero_this_track (ix_end, ring%n_ele_use)
        return
      endif
      orbit(ring%n_ele_use) = orbit(0)
      call track_back (ring%n_ele_use, ix_end+1)
    endif

  else
    print *, 'ERROR IN TRACK_MANY: BAD DIRECTION:', direction
    call err_exit
  endif

!--------------------------------------------------------------------------
! tracking forward is simple

contains

subroutine track_fwd (ix1, ix2)

  integer ix1, ix2

  do n = ix1, ix2

    call track1 (orbit(n-1), ring%ele_(n), ring%param, orbit(n))

! check for lost particles

    if (ring%param%lost) then
      ring%param%ix_lost = n
      if (ring%param%end_lost_at == exit_end$) then
        call zero_this_track (n+1, ix2)
      elseif (ring%param%end_lost_at == entrance_end$) then
        call zero_this_track (n, ix2)
      else
        call out_io (s_abort$, r_name, 'INTERNAL ERROR')
        call err_exit
      endif
      return
    endif

    if (debug) then
      print *, ring%ele_(n)%name
      print *, (orbit(n)%vec(i), i = 1, 6)
    endif

  enddo

end subroutine

!--------------------------------------------------------------------------
! contains

! reverse_ele is used to reverse an element for tracking backwards.
! However, a reversed element has a different coordinate system so
! we need to transform to the flipped coordinate system, then track, then
! flip back to the standard coord system.

subroutine track_back (ix1, ix2)

  type (ele_struct), save :: ele

  integer ix1, ix2, ix_last 

! flip to reversed coords

  orbit(ix1)%vec(2) = -orbit(ix1)%vec(2)
  orbit(ix1)%vec(4) = -orbit(ix1)%vec(4)
  orbit(ix1)%vec(5) = -orbit(ix1)%vec(5)

! track

  ix_last = ix2-1  ! last index we expect to track.

  do n = ix1, ix2, -1

    ele = ring%ele_(n)
    call reverse_ele (ele)
    call track1 (orbit(n), ele, ring%param, orbit(n-1))

! check for lost particles

    if (ring%param%lost) then
      if (ring%param%end_lost_at == exit_end$) then
        ring%param%end_lost_at = entrance_end$
        call zero_this_track (ix2-1, n-2)
      elseif (ring%param%end_lost_at == entrance_end$) then
        ring%param%end_lost_at = exit_end$
        call zero_this_track (ix2-1, n-1)
      else
        call out_io (s_abort$, r_name, 'INTERNAL ERROR')
        call err_exit
      endif
      ring%param%ix_lost = n 
      ix_last = n-1
      exit
    endif

    if (debug) then
      print *, ring%ele_(n)%name
      print *, (orbit(n)%vec(i), i = 1, 6)
    endif

  enddo

! flip back to normal coords

  orbit(ix_last:ix1)%vec(2) = -orbit(ix_last:ix1)%vec(2)
  orbit(ix_last:ix1)%vec(4) = -orbit(ix_last:ix1)%vec(4)
  orbit(ix_last:ix1)%vec(5) = -orbit(ix_last:ix1)%vec(5)

end subroutine

!--------------------------------------------------------------------------
! contains

subroutine zero_this_track (n1, n2)

  integer n, n1, n2

!

  do n = n1, n2
    if (n == ix_start) cycle  ! never zero starting coords.
    orbit(n)%vec = 0
  enddo

end subroutine
end subroutine
