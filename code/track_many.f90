!+
! Subroutine track_many (ring, orbit_, ix_start, ix_end, direction)
!
! Subroutine to track from one point in the ring to another.
!
! Modules Needed:
!   use bmad_struct
!   use bmad_interface
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
!
! Note: Starting and ending points are just after the elements with index
!   IX_START and IX_END. For example, if DIRECTION = +1 then the first element
!   tracked through is element ix_start+1. If DIRECTION = -1 then the first
!   element tracked through is element ix_start.
!
! Note: If needed the subroutine will track through from the end of the ring
!   to the beginning (or vice versa) to get to the end point. 
!   Also: If IX_START = IX_END then the subroutine will track 1 full turn.
!
! Note: Remember that if you are tracking backward and you have calculated
!   the orbit by tracking forward (with, for example, CLOSED_ORBIT_AT_START)
!   then you need to reverse velocity signs:
!           orbit_(ix_start)%x_vel = -orbit_(ix_start)%x_vel 
!           orbit_(ix_start)%y_vel = -orbit_(ix_start)%y_vel 
!
! Note: If x_limit (or y_limit) for an element is zero then TRACK_MANY will
!   take x_limit (or y_limit) as infinite (this is standard BMAD).
!-


subroutine track_many (ring, orbit_, ix_start, ix_end, direction)

  use bmad_struct
  use bmad_interface

  implicit none

  record / ring_struct / ring
  record / coord_struct / orbit_(0:*)
  record / ele_struct / ele            

  integer ix_start, ix_end, direction
  integer n, i

  real x_lim, y_lim

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
    type *, 'ERROR IN TRACK_MANY: BAD DIRECTION:', direction
    call err_exit
  endif

contains

!--------------------------------------------------------------------------

subroutine track_fwd (ix1, ix2)

  integer ix1, ix2

  do n = ix1, ix2

    if (ring%ele_(n)%key == coil$) then
      call coil_track (orbit_(n-1), n, ring, orbit_(n))
    else
      call track1 (orbit_(n-1), ring%ele_(n), ring%param, orbit_(n))
    endif

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
! TRACK1 thinks it is tracking "forward". However,
! tracking backward is equivalent to tracking forwards with the coordinate
! transformation s -> -s. Thus any longitudinal component of the
! magnetic field in an element must flipped in sign. 

subroutine track_back (ix1, ix2)

  integer ix1, ix2
  logical reverse
  real mat_save(6,6)

!

  do n = ix1, ix2, -1

    reverse = .false.
    if (ring%ele_(n)%key == solenoid$ .or. ring%ele_(n)%key == sol_quad$) &
                                                            reverse = .true.

    if (ring%ele_(n)%key == coil$) then
      call coil_track (orbit_(n), n, ring, orbit_(n-1))
    else
      if (reverse) then
        ring%ele_(n)%value(ks$) = -ring%ele_(n)%value(ks$)
        mat_save = ring%ele_(n)%mat6
        call make_mat6 (ring%ele_(n), ring%param, orbit_(n), orbit_(n))
      endif
      call track1 (orbit_(n), ring%ele_(n), ring%param, orbit_(n-1))
      if (reverse) then
        ring%ele_(n)%value(ks$) = -ring%ele_(n)%value(ks$)
        ring%ele_(n)%mat6 = mat_save
      endif
    endif

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
