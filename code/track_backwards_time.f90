! NOTE: THIS ROUTINE IS NOT CURRENTLY OPERATIONAL!
!       DO NOT USE!
!       IF YOU ARE INTERESTED IN USING THIS PLEASE CONTACT DAVID SAGAN.
!+
! Subroutine track_backwards_time (lat, orbit, ix_start, ix_end, direction, ix_branch, track_state)
!
! Subroutine to track from one point in the lattice to another going backwards in time.
! For forward time tracking use:
!   track_many
!   track_all
!
! Tracking with direction = -1 means the particle is traveling in the negative s-direction
! as well as traveling backwards in time.
!
! Note: The particle coordinates are independent of the direction of travel. 
! That is, it is always true that:
!   vec(2) = dx/ds    ! positive dx/ds means negative dx/dt for direction = -1
!   vec(4) = dy/ds
!   vec(5) = z        ! +z points in the +s direction.
!
! Note: Starting and ending points are at the downstream (+s) end of elements with index
!   ix_start and ix_end. Thus:
!     Direction  Starting-Orbit    First-Element-Tracked    Ending-Orbit    Last-Element-Tracked
!       +1       orbit(ix_start)   ix_start+1               orbit(ix_end)   ix_end
!       -1       orbit(ix_start)   ix_start                 orbit(ix_end)   ix_end+1
! Exception: 
!     If direction =  1 and ix_start = last-lattice-element -> First tracked = Element #1.
!     If direction = -1 and ix_end   = last-lattice-element -> Last tracked = Element #1.
!
! Note: If needed the subroutine will track through from the end of the lat
!   to the beginning (or vice versa) to get to the end point. 
!   Also: if ix_start = ix_end then the subroutine will track 1 full turn.
!
! Modules Needed:
!   use bmad
!
! Input:
!   lat             -- lat_struct: Lat to track through.
!   orbit(ix_start)  -- Coord_struct: Coordinates at start of tracking.
!   ix_start         -- Integer: Start index (See Note).
!   ix_end           -- Integer: End index (See Note).
!   direction        -- Integer: Direction to track.
!                            = +1  -> Track forward
!                            = -1  -> Track backward
!   ix_branch        -- Integer, optional: Branch to track. Default is 0 (main lattice).
!
! Output:
!   orbit(0:)    -- Coord_struct: Orbit.
!   track_state  -- Integer, optional: Set to moving_forward$ if everything is OK.
!                     Otherwise: set to index of element where particle was lost.
!-

subroutine track_backwards_time (lat, orbit, ix_start, ix_end, direction, ix_branch, track_state)

use bmad_struct
use bmad_interface, except_dummy => track_backwards_time
use bookkeeper_mod, only: control_bookkeeper

implicit none

type (lat_struct), target :: lat
type (coord_struct) orbit(0:)
type (branch_struct), pointer :: branch

integer ix_start, ix_end, direction, ix_br, n_ele_track, track_end_state
integer, optional :: ix_branch, track_state

logical :: debug = .false.
logical err

character(*), parameter :: r_name = 'track_backwards_time'

! init

if (bmad_com%auto_bookkeeper) call control_bookkeeper (lat)

ix_br = integer_option (0, ix_branch)
branch => lat%branch(ix_br)
if (present(track_state)) track_state = moving_forward$
track_end_state = moving_forward$

n_ele_track = branch%n_ele_track

! Track forward through the elements.

if (direction == +1) then

  if (ix_start < ix_end) then
    call track_fwd (ix_start+1, ix_end, track_end_state)
    return
  else
    call track_fwd (ix_start+1, n_ele_track, track_end_state)
    if (track_end_state /= moving_forward$) then
      call set_orbit_to_zero (orbit, 0, ix_end)
      return
    endif
    orbit(0) = orbit(n_ele_track) 
    call track_fwd (1, ix_end, track_end_state)
  endif

! Track backwards

elseif (direction == -1) then

  if (ix_start > ix_end) then
    call track_reverse_s (ix_start, ix_end+1, track_end_state)
    return
  else
    call track_reverse_s (ix_start, 1, track_end_state)
    if (track_end_state /= moving_forward$) then
      call set_orbit_to_zero (orbit, ix_end, n_ele_track)
      return
    endif
    orbit(n_ele_track) = orbit(0)
    call track_reverse_s (n_ele_track, ix_end+1, track_end_state)
  endif

else
  call out_io (s_fatal$, r_name, 'BAD DIRECTION: \i0\ ', i_array = [direction])
  if (global_com%exit_on_error) call err_exit
endif

!--------------------------------------------------------------------------
! tracking forward 

contains

subroutine track_fwd (ix1, ix2, track_end_state)

type (ele_struct), pointer :: ele
integer i, n, ix1, ix2, track_end_state

do n = ix1, ix2

  ele => branch%ele(n)
  call track1_backwards_time (orbit(n-1), ele, branch%param, orbit(n), err_flag = err)

  ! check for lost particles

  if (.not. particle_is_moving_forward(orbit(n)) .or. err) then
    track_end_state = n
    if (present(track_state)) track_state = n

    if (orbit(n)%location == upstream_end$) then
      call set_orbit_to_zero (orbit, n, ix2, ix_start)
    else
      call set_orbit_to_zero (orbit, n+1, ix2, ix_start)
    endif
    return
  endif

  if (debug) then
    print *, ele%name
    print *, (orbit(n)%vec(i), i = 1, 6)
  endif

enddo

end subroutine track_fwd

!--------------------------------------------------------------------------
! contains

subroutine track_reverse_s (ix1, ix2, track_end_state)

type (ele_struct), pointer :: ele

integer i, n, ix1, ix2, ix_last, track_end_state

! track

ix_last = ix2-1  ! last index we expect to track.

do n = ix1, ix2, -1

  ele => branch%ele(n)
  orbit(n)%direction = -1
  ele%orientation = -ele%orientation
  call track1_backwards_time (orbit(n), ele, branch%param, orbit(n-1), err_flag = err)
  ele%orientation = -ele%orientation

  ! check for lost particles

  if (.not. particle_is_moving_forward(orbit(n-1)) .or. err) then
    track_end_state = n
    if (present(track_state)) track_state = n

    if (orbit(n-1)%location == upstream_end$) then
      call set_orbit_to_zero (orbit, ix2-1, n-1, ix_start)
    else
      call set_orbit_to_zero (orbit, ix2-1, n-2, ix_start)
    endif
    exit
  endif

  if (debug) then
    print *, branch%ele(n)%name
    print *, (orbit(n)%vec(i), i = 1, 6)
  endif

enddo

end subroutine track_reverse_s

end subroutine
