!+
! Subroutine track_reverse (lat, orbit, ix_start, ix_end, ix_branch, track_state)
!
! Subroutine to track from one point in the lattice to another.
! The particle is still moving forward in time while moving in the -s direction.
!
! For tracking backwards in time, see the routine track_many.
!
! Note: The coordinates for tracking backward are the same as for tracking 
!   forward: That is, +z always points in the direction of increasing s, not 
!   in the actual direction that the particle is traveling.
!
! Note: Starting and ending points are just after the elements with index
!   IX_START and IX_END. The first element tracked through is element ix_start.
!
! Note: If needed, the subroutine will track through from the end of the lat
!   to the beginning (or vice versa) to get to the end point. 
!   Also: If IX_START = IX_END then the subroutine will track 1 full turn.
!
! Modules Needed:
!   use bmad
!
! Input:
!   lat             -- lat_struct: Lat to track through.
!     %branch(ix_branch)%param%aperture_limit_on 
!                     -- Logical: Sets whether track_reverse looks to
!                                 see whether a particle is lost or not
!   orbit(ix_start)  -- Coord_struct: Coordinates at start of tracking.
!   ix_start         -- Integer: Start index (See Note).
!   ix_end           -- Integer: End index (See Note).
!   ix_branch        -- Integer, optional: Branch to track. Default is 0 (main lattice).
!
! Output:
!   orbit(0:)    -- Coord_struct: Orbit. In particular orbit(ix_end) is
!                       the coordinates at the end of tracking. 
!   track_state  -- Integer, optional: Set to moving_forward$ if everything is OK.
!                     Otherwise: set to index of element where particle was lost.
!-

subroutine track_reverse (lat, orbit, ix_start, ix_end, ix_branch, track_state)

use bmad_interface, except_dummy => track_reverse
use bookkeeper_mod, only: control_bookkeeper

implicit none

type (lat_struct), target :: lat
type (coord_struct) orbit(0:)
type (branch_struct), pointer :: branch

integer ix_start, ix_end, ix_br, n_ele_track, track_end_state
integer, optional :: ix_branch, track_state

logical :: debug = .false.
logical err

character(16) :: r_name = 'track_reverse'

! init

if (bmad_com%auto_bookkeeper) call control_bookkeeper (lat)

ix_br = integer_option (0, ix_branch)
branch => lat%branch(ix_br)
if (present(track_state)) track_state = moving_forward$
track_end_state = moving_forward$

n_ele_track = branch%n_ele_track

if (orbit(ix_start)%state == not_set$) call init_coord(orbit(ix_start), orbit(ix_start)%vec, &
                                                       branch%ele(ix_start), .true., orbit(ix_start)%species)

! Track backwards

call convert_pc_to (branch%ele(ix_start)%value(p0c$) * (1 + orbit(ix_start)%vec(6)), &
                                                branch%param%particle, beta = orbit(ix_start)%beta)
orbit(ix_start)%p0c = branch%ele(ix_start)%value(p0c$)

if (ix_start > ix_end) then
  call track_back (ix_start, ix_end+1, track_end_state)
  return
else
  call track_back (ix_start, 1, track_end_state)
  if (track_end_state /= moving_forward$) then
    call zero_this_track (ix_end, n_ele_track)
    return
  endif
  orbit(n_ele_track) = orbit(0)
  call track_back (n_ele_track, ix_end+1, track_end_state)
endif

!--------------------------------------------------------------------------
contains

subroutine track_back (ix1, ix2, track_end_state)

type (ele_struct), pointer :: ele
integer i, n, ix1, ix2, ix_last, track_end_state

! track

ix_last = ix2-1  ! last index we expect to track.

do n = ix1, ix2, -1

  ele => branch%ele(n)
  orbit(n)%direction = -1 
  ele%orientation = -ele%orientation

  call track1 (orbit(n), ele, branch%param, orbit(n-1), err_flag = err)

  ele%orientation = -ele%orientation

  ! check for lost particles

  if (.not. particle_is_moving_forward(orbit(n-1)) .or. err) then
    track_end_state = n
    if (present(track_state)) track_state = n

    if (orbit(n-1)%location == upstream_end$) then
      call zero_this_track (ix2-1, n-1)
    else
      call zero_this_track (ix2-1, n-2)
    endif
    exit
  endif

  if (debug) then
    print *, branch%ele(n)%name
    print *, (orbit(n)%vec(i), i = 1, 6)
  endif

enddo

end subroutine

!--------------------------------------------------------------------------
! contains

subroutine zero_this_track (n1, n2)

integer n, n1, n2

!

do n = n1, n2
  if (n == ix_start) cycle  ! never zero starting coords.
  orbit(n)%vec = 0
  orbit(n)%state = not_set$
enddo

end subroutine

end subroutine
