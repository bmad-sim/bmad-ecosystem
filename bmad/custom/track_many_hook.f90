!+
! Subroutine track_many_hook (finished, lat, orbit, ix_start, ix_end, direction, ix_branch, track_state)
!
! Routine that can be customized to track from one point in the lat to another.
! To use, see the Bmad manual.
!
! Input:
!   lat              -- lat_struct: Lat to track through.
!   orbit(ix_start)  -- Coord_struct: Coordinates at start of tracking.
!   ix_start         -- Integer: Start index (See Note).
!   ix_end           -- Integer: End index (See Note).
!   direction        -- Integer: Direction to track.
!                            = +1 -> Track forward (+s)
!                            = -1 -> Track backward (-s)
!   ix_branch        -- Integer, optional: Branch to track. Default is 0 (main lattice).
!
! Output:
!   finished     -- logical: If set True then track_many will do nothing on the assumumption that this 
!                     routine has done the tracking.
!   orbit(0:)    -- Coord_struct: Orbit.
!   track_state  -- Integer, optional: Set to moving_forward$ if everything is OK.
!                     Otherwise: set to index of element where particle was lost.
!-

subroutine track_many_hook (finished, lat, orbit, ix_start, ix_end, direction, ix_branch, track_state)

use bmad_interface

implicit none

type (lat_struct), target :: lat
type (coord_struct) orbit(0:)
type (branch_struct), pointer :: branch

integer ix_start, ix_end, direction, ix_br, n_ele_track, track_end_state
integer, optional :: ix_branch, track_state

logical :: debug = .false.
logical finished, err

character(*), parameter :: r_name = 'track_many_hook'

!

finished = .false.

end subroutine track_many_hook
