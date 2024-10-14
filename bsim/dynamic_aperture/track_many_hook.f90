!+
! Subroutine track_many_hook (finished, lat, orbit, ix_start, ix_end, direction, ix_branch, track_state)
!
! Routine that can be customized to track from one point in the lat to another.
! This routine is always called by track_many. See the track_many documentation for details.
! 
! This particular version is used by the dynamic_aperture program
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

use da_program_mod

implicit none

type (lat_struct), target :: lat
type (coord_struct) orbit(0:), orb
type (ele_struct), pointer :: ele_start, ele0, ele1
type (branch_struct), pointer :: branch
type (ltt_params_struct), pointer :: lttp
type (ltt_com_struct), pointer :: ltt_com
type (probe) prb

integer ix_start, ix_end, direction, ix_br, n_ele_track, track_end_state
integer, optional :: ix_branch, track_state

logical finished, err

character(*), parameter :: r_name = 'track_many_hook'

!

branch => lat%branch(ix_branch)
lttp => ltt_params_global
ltt_com => ltt_com_global

call ltt_pointer_to_map_ends(lttp, lat, ele_start)

call ltt_setup_high_energy_space_charge(lttp, ltt_com, branch)

! Normal Bmad tracking?

finished = .false.
if (lttp%tracking_method == 'OLD' .or. ltt_com%track_bypass) return
finished = .true.

!

branch => lat%branch(integer_option(0, ix_branch))
orb = orbit(ele_start%ix_ele)

select case (lttp%tracking_method)
case ('BMAD')
  call ltt_track_bmad_single (lttp, ltt_com, branch%ele(ix_start), branch%ele(ix_end), orb)

case ('PTC')
  ele0 => ele_start
  do
    prb = orb%vec
    prb%q%x = [1, 0, 0, 0]  ! Unit quaternion
    ele1 => pointer_to_next_ele(ele0)
    call track_probe (prb, ltt_com%ptc_state, fibre1 = pointer_to_fibre(ele0), fibre2 = pointer_to_fibre(ele1))
    orb%vec = prb%x
    orb%spin = quat_rotate(prb%q%x, orb%spin)
    if (abs(orb%vec(1)) > lttp%ptc_aperture(1) .or. abs(orb%vec(3)) > lttp%ptc_aperture(2) .or. &
                                                   orbit_too_large(orb) .or. prb%u) orb%state = lost$
    if (orb%state /= alive$) exit
    if (ele1%ix_ele == ele_start%ix_ele) exit
    ele0 => ele1
  enddo

case ('MAP')
  call ltt_track_map(lttp, ltt_com, orb)

case default
  print '(a)', 'Unknown tracking_method: ' // lttp%tracking_method
  stop
end select

orbit(ix_end) = orb
if (present(track_state)) then
  if (orb%state == alive$) then
    track_state = moving_forward$
  else
    track_state = ix_end
  endif
endif

end subroutine track_many_hook
