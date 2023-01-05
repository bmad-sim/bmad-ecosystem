!+
! Subroutine track_from_s_to_s (lat, s_start, s_end, orbit_start, orbit_end, all_orb, ix_branch, track_state)
!
! Routine to track a particle between two s-positions.
! If the particle is lost in tracking, end_orb will hold the coordinates at the point of loss.
! To get the twiss parameters at the ending s-position use twiss_and_track_at_s.
!
! Input:
!   lat         -- lat_struct: Lattice to track through
!   s_start     -- Real(rp): Starting s-position.
!   s_end       -- Real(rp): Ending s-position. If <= s_start then will wrap
!   orbit_start -- coord_struct: Starting coordinates.
!   ix_branch   -- Integer, optional: Lattice branch index. Default is 0 (main branch).
!
! Output:
!   orbit_end   -- coord_struct: Ending coordinates.
!   all_orb(0:) -- coord_struct, allocatable, optional: If present then the orbit at the exit ends
!                   of the elements tracked through will be recorded in this structure. 
!   track_state -- Integer, optional: Set to moving_forward$ if everything is OK.
!                     Otherwise: set to index of element where particle was lost.
!-   

subroutine track_from_s_to_s (lat, s_start, s_end, orbit_start, orbit_end, all_orb, ix_branch, track_state)

use bmad_interface, dummy => track_from_s_to_s

implicit none

type (lat_struct), target :: lat
type (coord_struct) orbit_start, orbit_end
type (branch_struct), pointer :: branch
type (coord_struct), optional, allocatable :: all_orb(:)
type (ele_struct), pointer :: ele

real(rp) s_start, s_end
real(rp) s0

integer, optional :: ix_branch, track_state
integer ix_start, ix_end
integer ix_ele

logical err

character(40), parameter :: r_name = 'track_from_s_to_s'

! Easy case & error check

branch => lat%branch(integer_option(0, ix_branch))
if (present(track_state)) track_state = moving_forward$

if (s_start == s_end .and. branch%param%geometry == open$) then
  orbit_end = orbit_start
  return
endif

if (s_end < s_start .and. branch%param%geometry == open$) then
  call out_io (s_abort$, r_name, 'S_END < S_START WITH AN OPEN LATTICE.')
  if (global_com%exit_on_error) call err_exit
endif

! Find elements corresponding to s_start and s_stop

if (s_start == branch%ele(0)%s) then
  ix_start = 1
else
  ix_start = element_at_s (branch, s_start, .true.)
endif

s0 = branch%ele(ix_start)%s_start
ix_end = element_at_s (lat, s_end, .true., ix_branch)

! Track within a single element case

ele => branch%ele(ix_start)

if (s_end > s_start .and. ix_start == ix_end) then
  call twiss_and_track_intra_ele (ele, branch%param, s_start-s0, s_end-s0, &
                                               .true., .true., orbit_start, orbit_end)
  if (.not. particle_is_moving_forward(orbit_end) .and. present(track_state)) track_state = ix_start
  return
endif

! Track to end of current element

call twiss_and_track_intra_ele (ele, branch%param,  &
            s_start-s0, branch%ele(ix_start)%value(l$), .true., .true., orbit_start, orbit_end)

if (.not. particle_is_moving_forward(orbit_end)) then
  if (present(track_state)) track_state = ix_start
  return
endif

if (present(all_orb)) then
  call reallocate_coord(all_orb, branch%n_ele_max)
  all_orb(ix_start) = orbit_end
endif

! Track to ending element

ix_ele = modulo(ix_start, branch%n_ele_track) + 1
do
  if (ix_ele == ix_end) exit

  call track1 (orbit_end, branch%ele(ix_ele), branch%param, orbit_end)

  if (present(all_orb)) all_orb(ix_ele) = orbit_end
  if (.not. particle_is_moving_forward(orbit_end)) then
    if (present(track_state)) track_state = ix_ele
    return
  endif
  ix_ele = modulo(ix_ele, branch%n_ele_track) + 1
enddo

! Track to s_end

call twiss_and_track_intra_ele (branch%ele(ix_end), branch%param, 0.0_rp, s_end-branch%ele(ix_end)%s_start, &
                                                                      .true., .true., orbit_end, orbit_end)
if (.not. particle_is_moving_forward(orbit_end) .and. present(track_state)) track_state = ix_end

end subroutine
