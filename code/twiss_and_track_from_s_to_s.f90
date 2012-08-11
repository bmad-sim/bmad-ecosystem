!+
! Subroutine twiss_and_track_from_s_to_s (lat, s_start, s_end, track_entrance, 
!                track_exit, orbit_start, orbit_end, ele_start, ele_end, ix_branch, err)
!
! Routine to track a particle within an element.
!
! See also: 
!   twiss_and_track_at_s
!   track_from_s_to_s
!
! The track_entrance and track_exit arguments determine whether entrance and exit effects are 
! included. These effects are bend edge focusing and aperture checks. 
!
! Modules needed:
!   use bmad
!
! Input:
!   lat            -- lat_struct: Lattice to track through
!   s_start        -- Real(rp): Starting s-position.
!   s_end          -- Real(rp): Ending s-position. If <= s_start then will wrap
!   track_entrance -- Logical: If True then entrance effects are included in the tracking.
!                       But only if l_start = 0.
!   track_exit     -- Logical: If True then exit effects are included in the tracking but 
!                       only if s_end = ele%s (within bmad_com%significant_length tol).
!   orbit_start    -- Coord_struct, optional: Starting phase space coordinates at l_start.
!   ele_start      -- Ele_struct, optional: Holds the starting Twiss parameters at l_start.
!   ix_branch      -- Integer, optional: Lattice branch index. Default is 0 (main branch).
!
! Output:
!   orbit_end      -- Coord_struct, optional: End phase space coordinates. 
!                       If present then the orbit_start argument must also be present.
!   ele_end        -- Ele_struct, optional: Holds the ending Twiss parameters.
!                       If present then the ele_start argument must also be present.
!   err            -- Logical, optional: Set True if there is a problem like 
!                      the particle gets lost in tracking
!-   

subroutine twiss_and_track_from_s_to_s (lat, s_start, s_end, track_entrance, &
                       track_exit, orbit_start, orbit_end, ele_start, ele_end, ix_branch, err)

use bookkeeper_mod

implicit none

type (coord_struct), optional :: orbit_start, orbit_end
type (ele_struct), optional :: ele_start, ele_end
type (ele_struct) ele_here
type (lat_struct), target :: lat
type (ele_struct), pointer :: ele0, ele_track
type (branch_struct), pointer :: branch

real(rp) s_start, s_end
real(rp) s0

integer, optional :: ix_branch
integer ix_start, ix_end
integer ix_ele

logical track_entrance, track_exit, error
logical, optional :: err

character(40), parameter :: r_name = 'twiss_and_track_from_s_to_s'

! Easy case & error check

branch => lat%branch(integer_option(0, ix_branch))

if (s_start == s_end .and. branch%param%lattice_type == linear_lattice$) then
  if(present(orbit_end)) orbit_end = orbit_start
  if (present(ele_end)) ele_end = ele_start
  return
endif

if (s_end < s_start .and. branch%param%lattice_type == linear_lattice$) then
  call out_io (s_abort$, r_name, 'S_END < S_START WITH A LINEAR LATTICE.')
  if (bmad_status%exit_on_error) call err_exit
endif

! Find elements corresponding to s_start and s_stop

ix_start = element_at_s (lat, s_start, .true., ix_branch)
if (branch%ele(ix_start)%s == s_start) then
  ix_start = modulo (ix_start, branch%n_ele_track) + 1
endif
ele0 => branch%ele(ix_start)
s0 = branch%ele(ix_start-1)%s

ix_end = element_at_s (lat, s_end, .true., ix_branch)

! Track within a single element case

if (s_end > s_start .and. ix_start == ix_end) then
  call twiss_and_track_intra_ele (ele0, branch%param, s_start-s0, s_end-s0, &
                      track_entrance, track_exit, orbit_start, orbit_end, ele_start, ele_end, err)
  return
endif

! Track through multiple elements...
! First track to end of current element

call twiss_and_track_intra_ele (ele0, branch%param, s_start-s0, ele0%value(l$), &
                      track_entrance, .true., orbit_start, orbit_end, ele_start, ele_end, error)
if (present(err)) err = error
if (error) return
if (present(orbit_end)) then
  if (.not. particle_is_moving_forward(orbit_end, branch%param%particle)) return
endif

! Track to ending element

ix_ele = ix_start + 1
do
  if (ix_ele == ix_end) exit
  ele_track => branch%ele(ix_ele)

  if (present(orbit_end)) then
    call track1 (orbit_end, ele_track, branch%param, orbit_end)
    if (.not. particle_is_moving_forward(orbit_end, branch%param%particle)) then
      err = .true.
      return
    endif
  endif

  if (present(ele_end)) then
    call transfer_twiss (ele_end, ele_here)
    ele_here%mat6 = ele_end%mat6
    ele_here%vec0 = ele_end%vec0
    ele_end%mat6 = ele_track%mat6
    ele_end%map_ref_orb_in  = ele_track%map_ref_orb_in   ! Needed for dispersion calc.
    ele_end%map_ref_orb_out = ele_track%map_ref_orb_out  ! Needed for dispersion calc.
    call twiss_propagate1 (ele_here, ele_end, error)
    if (present(err)) err = error
    if (error) return
    ele_end%vec0 = matmul(ele_end%mat6, ele_here%vec0) + ele_track%vec0
    ele_end%mat6 = matmul(ele_end%mat6, ele_here%mat6)
  endif

  ix_ele = modulo(ix_ele, branch%n_ele_track) + 1
enddo

! Track to s_end

if (present(ele_end)) then
  ele_here%mat6 = ele_end%mat6
  ele_here%vec0 = ele_end%vec0
endif

call twiss_and_track_intra_ele (branch%ele(ix_end), branch%param, &
          0.0_rp, s_end-branch%ele(ix_end-1)%s, .true., track_exit, orbit_end, orbit_end, &
          ele_end, ele_end)

if (present(ele_end)) then
  ele_end%vec0 = matmul(ele_end%mat6, ele_here%vec0) + ele_end%vec0
  ele_end%mat6 = matmul(ele_end%mat6, ele_here%mat6)
endif

end subroutine
