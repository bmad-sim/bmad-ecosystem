!+
! Subroutine twiss_and_track_from_s_to_s (branch, orbit_start, s_end, orbit_end, ele_start, ele_end, 
!                                                                  err, compute_floor_coords, compute_twiss)
!
! Routine to track a particle from one location to another
!
! See also: 
!   twiss_and_track_at_s
!   track_from_s_to_s
!
! Note: orbit_start%ix_ele, orbit_start%location, only matter when orbit_start%s corresponds 
! to a position where there are zero length elements that affect the orbit like kickers or
! when the particle is at an element with edge fields like a bend.
!
! If s_end corresponds to the location of a zero length element or elements, the particle will
! stop just before any such elements.
!
! Note: Use element_at_s to properly initialize orbit_start.
!
! Input:
!   branch               -- branch_struct: Lattice branch to track through.
!   orbit_start          -- Coord_struct: Starting phase space coordinates at s_start.
!     %s                     -- Starting position.
!     %ix_ele                -- Starting element. 
!     %location              -- Location relative element.
!   s_end                -- real(rp): Ending position.
!   ele_start            -- Ele_struct, optional: Holds the starting parameters at s_start.
!   compute_floor_coords -- logical, optional: If present and True then the global "floor" coordinates will be 
!                             calculated and put in ele_end%floor.
!   compute_twiss        -- logical, optional: Default True. If False, to save a little time, do not 
!                             compute Twiss parameters.
!
! Output:
!   orbit_end      -- Coord_struct: End phase space coordinates. 
!   ele_end        -- Ele_struct, optional: Holds the ending Twiss parameters and the transfer matrix.
!                       If present then the ele_start argument must also be present.
!   err            -- Logical, optional: Set True if there is a problem like 
!                      the particle gets lost in tracking
!-   

subroutine twiss_and_track_from_s_to_s (branch, orbit_start, s_end, orbit_end, ele_start, ele_end, &
                                                                   err, compute_floor_coords, compute_twiss)

use bookkeeper_mod, dummy => twiss_and_track_from_s_to_s

implicit none

type (branch_struct), target :: branch
type (coord_struct) :: orbit_start, orbit_end
type (coord_struct) orbit
type (ele_struct), optional :: ele_start, ele_end
type (ele_struct) ele_here
type (ele_struct), pointer :: ele0, ele_track

real(rp) s_start, s_end
real(rp) s0, ds, s_true_end

integer ix_start, ix_end
integer ix_ele

logical, optional, intent(inout) :: err
logical, optional :: compute_floor_coords, compute_twiss
logical track_upstream_end, err_flag

character(*), parameter :: r_name = 'twiss_and_track_from_s_to_s'

! Easy case & error check

s_start = orbit_start%s
ix_end = element_at_s (branch, s_end, .true., err_flag, s_true_end)

if (s_start == s_true_end .and. branch%param%geometry == open$) then
  orbit_end = orbit_start
  if (present(ele_end)) ele_end = ele_start
  if (present(err)) err = .false.
  return
endif

if (s_true_end < s_start .and. branch%param%geometry == open$) then
  call out_io (s_error$, r_name, 'S_END < S_START WITH AN OPEN LATTICE.')
  if (present(err)) err = .true.
  return
endif

! 

orbit = orbit_start
if (orbit%location == downstream_end$) then
  orbit%ix_ele = modulo(orbit%ix_ele, branch%n_ele_track) + 1
  orbit%location = upstream_end$
  if (orbit%ix_ele == 1) then ! Have wrapped
    orbit%s = branch%ele(1)%s_start
    s_start = orbit%s
  endif
endif

ix_start = orbit%ix_ele
if (orbit%location == downstream_end$) then
  if (abs(orbit%s - branch%ele(ix_start)%s) > bmad_com%significant_length) then
    call out_io (s_error$, r_name, 'PARTICLE MARKED AT END OF ELEMENT BUT S-POSITION WRONG')
    if (present(err)) err = .false.
    return
  endif
  ix_start = modulo(ix_start, branch%n_ele_track) + 1
endif

if (present(err)) err = err_flag
if (err_flag) return

ele0 => branch%ele(ix_start)
s0 = branch%ele(ix_start)%s_start

track_upstream_end = (orbit%location == upstream_end$)

! Track within a single element case

if (s_true_end > s_start .and. ix_start == ix_end) then
  call twiss_and_track_intra_ele (ele0, branch%param, s_start-s0, s_true_end-s0, track_upstream_end, .true., &
                      orbit_start, orbit_end, ele_start, ele_end, err, compute_floor_coords, compute_twiss)
  return
endif

!--------------------------
! Track through multiple elements...
! First track to end of current element

call twiss_and_track_intra_ele (ele0, branch%param, s_start-s0, ele0%value(l$), track_upstream_end, .true., &
                      orbit_start, orbit_end, ele_start, ele_end, err_flag, compute_floor_coords, compute_twiss)

if (present(err)) err = err_flag
if (err_flag) return
if (.not. particle_is_moving_forward(orbit_end)) return

! Track to ending element

ix_ele = modulo(ix_start, branch%n_ele_track) + 1
do
  if (ix_ele == ix_end) exit
  ele_track => branch%ele(ix_ele)

  call track1 (orbit_end, ele_track, branch%param, orbit_end)
  if (.not. particle_is_moving_forward(orbit_end)) then
    err = .true.
    return
  endif

  ix_ele = modulo(ix_ele, branch%n_ele_track) + 1

  if (present(ele_end)) then
    if (logic_option(.false., compute_floor_coords)) call ele_geometry (ele_end%floor, ele_track, ele_end%floor)

    if (orbit_start%species == photon$) cycle

    ele_end%map_ref_orb_out = ele_track%map_ref_orb_out  ! Needed for dispersion calc.
    if (logic_option(.true., compute_twiss)) then
      call transfer_twiss (ele_end, ele_here)
      ele_here%mat6 = ele_end%mat6
      ele_here%vec0 = ele_end%vec0
      ele_end%key             = ele_track%key
      ele_end%value           = ele_track%value
      ele_end%mat6            = ele_track%mat6
      ele_end%map_ref_orb_in  = ele_track%map_ref_orb_in   ! Needed for dispersion calc.
      call twiss_propagate1 (ele_here, ele_end, err_flag)
      if (present(err)) err = err_flag
      if (err_flag) return
    endif
    ele_end%vec0 = matmul(ele_end%mat6, ele_here%vec0) + ele_track%vec0
    ele_end%mat6 = matmul(ele_end%mat6, ele_here%mat6)
  endif

enddo

! Track to s_true_end

if (present(ele_end) .and. orbit_start%species /= photon$) then
  ele_here%mat6 = ele_end%mat6
  ele_here%vec0 = ele_end%vec0
endif

ds = s_true_end-branch%ele(ix_end)%s_start
if (ds /= 0) then
  call twiss_and_track_intra_ele (branch%ele(ix_end), branch%param, 0.0_rp, ds, .true., .true., &
                                  orbit_end, orbit_end, ele_end, ele_end, err, compute_floor_coords)

  if (present(ele_end) .and. orbit_start%species /= photon$) then
    ele_end%vec0 = matmul(ele_end%mat6, ele_here%vec0) + ele_end%vec0
    ele_end%mat6 = matmul(ele_end%mat6, ele_here%mat6)
  endif
endif

end subroutine
