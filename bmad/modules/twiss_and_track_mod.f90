module twiss_and_track_mod

use element_at_s_mod

!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!+
! Subroutine twiss_and_track
!
! Routine to calculate the twiss parameters, transport matrices and orbit.
!
! This routine is an overloaded name for:
!   twiss_and_track_branch (lat, orb, status, ix_branch, use_particle_start, print_err)
!   twiss_and_track_all (lat, orb_array, status, print_err)
!
! The essential difference between these two procedures is that
! twiss_and_track_branch only does the main branch while twiss_and_track_all
! does everything but the photon_fork elements.
!
! Note: This is not necessarily the fastest way to do things since this
! routine does the entire calculation from scratch.
!
! For a circular ring: If the RF is on, the computed orbit will be the 6D closed orbit.
! If the RF is off, the 4D transverse closed orbit using orbi(0)%vec(6) is computed.
!
! For an open lattice, the orbit will be computed using orb(0) as 
! starting conditions.
! 
! If there is a problem in a closed geometry branch, status argument settings are: in_stop_band$, 
! unstable$, non_symplectic$, in_stop_band$, unstable$, non_symplectic$, xfer_mat_clac_failure$, 
! twiss_propagate_failure$, or no_closed_orbit$. Note: in_stop_band$, unstable$, and non_symplectic$ 
! refer to the 1-turn matrix which is computed with closed lattices. A negative sign is used to 
! differentiate an error occuring in the first call to twiss_at_start from the second call to twiss_at_start.
!
! If there is a problem in an open geometry branch, status argument setting is -N where N is the element 
! where the particle was lost in tracking (negative numbers are used here to avoid confusion with ok$
! which is mapped to 1.
!
! Input:
!   lat                 -- lat_struct: lattice.
!     %param%geometry      -- Used to determine if lattice is open or closed.
!   orb(0:)             -- Coord_struct, allocatable: Orbit to be computed
!     orb(0)            -- Initial conditions to be used for an open geometry lattices.
!     orb(0)%vec(6)     -- For a closed lat: Energy at which the closed orbit is computed.
!   orb_array(0:)       -- Coord_array_struct, allocatable: Array of orbit arrays.
!     orb_array(0)%orbit(0) -- Used as the starting point for an open lattice.
!   ix_branch           -- Integer, optional: Branch to track.
!   use_particle_start  -- logical, optional: If True, use lat%particle_start instead of orb(0)
!                            as the initial coords for open geometry lattices. Default is False.
!   print_err           -- logical, optional: Default is True. If False, suppress error messages.
!
! Output:
!   lat                -- lat_struct: Lat with computed twiss parameters.
!     %param%stable          -- Set true or false.
!     %param%unstable_factor -- unstable growth rate (= 0 if stable)
!   orb(0:)            -- Coord_struct: Computed orbit.
!   orb_array(0:)      -- Coord_array_struct: Array of orbit arrays.
!   status             -- integer, optional: Set ok$ if everything is OK and set to something else otherwise.
!                           See above for more details.
!-

interface twiss_and_track
  module procedure twiss_and_track_branch
  module procedure twiss_and_track_all
end interface

private twiss_and_track_branch, twiss_and_track_all, twiss_and_track1

contains

!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!+
! Subroutine twiss_and_track_branch (lat, orb, status, ix_branch, use_particle_start, print_err)
!
! Subroutine to calculate the twiss parameters, transport matrices and orbit.
!
! This routine is overloaded by twiss_and_track.
! See twiss_and_track for more details.
!-

subroutine twiss_and_track_branch (lat, orb, status, ix_branch, use_particle_start, print_err)

implicit none

type (lat_struct) lat
type (coord_struct), allocatable :: orb(:)


integer, optional :: status, ix_branch
integer ib, status2

logical, optional :: use_particle_start, print_err

!

ib = integer_option(0, ix_branch)
call reallocate_coord (orb, lat%branch(ib)%n_ele_max)
if (logic_option(.false., use_particle_start)) &
                    call init_coord (orb(0), lat%particle_start, lat%branch(ib)%ele(0), downstream_end$)
call twiss_and_track1 (lat, orb, ib, status2, print_err)
if (present(status)) status = status2

end subroutine twiss_and_track_branch

!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!+
! Subroutine twiss_and_track_all (lat, orb_array, status, print_err)
!
! Subroutine to calculate the twiss parameters, transport matrices and orbit.
! Note: photon branches are currently ignored.
! 
! This routine is overloaded by twiss_and_track.
! See twiss_and_track for more details.
!-

subroutine twiss_and_track_all (lat, orb_array, status, print_err)

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (coord_array_struct), allocatable :: orb_array(:)

integer, optional :: status
integer i, status2
logical, optional :: print_err

!

call reallocate_coord (orb_array, lat)

do i = 0, ubound(lat%branch, 1)
  branch => lat%branch(i)
  if (branch%param%particle == photon$) cycle
  if (i /= 0 .and. branch%param%geometry == open$) then
    orb_array(i)%orbit(0) = orb_array(branch%ix_from_branch)%orbit(branch%ix_from_ele)
    call transfer_twiss (lat%branch(branch%ix_from_branch)%ele(branch%ix_from_ele), branch%ele(0))
  endif
  call twiss_and_track1 (lat, orb_array(i)%orbit, i, status2, print_err)
  if (present(status)) status = status2
  if (status2 /= ok$) return
enddo 

end subroutine twiss_and_track_all

!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!+
! Subroutine twiss_and_track1 (lat, ix_branch, orb, status, print_err)
!
! Subroutine to calculate the twiss parameters, transport matrices and orbit.
!
! This routine is private and used by:
!   twiss_and_track_branch
!   twiss_and_track_all
!-

subroutine twiss_and_track1 (lat, orb, ix_branch, status, print_err)

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele
type (coord_struct), allocatable :: orb(:)

integer i, ix_branch, status, stat

logical, optional :: print_err
logical err_flag, err

character(20) :: r_name = 'twiss_and_track1'

! We need to know the orbit first before we can compute the linear
! transfer matrices for the elements.
! However closed_orbit_calc needs some crude notion of the 1-turn transfer
! matrix in order for it to do the calculation.

branch => lat%branch(ix_branch)

! A match element may complicate things since in order to track correctly we may
! need to know the Twiss parameters. This situation is only allowed for open lattices.

if (branch%param%geometry == closed$) then
  call lat_make_mat6 (lat, -1, ix_branch = ix_branch)
  ! It can happen that the lattice is unstable about the zero orbit but not unstable about the closed orbit.
  ! So don't give up yet if twiss_at_start returns status /= ok$.
  call twiss_at_start (lat, status, ix_branch, .false.)

  if (rf_is_on(branch)) then
    call closed_orbit_calc (lat, orb, 6, 1, ix_branch, err, print_err)
  else
    call closed_orbit_calc (lat, orb, 4, 1, ix_branch, err, print_err)
  endif

  if (err) then
    status = no_closed_orbit$
    return
  endif

  call lat_make_mat6 (lat, -1, orb, ix_branch = ix_branch, err_flag = err)
  if (err) then
    status = xfer_mat_calc_failure$
    return
  endif

  call twiss_at_start (lat, status, branch%ix_branch, print_err)
  if (status /= ok$) then
    status = -status  ! To differentiate a failure here from the first call to twiss_at_start.
    return
  endif

  call twiss_propagate_all (lat, ix_branch, err)
  if (err) then
    status = twiss_propagate_failure$
    return
  endif

! geometry = open
else
  do i = 1, branch%n_ele_track
    ele => branch%ele(i)
    call make_mat6(ele, branch%param, orb(i-1), orb(i), err_flag)
    if (err_flag) then
      status = xfer_mat_calc_failure$
      return
    endif

    call twiss_propagate1(branch%ele(i-1), ele, err_flag)
    if (err_flag) then
      status = twiss_propagate_failure$
      return
    endif
  enddo

  status = ok$
endif

end subroutine twiss_and_track1

!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!+
! Subroutine twiss_and_track_at_s (lat, s, ele_at_s, orb, orb_at_s, ix_branch, err, use_last, compute_floor_coords)
! 
! Subroutine to return the twiss parameters and particle orbit at a 
! given longitudinal position.
!
! When calculating the Twiss parameters, this routine assumes 
! that the lattice elements already contain the Twiss parameters calculated
! for the ends of the elements. 
!
! Additionally, the orbit at the ends of the elements (contained in orb(:)) must be 
! precomputed when orb_at_s is present.
!
! Precomputation of Twiss and orbit at the element ends may be done with the twiss_and_track routine.
!
! See also:
!   twiss_and_track_from_s_to_s
!   twiss_and_track_intra_ele
!
! Input:
!   lat             -- lat_struct: Lattice.
!   s               -- Real(rp): Longitudinal position. If s is negative the
!                        the position is taken to be lat%param%total_length - s.
!   ele_at_s        -- Ele_struct, optional: If the use_last argument is True, ele_at_s is
!                        taken to contain valid Twiss parameters stored from a previous call
!                        to this routine.
!   orb(0:)         -- Coord_struct, optional: Orbit through the Lattice.
!   orb_at_s        -- Coord_struct, optional: If the use_last argument is True, orb_at_s is
!                        taken to contain the valid orbit stored from a previous call.
!   ix_branch       -- Integer, optional: Branch index, Default is 0 (main lattice).
!   use_last        -- logical, optional: If present and True, and if ele_at_s%s < s,
!                        then use ele_at_s and orb_at_s as the starting point for the present calculation.
!                        This can speed things up when the present s-position is in the middle 
!                        of a long complicated element and the tracking (EG: Runge-Kutta) is slow. 
!   compute_floor_coords
!                   -- logical, optional: If present and True then the global "floor" coordinates (without
!                        misalignments) will be calculated and put in ele_at_s%floor.
!
! Output:
!   ele_at_s  -- Ele_struct, optional: Element structure holding the Twiss parameters.
!                  if orb is not given then the Twiss parameters are calculated
!                  with respect to the zero orbit.
!   orb_at_s  -- Coord_struct, optional: Particle position at the position s.
!                  If orb_at_s is present then this routine assumes that orb is
!                  present.
!   err       -- Logical, optional: Set True if there is a problem in the 
!                 calculation, False otherwise.
!-

subroutine twiss_and_track_at_s (lat, s, ele_at_s, orb, orb_at_s, ix_branch, err, use_last, compute_floor_coords)

implicit none

type (lat_struct), target :: lat
type (ele_struct), optional :: ele_at_s
type (coord_struct), optional :: orb(0:)
type (coord_struct), optional :: orb_at_s
type (branch_struct), pointer :: branch

real(rp) s, s_use, s0, s_saved

integer, optional :: ix_branch
integer ie_at_s, i_branch

logical err_flag, use_l
logical :: init_needed = .true.
logical, optional :: err, use_last, compute_floor_coords

character(*), parameter :: r_name = 'twiss_and_track_at_s'

! If close enough to edge of element just use element info.

i_branch = integer_option(0, ix_branch)
branch => lat%branch(i_branch)

ie_at_s = element_at_s (lat, s, .false., i_branch, err_flag, s_use)
if (err_flag) then
  call out_io (s_error$, r_name, 'Bad S-position. Cannot compute Twiss parameters.')
  if (present(err)) err = .true. 
  return
endif

if (abs(s_use - branch%ele(ie_at_s)%s) < bmad_com%significant_length) then
  if (present(ele_at_s)) call transfer_ele(branch%ele(ie_at_s), ele_at_s, .true.)
  if (present(orb_at_s)) orb_at_s = orb(ie_at_s)
  if (present(err)) err = .false.
  return
endif

! Normal case where we need to partially track through an element.

s0 = branch%ele(ie_at_s)%s_start
use_l = logic_option(.false., use_last)
if (use_l) then
  if (present(ele_at_s)) then
    s_saved = ele_at_s%s
  else
    s_saved = orb_at_s%s
  endif
endif

if (use_l .and. s_saved < s_use .and. s_saved > s0) then
  call twiss_and_track_intra_ele (branch%ele(ie_at_s), branch%param, s_saved-s0, s_use-s0, &
                              .true., .true., orb_at_s, orb_at_s, ele_at_s, ele_at_s, err, compute_floor_coords)

else
  if (present(orb)) then
    call twiss_and_track_intra_ele (branch%ele(ie_at_s), branch%param, 0.0_rp, s_use-s0, &
            .true., .true., orb(ie_at_s-1), orb_at_s, branch%ele(ie_at_s-1), ele_at_s, err, compute_floor_coords)
  else
    call twiss_and_track_intra_ele (branch%ele(ie_at_s), branch%param, 0.0_rp, s_use-s0, &
            .true., .true., ele_start = branch%ele(ie_at_s-1), ele_end = ele_at_s, err = err, &
            compute_floor_coords = compute_floor_coords)
  endif
endif

end subroutine twiss_and_track_at_s

end module
