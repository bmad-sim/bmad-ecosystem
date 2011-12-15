module twiss_and_track_mod

use bmad_struct
use bmad_interface
use lat_geometry_mod

logical, private, save :: ele_is_saved, orb_is_saved

!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!+
! Subroutine twiss_and_track
!
! Routine to calculate the twiss parameters, transport matrices and orbit.
!
! This routine is an overloaded name for:
!   twiss_and_track_main (lat, orb, ok)
!   twiss_and_track_all (lat, orb_array, ok)
!
! The essential difference between these two procedures is that
! twiss_and_track_main only does the main branch while twiss_and_track_all
! does everything but the photon_branches.
!
! Note: This is not necessarily the fastest way to do things since this
! routine does the entire calculation from scratch.
!
! For a circular ring the computed orbit will be the closed orbit.
! For an open lattice, the orbit will be computed using orb(0) as 
! starting conditions.
!
! Modules Needed:
!   use bmad
!
! Input:
!   lat                  -- lat_struct: Input lat holding the lattice.
!     %param%lattice_type -- Used to determine if lattice is open or closed.
!                                = circular_lattice$ implies a closed lattice.
!                                all others imply an open lattice.
!   orb(0:)             -- Coord_struct, allocatable: Orbit to be computed
!     orb(0)            -- Initial conditions to be used for an open lat.
!     orb(0)%vec(6)     -- For a closed lat: Energy at which the closed orbit 
!                             is computed.
!   orb_array(0:)       -- Coord_array_struct, allocatable: Array of orbit arrays.
!     orb_array(0)%orb(0) -- Used as the starting point for a linear lattice.
!
! Output:
!   lat                -- lat_struct: Lat with computed twiss parameters.
!   orb(0:)            -- Coord_struct: Computed orbit.
!   orb_array(0:)      -- Coord_array_struct: Array of orbit arrays.
!   ok                 -- Logical, optional: Set True if everything ok. 
!                           False otherwise
!-

interface twiss_and_track
  module procedure twiss_and_track_main
  module procedure twiss_and_track_all
end interface

private twiss_and_track_main, twiss_and_track_all, twiss_and_track1

contains

!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!+
! Subroutine twiss_and_track_main (lat, orb, ok)
!
! Subroutine to calculate the twiss parameters, transport matrices and orbit.
!
! This routine is overloaded by twiss_and_track.
! See twiss_and_track for more details.
!-

subroutine twiss_and_track_main (lat, orb, ok)

implicit none

type (lat_struct) lat
type (coord_struct), allocatable :: orb(:)

logical, optional :: ok

!

call reallocate_coord (orb, lat%n_ele_max)
call twiss_and_track1 (lat, orb, 0, ok)

end subroutine twiss_and_track_main

!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!+
! Subroutine twiss_and_track_all (lat, orb_array, ok)
!
! Subroutine to calculate the twiss parameters, transport matrices and orbit.
! Note: photon branches are currently ignored.
! 
! This routine is overloaded by twiss_and_track.
! See twiss_and_track for more details.
!-

subroutine twiss_and_track_all (lat, orb_array, ok)

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (coord_array_struct), allocatable :: orb_array(:)

integer i

logical, optional :: ok

!

call reallocate_coord_array (orb_array, lat)

do i = 0, ubound(lat%branch, 1)
  branch => lat%branch(i)
  if (branch%param%particle == photon$) cycle
  if (i /= 0 .and. branch%param%lattice_type == linear_lattice$) then
    orb_array(i)%orb(0) = orb_array(branch%ix_from_branch)%orb(branch%ix_from_ele) 
    call transfer_twiss (lat%branch(branch%ix_from_branch)%ele(branch%ix_from_ele), branch%ele(0))
  endif
  call twiss_and_track1 (lat, orb_array(i)%orb, i, ok)
  if (.not. ok) return
enddo 

end subroutine twiss_and_track_all

!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!+
! Subroutine twiss_and_track1 (lat, ix_branch, orb, ok)
!
! Subroutine to calculate the twiss parameters, transport matrices and orbit.
!
! This routine is private and used by:
!   twiss_and_track_main
!   twiss_and_track_all
!-

subroutine twiss_and_track1 (lat, orb, ix_branch, ok)

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (coord_struct), allocatable :: orb(:)

integer i, ix_branch

logical, optional :: ok

character(20) :: r_name = 'twiss_and_track1'

! We need to know the orbit first before we can compute the linear
! transfer matrices for the elements.
! However closed_orbit_calc needs some crude notion of the 1-turn transfer
! matrix in order for it to do the calculation.

if (present(ok)) ok = .false.
branch => lat%branch(ix_branch)

! A match with match_end$ complicates things since in order to track correctly we
! need to know the Twiss parameters. This situation is only allowed for linear lattices.

if (branch%param%lattice_type == circular_lattice$) then
  if (ix_branch /= 0) then
    call out_io (s_fatal$, r_name, 'CIRCULAR NON-MAIN BRANCHES NOT YET IMPLEMENTED!')
    call err_exit
  endif

  call lat_make_mat6 (lat, -1, ix_branch = ix_branch)
  call twiss_at_start (lat)
  if (.not. bmad_status%ok) return
  call closed_orbit_calc (lat, orb, 4)
  if (.not. bmad_status%ok) return
else
  do i = 1, branch%n_ele_track
    if (branch%ele(i)%key == match$ .and. branch%ele(i)%value(match_end$) /= 0) then
      call twiss_propagate_all (lat, ix_branch)
      exit
    endif
  enddo
  call track_all (lat, orb, ix_branch)
endif

! Now we can compute the Twiss parameters.

call lat_make_mat6 (lat, -1, orb, ix_branch = ix_branch)

if (lat%param%lattice_type == circular_lattice$) then
  call twiss_at_start (lat)
  if (.not. bmad_status%ok) return
endif

call twiss_propagate_all (lat, ix_branch)
if (.not. bmad_status%ok) return

if (present(ok)) ok = .true.

end subroutine twiss_and_track1

!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!+
! Subroutine twiss_and_track_nullify_saved_data
!
! Routine to nullify the saved data which is used in the routine
!   twiss_and_track_at_s
!-

subroutine twiss_and_track_nullify_saved_data

ele_is_saved = .false.
orb_is_saved = .false.

end subroutine twiss_and_track_nullify_saved_data

!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!+
! Subroutine twiss_and_track_at_s (lat, s, ele, orb, orb_at_s, ix_branch, err, use_saved_data)
! 
! Subroutine to return the twiss parameters and particle orbit at a 
! given longitudinal position.
!
! Note: When calculating the Twiss parameters, this routine assumes 
! that the lattice elements already contain the Twiss parameters calculated
! for the ends of the elements. 
!
! Note: The use_saved argument is only effective if this routine is called sucessively with
! use_saved = True. Specifically:
!   1) If use_saved is not present, or is false, then the saved data is nullified.
!   2) If this routine is called with use_save = True and the data has been nullified then
!      the data will be initiated but the calculation is not affected.
! The twiss_and_track_nullify_saved_data routine can be used to nullify the saved data.
!
! See also:
!   twiss_and_track_from_s_to_s
!   twiss_and_track_intra_ele
!
! Modules Needed:
!   use bmad
!
! Input:
!   lat       -- lat_struct: Lattice.
!   s         -- Real(rp): Longitudinal position. If s is negative the
!                  the position is taken to be lat%param%total_length - s.
!   orb(0:)   -- Coord_struct, optional: Orbit through the Lattice.
!   ix_branch -- Integer, optional: Branch index, Default is 0 (main lattice).
!   use_saved_data 
!             -- Logical, optional: If present and True then use the result from the previous 
!                  calculation as a starting point for the present calculation.
!                  This can speed things up when the given s-position is in the middle 
!                  of a long complicated element and the tracking (EG: Runge Kutta) is slow. 
!
! Output:
!   ele      -- Ele_struct, optional: Element structure holding the Twiss parameters.
!                  if orb is not given then the Twiss parameters are calculated
!                  with respect to the zero orbit.
!   orb_at_s -- Coord_struct, optional: Particle position at the position s.
!             If orb_at_s is present then this routine assumes that orb is
!             present.
!   err      -- Logical, optional: Set True if there is a problem in the 
!                 calculation, False otherwise.
!-

subroutine twiss_and_track_at_s (lat, s, ele, orb, orb_at_s, ix_branch, err, use_saved_data)

implicit none

type (lat_struct), target :: lat
type (ele_struct), optional :: ele
type (coord_struct), optional :: orb(0:)
type (coord_struct), optional :: orb_at_s
type (branch_struct), pointer :: branch
type (coord_struct), save :: orb_saved
type (ele_struct), save :: ele_saved

real(rp) s, s_use, s0
real(rp), save :: s_saved

integer, optional :: ix_branch
integer i, i_branch

logical err_flag, saved_flag
logical :: init_needed = .true.
logical, optional :: err, use_saved_data

character(20), parameter :: r_name = 'twiss_and_track_at_s'

!

if (init_needed) then
  call init_ele(ele_saved)
  ele_is_saved = .false.
  orb_is_saved =  .false.
  init_needed = .false.
endif

saved_flag = logic_option(.false., use_saved_data)
if (.not. saved_flag) then
  ele_is_saved = .false.
  orb_is_saved = .false.
endif

! If close enough to edge of element just use element info.

i_branch = integer_option(0, ix_branch)
branch => lat%branch(i_branch)

!

call ele_at_s (lat, s, i, ix_branch, err_flag, s_use)
if (err_flag) then
  if (present(err)) err = .true. 
  return
endif

if (abs(s_use - branch%ele(i)%s) < bmad_com%significant_length) then
  if (present(ele)) ele = branch%ele(i)
  if (present(orb_at_s)) orb_at_s = orb(i)
  if (present(err)) err = .false.
  return
endif

! Normal case where we need to partially track through an element.

s0 = branch%ele(i-1)%s

if (saved_flag .and. s_saved < s_use .and. s_saved > s0) then
  if (present(orb)) then
    call twiss_and_track_intra_ele (branch%ele(i), branch%param, s_saved-s0, s_use-s0, &
                              .true., .true., orb_saved, orb_at_s, ele_saved, ele, err)
  else
    call twiss_and_track_intra_ele (branch%ele(i), branch%param, s_saved-s0, s_use-s0, &
                              .true., .true., ele_start = ele_saved, ele_end = ele, err = err)
  endif

else
  if (present(orb)) then
    call twiss_and_track_intra_ele (branch%ele(i), branch%param, 0.0_rp, s_use-s0, &
                              .true., .true., orb(i-1), orb_at_s, branch%ele(i-1), ele, err)
  else
    call twiss_and_track_intra_ele (branch%ele(i), branch%param, 0.0_rp, s_use-s0, &
                              .true., .true., ele_start = branch%ele(i-1), ele_end = ele, err = err)
  endif
endif

call ele_geometry (branch%ele(i-1)%floor, ele, ele%floor)

if (saved_flag) then
  s_saved = s_use
  if (present(orb)) then
    orb_saved = orb_at_s
    orb_is_saved = .true.
  endif
  if (present(ele)) then
    ele_saved = ele
    ele_is_saved = .true.
  endif
endif

end subroutine twiss_and_track_at_s

end module
