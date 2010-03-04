module twiss_and_track_mod

use bmad_struct
use bmad_interface

!+
! Sburoutine twiss_and_track
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
!   orb_array(0:)       -- Coord_array_struct: Array of orbit arrays.
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
  if (branch%key == photon_branch$) cycle
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

end module
