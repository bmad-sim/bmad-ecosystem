!+
! Subroutine twiss_and_track (lat, orb, ok)
!
! Subroutine to calculate the twiss parameters, transport matrices and orbit.
!
! Note: This is not necessarily the fastest way to do things since this
! routine does the entire calculation from scratch.
!
! For a circular ring the computed orbit will be the closed orbit.
! For an open lat the orbit will be computed using orb(0) as 
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
!
! Output:
!   lat                -- lat_struct: Lat with computed twiss parameters.
!     %mat6(6,6)       -- Transport matrices.
!   orb(0:)            -- Coord_struct: Computed orbit.
!   ok                 -- Logical, optional: Set True if everything ok. 
!                           False otherwise
!-

#include "CESR_platform.inc"

subroutine twiss_and_track (lat, orb, ok)

use bmad_struct
use bmad_interface, except_dummy => twiss_and_track

implicit none

type (lat_struct) lat
type (coord_struct), allocatable :: orb(:)

integer i

logical, optional :: ok

!

call reallocate_coord (orb, lat%n_ele_max)

! We need to know the orbit first before we can compute the linear
! transfer matrices for the elements.
! However closed_orbit_calc needs some crude notion of the 1-turn transfer
! matrix in order for it to do the calculation.

if (present(ok)) ok = .false.

! A match with match_end$ complicates things since in order to track correctly we
! need to know the Twiss parameters. This situation is only allowed for linear lattices.

if (lat%param%lattice_type == circular_lattice$) then
  call lat_make_mat6 (lat, -1)
  call twiss_at_start (lat)
  if (.not. bmad_status%ok) return
  call closed_orbit_calc (lat, orb, 4)
  if (.not. bmad_status%ok) return
else
  do i = 1, lat%n_ele_track
    if (lat%ele(i)%key == match$ .and. lat%ele(i)%value(match_end$) /= 0) then
      call twiss_propagate_all (lat)
      exit
    endif
  enddo
  call track_all (lat, orb)
endif

! now we can compute the Twiss parameters.

call lat_make_mat6 (lat, -1, orb)

if (lat%param%lattice_type == circular_lattice$) then
  call twiss_at_start (lat)
  if (.not. bmad_status%ok) return
endif

call twiss_propagate_all (lat)
if (.not. bmad_status%ok) return

if (present(ok)) ok = .true.

end subroutine
