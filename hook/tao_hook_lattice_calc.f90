!+
! Subroutine tao_hook_lattice_calc (universe, tao_lat, used, calc_ok)
!
!  This hook is used to do custom lattice calculations. This routine is only called if
! the lattice needs to be recalculated. If this hook is used then set used to
! TRUE or the standard lattice calculation will overwrite whatever you do here!
! You have the option of only doing a custum lattice calculation for a
! particular universe. In that case, only set used to TRUE when that universe is
! passed.
!
!  Tao_lattice_calc finds the lattice parameters and then calls
! tao_load_data_array for every element as the beam is being tracked. If you do
! custom tracking then be sure to also calculate the lattice parameters and
! load the data arrays or TAO is liable to break!
!
!  See tao/code/tao_lattice_calc_mod.f90 for how the standard lattice
! calculation is performed. 
!
! Input:
!   universe   -- tao_universe_struct: universe to do calculation in
!   tao_lat    -- Tao_lattice_struct: Lattice and orbit to use (i.e. model, design, etc...)
!
! Output:
!   universe   -- tao_universe_struct: universe to do calculation in
!   tao_lat    -- Tao_lattice_struct: Lattice and orbit to use (i.e. model, design, etc...)
!   used       -- Logical: is this hook being used for this universe
!   calc_ok    -- Logical: Set False if there was an error in the 
!                   calculation like a particle was lost or a ring is unstable.
!-

subroutine tao_hook_lattice_calc (universe, tao_lat, used, calc_ok)

use tao_mod
use tao_lattice_calc_mod
use tao_data_mod

implicit none

type (tao_universe_struct) :: universe
type (tao_lattice_struct) :: tao_lat
logical used, calc_ok

integer i

! For example here's a simple lattice calculation:
! call twiss_and_track (lattice, orbit)

! change this to TRUE if a custom lattice calculation is performed
  used = .false.

end subroutine tao_hook_lattice_calc
