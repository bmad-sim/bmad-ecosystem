!+
! Subroutine tao_hook_lattice_calc (universe, lattice, orbit, used)
!
!  This hook is used to do custom lattice calculations. This routine is only called if
! the lattice needs to be recalculated. If this hook is used then set used to
! TRUE or the standard lattice calculation will overwrite whatever you do here!
! You have the option of only doing a custum lattice calculation for a
! particular universe. In that case, only set used to TRUE when that universe is
! passed.
!
!  Keep in mind that the standard Tao data types assume twiss_and_track
! is called. If you do custom tracking then somehow the lattice orbit and
! twiss parameters need to be updated or the standard TAO data types are liable to break! 
!
!  See tao/code/tao_lattice_calc_mod.f90 for how the standard lattice calculation is
! performed. All bookkeeping is performed in tao_lattice_calc so all you need to do
! here is the actual lattice calculation for whatever universe, lattice and
! orbit is passed to it. 
!
!  However, if coupled lattices are used then the injection needs to be taken 
! care of in here. Coupled lattice injection is complex so look through 
! tao_lattice_calc_mod.f90 thoroughly, expecially if using macroparticles!
! The routines in tao_lattice_calc_mod can be useful.
!
! Input:
!  universe   -- tao_universe_struct: universe to do calculation in
!  lattice    -- ring_struct: lattice to use (i.e. model, design, etc...)
!  orbit      -- coord_struct(0:): orbit structure (i.e. model, design, etc...)
!
! Output:
!  universe   -- tao_universe_struct: universe to do calculation in
!  lattice    -- ring_struct: lattice to use (i.e. model, design, etc...)
!  orbit      -- coord_struct(0:): orbit structure (i.e. model, design, etc...)
!  used       -- Logical: is this hook being used for this universe
!-

subroutine tao_hook_lattice_calc (universe, lattice, orbit, used)

use tao_mod
use tao_lattice_calc_mod
use tao_data_mod

implicit none

type (tao_universe_struct) :: universe
type (ring_struct) :: lattice
type (coord_struct) :: orbit(0:)
logical used

integer i

! For example here's a simple lattice calculation:
! call twiss_and_track (lattice, orbit)

! change this to TRUE if a custom lattice calculation is performed
  used = .false.

end subroutine tao_hook_lattice_calc
