!+
! Subroutine tao_hook_lattice_calc (universe, lattice, orbit)
!
!  This hook is used to do custom lattice calculations. If you need to do something 
! more complex then add it here. s%global%lattice_recalc is set to true whenever the 
! model lattice has changed.
!
!  The optional types passed are when not doing a typical lattice calculation
!  but needing to do the calculations for a specific lattice. This is used in
!  particular during initialization where the design lattice and orbit is calculated. 
!  So, make sure you have the logic to take care of this special case!
!
!  Keep in mind that the standard Tao data types assume twiss_and_track (...)
! is called. If you do custom tracking then somehow s%u(:)%model and
! s%u(:)%model_orb needs to be updated or the standard Tao data types are liable to break! 
! If s%global%lattice_recalc if FALSE after this subrotuine is
! complete then the normal TAO lattice calculation will be performed.
!
!  See tao/code/tao_lattice_calc.f90 for how the standard lattice calculation is
! performed. All bookkeeping is performed in that routine so all you need to do
! here is the actual lattice calculation.
!
!-

subroutine tao_hook_lattice_calc (universe, lattice, orbit)

use tao_mod
use tao_data_mod

implicit none

type (tao_universe_struct), optional :: universe
type (ring_struct), optional :: lattice
type (coord_struct), optional :: orbit(0:)

integer i

! if you do a custom lattice calculation here then make sure
! s%global%lattice_racalc is set to FALSE or the normal TAO lattice calculation
! will be performed

if (s%global%lattice_recalc) then
  !do your custom calculation here and then uncomment the following line:

! s%global%lattice_recalc = .false.

endif

end subroutine tao_hook_lattice_calc
