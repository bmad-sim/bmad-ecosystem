!+
! Subroutine tao_hook_lattice_calc (used)
!
!  This hook is used to do custom lattice calculations. If you need to do something 
! more complex then add it here. s%global%lattice_recalc is set to true whenever the 
! model lattice has changed.
!
!  Keep in mind that the standard Tao data types assume twiss_and_track (...)
! is called. If you do custom tracking then somehow s%u(:)%model and
! s%u(:)%model_orb needs to be update or the standard Tao data types are liable to break! 
! If s%global%lattice_recalc if FALSE after this subrotuine is
! complete then the normal TAO lattice calculation will be performed.
!
!-

subroutine tao_hook_lattice_calc ()

use tao_mod
use tao_data_mod

implicit none

integer i

! if you do a custom lattice calculation here then make sure
! s%global%lattice_racalc is set to FALSE or the normal TAO lattice calculation
! will be performed

if (s%global%lattice_recalc) then
  !do your custom calculation here and then uncomment the following line:

! s%global%lattice_recalc = .false.

endif

end subroutine tao_hook_lattice_calc
