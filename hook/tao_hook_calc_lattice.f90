!+
! Subroutine tao_hook_lattice_calc (used)
!
!  This hook is used to do the lattice calculations. Vanilla Tao just calls
! twiss_and_track when necessary. If you need to do something more complex then
! add it here. s%global%lattice_recalc is set to true whenever the model lattice
! has changed.
!
!  Keep in mind that the standard Tao data types assume twiss_and_track (...)
! is called here. If you remove this then be sure the model lattice and orbit
! is recalculated some other way!
!
!-

subroutine tao_hook_lattice_calc ()

use tao_mod
use tao_data_mod

implicit none

integer i

! Closed orbit and Twiss calculation.
! This can be slow for large lattices so only do it if the lattice changed.
if (s%global%lattice_recalc) then
  do i = 1, size(s%u)
    call twiss_and_track (s%u(i)%model, s%u(i)%model_orb)
  enddo
  s%global%lattice_recalc = .false.
endif

end subroutine tao_hook_lattice_calc
