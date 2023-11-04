!+
! Subroutine tao_hook_lattice_calc (calc_ok)
!
! This hook is used to do custom lattice calculations. This routine is only called if
! the lattice needs to be recalculated. If this hook is used then set used to
! TRUE or the standard lattice calculation will overwrite whatever you do here!
! You have the option of only doing a custum lattice calculation for a
! particular universe. In that case, only set used to TRUE when that universe is
! passed.
!
! Tao_lattice_calc finds the lattice parameters and then calls
! tao_load_data_array for every element as the beam is being tracked. If you do
! custom tracking then be sure to also calculate the lattice parameters and
! load the data arrays or TAO is liable to break!
!
! See tao/code/tao_lattice_calc_mod.f90 for how the standard lattice
! calculation is performed. 
!
! Output:
!   s%u(i)%universe_recalc 
!              -- Logical: Set this to False to suppress tao_lattice_recalc.
!   calc_ok    -- Logical: Set False if there was an error in the 
!                   calculation like a particle was lost or a lat is unstable.
!-

subroutine tao_hook_lattice_calc (calc_ok)

use tao_interface

implicit none

logical calc_ok

! 

end subroutine tao_hook_lattice_calc
