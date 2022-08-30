!+
! Subroutine tao_spin_new_polarization_calc (branch, tao_branch)
!
! Routine to calculate the spin equalibrium polarization in a ring along with the polarization rate and
! the depolarization rate due to emission of synchrotron radiation photons.
!
! From the Handbook of Accelerator Physics
!
! Input:
!   branch        -- branch_struct: Lattice branch to analyze.
!   tao_branch    -- tao_lattice_branch_struct: Contains %orbit
!
! Output:
!   tao_branch    -- tao_lattice_branch_struct: Calculated is:
!     %dn_dpz(:) 
!     %spin
!-

subroutine tao_spin_new_polarization_calc (branch, tao_branch)

use tao_data_and_eval_mod, dummy => tao_spin_polarization_calc
use radiation_mod
use ptc_interface_mod

implicit none

type (branch_struct), target :: branch
type (tao_lattice_branch_struct), target :: tao_branch

!

end subroutine tao_spin_new_polarization_calc
