!+
! Subroutine tao_hook_universe_calc_post_process (u, tao_lat, calc_ok)
!
! This hook is used to do custom lattice calculations at the end of the tao_lattice_calc routine.
!
! See tao/code/tao_lattice_calc.f90 for how the standard lattice 
! calculation is performed and how this routine is called.
!
! Also consider:
!   tao_hook_lattice_calc   ! Called at the beginning of tao_lattice_calc 
!   tao_hook_branch_calc
!
! Input:
!   u           -- tao_universe_struct: Universe under consideration.
!   tao_lat     -- tao_lattice_struct: Structure containing the lattice.
!   calc_ok     -- logical: Current state of the lattice calculation.
!
! Output:
!   u           -- tao_universe_struct: Universe under consideration.
!   tao_lat     -- tao_lattice_struct: Structure containing the lattice.
!   calc_ok     -- logical: Set False if there was an error in the calculation like a particle 
!                   was lost or a lattice is unstable. Note to programmers: Do not set True.
!-

subroutine tao_hook_universe_calc_post_process (u, tao_lat, calc_ok)

use tao_interface
implicit none

type (tao_universe_struct), target :: u
  type (tao_lattice_struct), target :: tao_lat
logical calc_ok

! 

end subroutine tao_hook_universe_calc_post_process
