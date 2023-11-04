!+
! Subroutine tao_hook_branch_calc (u, tao_lat, branch)
!
! Subroutine to do custom calculations on a lattice branch.
!
! Input:
!   u       -- tao_universe_struct:
!   tao_lat -- tao_lattice_struct: 
!   branch  -- branch_struct: Lattice branch to work on.
!-

subroutine tao_hook_branch_calc (u, tao_lat, branch)

use tao_interface

implicit none

type (tao_universe_struct), target :: u
type (tao_lattice_struct), target :: tao_lat
type (branch_struct), target :: branch

! This dummy routine does nothing.

end subroutine
