!+
! Subroutine tao_hook_init_lattice_post_parse (u)
!
! Subroutine to do custom post processing on the design lattice
!
! Input:
!   u   -- Tao_universe_struct: Universe holding the lattice
!
! Output:
!   u   -- Tao_universe_struct: Universe holding the lattice
!-

subroutine tao_hook_init_lattice_post_parse (u)

use tao_interface

implicit none

type (tao_universe_struct) u

character(*), parameter :: r_name = 'tao_hook_init_lattice_post_parse'
  
!
  
end subroutine 
