!+
! Subroutine tao_hook_init_read_lattice_info (file_name) 
!
! Hook routine to initialize lattice information.
! Note: use the logical s%com%init_read_lat_info to control 
! the standard lattice initialization.
!
! Input:
!   file_name -- Character(*): File containing lattice information.
!-

subroutine tao_hook_init_read_lattice_info (file_name)

use tao_interface

implicit none

character(*) file_name

!

end subroutine
