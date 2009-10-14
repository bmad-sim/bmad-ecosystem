!+
! Subroutine tao_hook_init_read_lattice_info (is_set) 
!
! Hook routine to initialize lattice information
!
! Output:
!   is_set -- Logical: Set True if the lattice info has been initialized by this 
!                routine. False otherwise. This tells the calling routine whether 
!                to skip the standard lattice initialization.
!-

subroutine tao_hook_init_read_lattice_info (file_name, is_set)

implicit none

character(*) file_name
logical is_set

!

is_set = .false.

end subroutine
