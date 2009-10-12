!+
! Subroutine tao_hook_init_global (init_file, global)
!
! This routine is part of a collection of hook routines 
! used to bypass the use of an initialization file.
!-

subroutine tao_hook_init_global (init_file, global)

use tao_mod

use tao_struct

implicit none

type (tao_global_struct) global
character(*) init_file

!

end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine tao_hook_init_connected_uni (is_set)
!
! Hook routine to initialize connected universe information.
!
! Output:
!   is_set -- Logical: Set True if connected universe info has been initialized 
!                by this routine. False otherwise. This tells the calling routine 
!                whether to skip the standard connected universe initialization.
!-

subroutine tao_hook_init_connected_uni (is_set)

implicit none

logical is_set

!

is_set = .false.

end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine tao_hook_init_data (is_set) 
!
! Hook routine to initialize data.
!
! Output:
!   is_set -- Logical: Set True if data has been initialized by this routine.
!                False otherwise. This tells the calling routine whether to skip
!                the standard data initialization.
!-

subroutine tao_hook_init_data (is_set) 

use tao_mod

implicit none

logical is_set

!

is_set = .false.

end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine tao_hook_init_var (is_set) 
!
! Hook routine to initialize Tao variables.
!
! Output:
!   is_set -- Logical: Set True if variables have been initialized by this routine.
!                False otherwise. This tells the calling routine whether to skip
!                the standard variable initialization.
!-

subroutine tao_hook_init_var (is_set) 

use tao_mod

implicit none

logical is_set

!

is_set = .false.

end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine tao_hook_init_beam (is_set) 
!
! Hook routine to initialize the beam for beam tracking
!
! Output:
!   is_set -- Logical: Set True if the beam has been initialized by this routine.
!                False otherwise. This tells the calling routine whether to skip
!                the standard beam initialization.
!-

subroutine tao_hook_init_beam (is_set)

implicit none

logical is_set

!

is_set = .false.

end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
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

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine tao_hook_parse_command_args (is_set) 
!
! Hook routine to read the command arguments
!
! Output:
!   is_set -- Logical: Set True if the command arguments are parsed by this routine.
!                False otherwise. This tells the calling routine whether to skip
!                the standard command arg parsing.
!-

subroutine tao_hook_parse_command_args (is_set)

implicit none

logical is_set

!

is_set = .false.

end subroutine

