!+
! Subroutine tao_hook_init_global (init_file, global)
!
! This routine is part of a collection of hook routines 
! used to bypass the use of an initialization file.
!-

subroutine tao_hook_init_global (init_file, global)

use tao_interface, dummy => tao_hook_init_global

implicit none

type (tao_global_struct) global
character(*) init_file

!

end subroutine
