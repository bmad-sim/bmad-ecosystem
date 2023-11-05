!+
! Program tao_program
!
! Program wrapper
!-

program tao_program

use tao_interface

implicit none

procedure(tao_hook_command_def) :: tao_hook_command

! Call tao_top_level and we are done.

tao_hook_command_ptr => tao_hook_command

call tao_top_level()

end program
