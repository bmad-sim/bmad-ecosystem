!+
! Subroutine tao_call_cmd (file_name, cmd_arg)
!
! Routine to open a tao command file. If not found in the current director
! than the TAO_INIT_DIR will be searched.
! 
! Input:
!   file_name  -- Name of the tao command file.
!   cmd_arg(9) -- Command file arguments.
!
! Output:
!   s%global%lun_command_file -- Integer: Logical unit number of the 
!                                   command file.
!-

subroutine tao_call_cmd (file_name, cmd_arg)

use tao_mod

implicit none


character(*) file_name
character(*) cmd_arg(:)
character(200) full_name

integer iu

! Open the command file and store the unit number

iu = lunget()
call tao_open_file ('TAO_COMMAND_DIR', file_name, iu, full_name)
if (iu == 0) return
s%global%lun_command_file = iu

tao_com%cmd_arg = cmd_arg

end subroutine 
