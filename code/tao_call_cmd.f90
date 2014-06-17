!+
! Subroutine tao_call_cmd (file_name, cmd_arg)
!
! Routine to open a tao command file. If not found in the current director
! than the TAO_INIT_DIR will be searched.
! 
! Input:
!   file_name  -- Character(*): Name of the tao command file.
!   cmd_arg(9) -- Character(*), optional: Command file arguments.
!
! Output:
!   s%global%cmd_file -- Integer: Logical unit number of the 
!                                   command file.
!-

subroutine tao_call_cmd (file_name, cmd_arg)

use tao_mod, dummy => tao_call_cmd

implicit none


character(*) file_name
character(*), optional :: cmd_arg(:)
character(200) full_name
character(16) :: r_name = 'tao_call_cmd'

integer iu, nl
type (tao_command_file_struct) :: cmd_file(0:s%com%cmd_file_level)

! Open the command file and store the unit number

nl = s%com%cmd_file_level + 1
s%com%cmd_file_level = nl

! reallocate cmd_file array

if (nl > 1) cmd_file = s%com%cmd_file

if (allocated (s%com%cmd_file)) deallocate (s%com%cmd_file)
allocate (s%com%cmd_file(0:nl))

if (nl > 1) s%com%cmd_file(0:nl-1) = cmd_file
  
iu = lunget()
call tao_open_file (file_name, iu, full_name, s_error$)
if (iu == 0) then ! open failed
  s%com%cmd_file_level = nl - 1
  return
endif

s%com%cmd_file(nl)%ix_unit = iu
s%com%cmd_file(nl)%name = file_name

! Save command arguments.

if(present(cmd_arg)) then
  s%com%cmd_file(nl)%cmd_arg = cmd_arg
else
  s%com%cmd_file(nl)%cmd_arg = ' '
endif

end subroutine 
