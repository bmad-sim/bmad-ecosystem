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
!   s%global%lun_command_file -- Integer: Logical unit number of the 
!                                   command file.
!-

subroutine tao_call_cmd (file_name, cmd_arg)

use tao_mod

implicit none


character(*) file_name
character(*), optional :: cmd_arg(:)
character(200) full_name
character(16) :: r_name = 'tao_call_cmd'

integer iu
integer, automatic :: lun_save(tao_com%cmd_file_level)

! Open the command file and store the unit number

tao_com%cmd_file_level = tao_com%cmd_file_level + 1

! reallocate lun_command_file array
if (tao_com%cmd_file_level .gt. 1) &
  lun_save = tao_com%lun_command_file(1:tao_com%cmd_file_level-1)

if (associated (tao_com%lun_command_file)) deallocate (tao_com%lun_command_file)
allocate (tao_com%lun_command_file(tao_com%cmd_file_level))

if (tao_com%cmd_file_level .gt. 1) &
  tao_com%lun_command_file(1:tao_com%cmd_file_level-1) = lun_save
  

! check for nested command files
!if (s%global%lun_command_file /= 0) then
!  call out_io (s_abort$, r_name, 'NESTED COMMAND FILES NOT ALLOWED!')
!  call err_exit
!endif

iu = lunget()
call tao_open_file ('TAO_COMMAND_DIR', file_name, iu, full_name)
if (iu == 0) then ! open failed
  tao_com%cmd_file_level = tao_com%cmd_file_level - 1
  return
endif

tao_com%lun_command_file(tao_com%cmd_file_level) = iu

! command arguments only valid for first level
if (tao_com%cmd_file_level .ne. 1 .and. present(cmd_arg)) then
  if (cmd_arg(1) .ne. ' ') then
    call out_io (s_warn$, r_name, &
                 "Command arguments only valid for first level command file")
    call out_io (s_blank$, r_name, &
                 "command arguments for this call will be ignored!")
  endif
elseif(present(cmd_arg)) then
  tao_com%cmd_arg = cmd_arg
else
  tao_com%cmd_arg = ' '
endif

end subroutine 
