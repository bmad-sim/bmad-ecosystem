!+
! Subroutine tao_abort_command_file(force_abort)
!
! Routine to close all open command files.
!
! Input:
!   force_abort     -- Logical, optional :: If present and True, ignore s%global%cmd_file_abort_on_error
!                        and abort any open command files.
!-

subroutine tao_abort_command_file(force_abort)

use tao_struct

implicit none

integer i
character(40) :: r_name = 'tao_abort_command_file'
logical, optional :: force_abort

!

if (.not. (logic_option(.false., force_abort) .or. s%global%cmd_file_abort_on_error)) return

do i = s%com%cmd_file_level, 1, -1
  if (s%com%cmd_file(i)%paused .or. s%global%single_step) return
  call out_io (s_error$, r_name, 'ABORTING COMMAND FILE DUE TO DETECTED ERROR: ' // s%com%cmd_file(i)%full_name)
  call tao_close_command_file()
enddo

s%com%saved_cmd_line = ''
call tao_quiet_set('cmd-file-end')

end subroutine
