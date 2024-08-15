!+
! Subroutine tao_abort_command_file()
!
! Routine to close all open command files.
!-

subroutine tao_abort_command_file()

use tao_struct

implicit none

integer i
character(40) :: r_name = 'tao_abort_command_file'

!

call tao_quiet_set('cmd-file-end')

do i = s%com%cmd_file_level, 1, -1
  if (s%com%cmd_file(i)%paused .or. s%global%single_step) return
  call out_io (s_warn$, r_name, 'ABORTING COMMAND FILE DUE TO DETECTED ERROR: ' // s%com%cmd_file(i)%full_name)
  call tao_close_command_file()
enddo

s%com%saved_cmd_line = ''

end subroutine
