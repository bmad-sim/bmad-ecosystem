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

call tao_quiet_set('off')

do i = s%com%cmd_file_level, 1, -1
  if (s%com%cmd_file(i)%paused .or. s%global%single_step) return
  call out_io (s_warn$, r_name, 'ABORTING COMMAND FILE: ' // s%com%cmd_file(i)%full_name)
  close (s%com%cmd_file(i)%ix_unit)
  s%com%cmd_file(i)%ix_unit = 0 
  s%com%cmd_file_level = s%com%cmd_file_level - 1
enddo

s%com%saved_cmd_line = ''

end subroutine
