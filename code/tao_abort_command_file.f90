subroutine tao_abort_command_file

use tao_struct

integer i
character(40) :: r_name = 'tao_abort_command_file'

!

do i = s%com%cmd_file_level, 1, -1
  if (s%com%cmd_file(i)%paused) return
  call out_io (s_warn$, r_name, 'ABORTING COMMAND FILE: ' // s%com%cmd_file(i)%name)
  close (s%com%cmd_file(i)%ix_unit)
  s%com%cmd_file(i)%ix_unit = 0 
  s%com%cmd_file_level = s%com%cmd_file_level - 1
enddo

end subroutine
