subroutine tao_abort_command_file

use tao_struct

integer i
character(40) :: r_name = 'tao_abort_command_file'

!

do i = tao_com%cmd_file_level, 1, -1
  if (tao_com%cmd_file(i)%paused) return
  call out_io (s_warn$, r_name, 'ABORTING COMMAND FILE: ' // tao_com%cmd_file(i)%name)
  close (tao_com%cmd_file(i)%ix_unit)
  tao_com%cmd_file(i)%ix_unit = 0 
  tao_com%cmd_file_level = tao_com%cmd_file_level - 1
enddo

end subroutine
