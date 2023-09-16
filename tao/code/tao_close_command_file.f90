!+
! Subroutine tao_close_command_file()
!
! Routine to close a command file
!-

subroutine tao_close_command_file()

use tao_interface, dummy => tao_close_command_file

implicit none

integer n

!

n = s%com%cmd_file_level

close (s%com%cmd_file(n)%ix_unit)
s%com%cmd_file(n)%ix_unit = 0 

if (s%com%cmd_file(n)%reset_at_end) then
  s%global%lattice_calc_on = s%com%cmd_file(n)%lattice_calc_save
  s%global%plot_on = s%com%cmd_file(n)%plot_save
endif

s%com%cmd_file_level = n - 1

end subroutine
