!+
! Program generalized_gradient_fit
!
! Program to fit generalized gradiants to a field table.
!-

program generalized_gradient_fit

use gg_fit_mod
implicit none

integer :: n_header_lines = 0

character(100) file_name
character(16) mode

namelist /params / field_file, every_n_th_plane, n_deriv_max, m_cos, m_sin, n_cycles, &
              z_min, z_max, sym_x, sym_y, mode, lmdif_eps, printit, out_file, &
              n_planes_add, n_deriv_keep, optimizer, x_pos_plot, y_pos_plot

!

file_name = 'gg_fit.init'
if (command_argument_count() > 0) then
  call get_command_argument(1, file_name)
endif

print *, 'Input file: ', trim(file_name)
open (1, file = file_name, status = 'old')
read (1, nml = params)
close (1)

!

call read_field_table(field_file)

!

select case (mode)
case ('ascii')
  call write_ascii_field_table(field_file)

case ('binary')
  call write_binary_field_table(field_file)

case ('fit')
  call fit_field()
  call write_gg()

case default
  print *, 'I do not understand this mode: ', trim(mode)
end select

end program
