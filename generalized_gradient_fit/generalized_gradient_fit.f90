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
              ix_z_min, ix_z_max, Bx_sym_x, Bx_sym_y, By_sym_x, By_sym_y, Bz_sym_x, Bz_sym_y, &
              mode, lmdif_eps, printit

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

select case (mode)
case ('binary')
  call read_field_table(field_file)
  call write_binary_field_table(field_file)

!

case ('fit')
  call read_field_table(field_file)
  call fit_field_table()

!

case default
  print *, 'I do not understand this mode: ', trim(mode)
end select

end program
