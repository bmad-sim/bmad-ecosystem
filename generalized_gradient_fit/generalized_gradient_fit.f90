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
              z_min, z_max, sym_x, sym_y, mode, lmdif_eps, printit, out_file, n_deriv_max_cos, n_deriv_max_sin, &
              n_planes_add, optimizer, x_pos_plot, y_pos_plot, core_weight, ele_anchor_pt, &
              Nx_min, Nx_max, Ny_min, Ny_max, Nz_min, Nz_max, del_grid, r0_grid, field_scale, length_scale

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
! Write ASCII version of field table from binary version.
case ('ascii')
  call read_field_table(field_file)
  call write_ascii_field_table(out_file)

! Write binaray version of field talbe from ASCII version.
case ('binary')
  call read_field_table(field_file)
  call write_binary_field_table(out_file)

! Fit field table and output GG.
case ('fit')
  call read_field_table(field_file)
  call fit_field()
  call write_gg_info()
  call write_gg_bmad()

! Create field table from GG and write it (ASCII)
case ('output_table')
  call read_gg(.true.)
  call write_ascii_field_table(out_file)

! Read GG and output truncated version.
case ('truncate_gg')
  call read_gg(.false.)
  call write_gg_bmad()
  
case default
  print *, 'I do not understand this mode: ', trim(mode)
end select

end program
