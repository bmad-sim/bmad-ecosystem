program write_foreign_test

use write_lattice_file_mod

type (lat_struct) lat

integer nargs
logical debug_mode

character(200) lat_file, out_file

!

bmad_com%auto_bookkeeper = .false.
global_com%exit_on_error = .false.
write_lat_debug_flag = .true. ! Limit output precision to prevent output shifts when running with different compilers
lat_file = 'write_foreign_test.bmad'

debug_mode = .false.
nargs = command_argument_count()

if (nargs > 0) then
  call get_command_argument(1, lat_file)
  print *, 'Using ', trim(lat_file)
  debug_mode = .true.
endif

if (.not. debug_mode) call output_direct(-1, .false., s_info$, s_error$)
call bmad_parser (lat_file, lat, .false.)

call file_suffixer(lat_file, out_file, 'mad8.now', .true.)
call write_lattice_in_foreign_format ('MAD-8', out_file, lat)

call file_suffixer(lat_file, out_file, 'madx.now', .true.)
call write_lattice_in_foreign_format ('MAD-X', out_file, lat)

!call file_suffixer(lat_file, out_file, 'sad.now', .true.)
!call write_lattice_in_foreign_format ('SAD', out_file, lat)

call file_suffixer(lat_file, out_file, 'lte.now', .true.)
call write_lattice_in_foreign_format ('ELEGANT', out_file, lat)

call file_suffixer(lat_file, out_file, 'julia.now', .true.)
call write_lattice_in_foreign_format ('JULIA', out_file, lat)

!call file_suffixer(lat_file, out_file, 'opal.now', .true.)
!call write_lattice_in_foreign_format ('OPAL-T', out_file, lat)

end program
