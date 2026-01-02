program write_foreign_test

use write_lattice_file_mod
use twiss_and_track_mod

type (lat_struct) lat
type (coord_struct), allocatable :: orbit(:)

integer nargs, status
logical debug_mode

character(200) lat_file, out_file

!

bmad_com%auto_bookkeeper = .false.
global_com%exit_on_error = .false.
lat_file = 'write_foreign_test.bmad'

debug_mode = .false.
nargs = command_argument_count()

if (nargs > 0) then
  call get_command_argument(1, lat_file)
  print *, 'Using ', trim(lat_file)
  debug_mode = .true.
endif

if (.not. debug_mode) call output_direct(-1, .false., s_info$, s_error$)
call bmad_parser (lat_file, lat)
call twiss_and_track (lat, orbit, status, orb_start = lat%particle_start)

write_lat_debug_flag = .true. ! Limit output precision to prevent output shifts when running with different compilers
call write_lattice_in_foreign_format ('MAD-8', 'mad8.now', lat)
call write_lattice_in_foreign_format ('MAD-X', 'madx.now', lat)
call write_lattice_in_foreign_format ('SAD', 'sad.now', lat)
call write_lattice_in_foreign_format ('ELEGANT', 'lte.now', lat)
call write_lattice_in_foreign_format ('SCIBMAD', 'scibmad.now', lat)

! OPAL testing needs some work. Specifically there should be a separate lattice for testing.
!call file_suffixer(lat_file, out_file, 'opal.now', .true.)
!call write_lattice_in_foreign_format ('OPAL-T', out_file, lat)

end program
