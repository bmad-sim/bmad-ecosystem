program bmad_to_scibmad

use bmad_routine_interface

implicit none

type (lat_struct) lat
character(200) bmad_name, scibmad_name
logical err_flag

!

call get_command_argument (1, bmad_name)
call get_command_argument (2, scibmad_name)

call bmad_parser(bmad_name, lat)

if (scibmad_name == '') then
  call file_suffixer(bmad_name, scibmad_name, '.jl', .true.)
endif

call write_lattice_in_scibmad(scibmad_name, lat, err_flag)
print *, 'Scibmad file: ' // trim(scibmad_name)

end program
