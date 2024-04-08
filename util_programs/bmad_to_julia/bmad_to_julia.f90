program bmad_to_julia

use bmad_routine_interface

implicit none

type (lat_struct) lat
character(200) bmad_name, julia_name

!

call get_command_argument (1, bmad_name)
call get_command_argument (2, julia_name)
call bmad_parser(bmad_name, lat)
if (julia_name == '') julia_name = bmad_name

call write_lattice_in_julia(julia_name, lat, julia_name)
print *, 'Julia file: ' // trim(julia_name)

end program
