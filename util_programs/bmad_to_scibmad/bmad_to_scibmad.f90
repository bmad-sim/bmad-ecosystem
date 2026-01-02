program bmad_to_scibmad

use bmad_routine_interface

implicit none

type (lat_struct) lat
character(200) :: bmad_name = '', scibmad_name = '', arg
integer i, i_arg
logical err_flag, force

!

force = .false.

i_arg = 0

do i = 1, 3
  i_arg = i_arg + 1
  call get_command_argument (i_arg, arg)

  if (arg == '') exit

  if (arg == '-force') then
    force = .true.
  elseif (bmad_name == '') then
    bmad_name = arg
  else
    scibmad_name = arg
    exit
  endif
enddo
  
call bmad_parser(bmad_name, lat)

if (.not. force) then
  call twiss_and_track (lat, orbit, status, orb_start = lat%particle_start)
  if (status /= ok$) then
    call out_io (s_error$, r_name, 'PROBLEM TRACKING. NO OUTPUT GENERATED!', &
                                   'USE THE "-force" OPTION TO FORCE TRANSLATION.')
    stop
  endif
endif

if (scibmad_name == '') then
  call file_suffixer(bmad_name, scibmad_name, '.jl', .true.)
endif

call write_lattice_in_scibmad(scibmad_name, lat, err_flag)
print *, 'SciBmad file: ' // trim(scibmad_name)

end program
