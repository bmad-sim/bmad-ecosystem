!+
! Example program showing how to do multi-turn tracking
!-

program multi_turn_tracking_example

use bmad

implicit none

type (lat_struct) lat
type (coord_struct), allocatable :: orbit(:)

integer i_turn, n_turn, track_state
logical err_flag

character(120) lat_file_name, input_file_name

namelist / input_params / lat_file_name, n_turn

! Read in parameters

if (command_argument_count() > 1) then
  print *, 'Usage: multi_turn_tracking_example <input_file_name>'
  stop
endif

input_file_name = 'multi_turn_tracking.init'
if (command_argument_count() == 1) call get_command_argument(1, input_file_name)

print *, 'Input file name: ', trim(input_file_name)
open (1, file = input_file_name)
read (1, nml = input_params)
close (1)

! Get lattice and track

call bmad_parser (lat_file_name, lat)
call reallocate_coord (orbit, lat)
call init_coord (orbit(0), lat%particle_start%vec, lat%ele(0), downstream_end$, lat%param%particle)

do i_turn = 1, n_turn
  call track_all (lat, orbit, 0, track_state, err_flag)
  if (track_state /= moving_forward$) then
    print *, 'Lost particle at element: ', trim(lat%ele(track_state)%name), track_state
    print *, '     At turn:', i_turn
    exit
  endif
  if (i_turn /= n_turn) orbit(0) = orbit(lat%n_ele_track)
enddo

end program
