program synrad

use bmad
use synrad_mod
use synrad_write_power_mod

implicit none

type (lat_struct) :: lat
type (coord_struct), allocatable :: orb(:)
type (walls_struct), target :: walls
type (wall_struct),pointer :: negative_x_wall, positive_x_wall
type (synrad_param_struct) :: sr_param
type (synrad_mode_struct) :: synrad_mode

type (ele_power_struct), allocatable :: fwd_power(:), back_power(:)

integer i, n, ix, n_arg, n_wall, ios, beam_direction
integer use_ele_ix

character(100) this_lat, line, temp
character(100) lat_file, in_file, wall_file
character(16) forward_beam, backward_beam

real(rp) end_s, wall_offset, s, x_in, x_out, seg_len

logical err_flag

namelist / synrad_params / sr_param, seg_len, wall_file, wall_offset, beam_direction, &
                           forward_beam, backward_beam, use_ele_ix

! set pointers
positive_x_wall => walls%positive_x_wall
negative_x_wall => walls%negative_x_wall

! get parameters

n_arg = cesr_iargc()
if (n_arg > 1) then
  print *, 'Usage: synrad <input_file>'
  print *, 'Default: <input_file> = synrad.in'
  stop
endif

in_file = 'synrad.in'
if (n_arg == 1) call cesr_getarg(1, in_file)

! Defaults

beam_direction = 0    
seg_len = 0.1
wall_file = 'NONE'
wall_offset = 0.045
sr_param%i_beam = 0.1
sr_param%epsilon_y = 10e-12
sr_param%n_slice = 20
forward_beam = "POSITRON"
backward_beam = "ELECTRON"
use_ele_ix = 0


! Read file

print *, 'Opening: ', trim(in_file)
open (1, file = in_file, status = "old")
read (1, nml = synrad_params)
close (1)

if (use_ele_ix .ne. 0) print *, "Only calculating power for element #",use_ele_ix

!

if (sr_param%lat_file(1:6) == 'xsif::') then
  call xsif_parser(sr_param%lat_file(7:), lat)
else
  call bmad_parser(sr_param%lat_file, lat)
endif

call reallocate_coord (orb, lat%n_ele_max)
allocate(back_power(lat%n_ele_max), fwd_power(lat%n_ele_max))

! create a wall outline and break into segments

end_s = lat%ele(lat%n_ele_track)%s
n = 2 * end_s / seg_len + 2
allocate (positive_x_wall%seg(n), negative_x_wall%seg(n))

negative_x_wall%side = negative_x$
positive_x_wall%side = positive_x$

if (wall_file == 'NONE') then

  allocate (positive_x_wall%pt(0:1), negative_x_wall%pt(0:1))

  positive_x_wall%pt(0)%s = 0.0
  positive_x_wall%pt(0)%ix_pt = 0
  positive_x_wall%pt(1)%s = end_s
  positive_x_wall%pt(1)%ix_pt = 1

  positive_x_wall%n_pt_tot = 1
  positive_x_wall%pt(:)%type = no_alley$
  positive_x_wall%pt(:)%name = 'POSITIVE_X_WALL'
  positive_x_wall%pt(:)%x = wall_offset
  positive_x_wall%pt(:)%phantom = .false.

  negative_x_wall%pt(0)%s = 0.0
  negative_x_wall%pt(0)%ix_pt = 0
  negative_x_wall%pt(1)%s = end_s
  negative_x_wall%pt(1)%ix_pt = 1

  negative_x_wall%n_pt_tot = 1
  negative_x_wall%pt(:)%type = no_alley$
  negative_x_wall%pt(:)%name = 'NEGATIVE_X_WALL'
  negative_x_wall%pt(:)%x = -wall_offset
  negative_x_wall%pt(:)%phantom = .false.

else
  open (1, file = wall_file, status = 'old')

  call skip_header (1, err_flag)

  ! count lines

  i = -1
  do     
    read (1, '(a)', iostat = ios) line
    if (ios < 0) exit
    if (ios > 0) then
      print *, 'READ ERROR IN FILE: ', trim(wall_file)
      call err_exit
    endif
    if (line == '') cycle
    i = i + 1
  enddo

  ! Allocate arrays read in data

  n_wall = i
  allocate (positive_x_wall%pt(0:n_wall), negative_x_wall%pt(0:n_wall))
  rewind (1)
  call skip_header (1, err_flag)
  i = -1
  do 
    read (1, '(a)', iostat = ios) line
    if (ios < 0) exit
    if (line == '') cycle
    i = i + 1
    read (line, *) s, x_in, x_out
    positive_x_wall%pt(i)%s = s
    positive_x_wall%pt(i)%x = x_out
    positive_x_wall%pt(i)%name = 'POSITIVE_X_WALL'
    positive_x_wall%pt(i)%type = no_alley$
    positive_x_wall%pt(i)%phantom = .false.
    positive_x_wall%pt(i)%ix_pt = i

    negative_x_wall%pt(i)%s = s
    negative_x_wall%pt(i)%x = x_in
    negative_x_wall%pt(i)%name = 'NEGATIVE_X_WALL'
    negative_x_wall%pt(i)%type = no_alley$
    negative_x_wall%pt(i)%phantom = .false.
    negative_x_wall%pt(i)%ix_pt = i
  enddo
  close (1)

  positive_x_wall%n_pt_tot = i
  negative_x_wall%n_pt_tot = i

  negative_x_wall%pt(i)%s  = lat%ele(lat%n_ele_track)%s
  positive_x_wall%pt(i)%s = lat%ele(lat%n_ele_track)%s
 
endif

do i = 0, negative_x_wall%n_pt_tot

enddo

!

call delete_overlapping_wall_points(positive_x_wall)
call delete_overlapping_wall_points(negative_x_wall)

call break_wall_into_segments(negative_x_wall, seg_len)
call break_wall_into_segments(positive_x_wall, seg_len)

! calculate power densities

call init_wall(positive_x_wall)
call init_wall(negative_x_wall)

call init_wall_ends(walls)

! Synch calculation

if (beam_direction < -1 .or. beam_direction > 1) then
  print *, 'INVALID BEAM DIRECITON:', beam_direction
  stop
endif

if (beam_direction == 0 .or. beam_direction == 1) then
  call synch_calc (1, forward_beam, fwd_power, synrad_mode) 
endif

if (beam_direction == 0 .or. beam_direction == -1) then
  call synch_calc (-1, backward_beam, back_power, synrad_mode) 
endif

! write out results
! set lat elements and twiss at wall segments

call write_power_results(positive_x_wall, lat, sr_param, use_ele_ix, synrad_mode)
call write_power_results(negative_x_wall, lat, sr_param, use_ele_ix, synrad_mode)

call write_results(positive_x_wall, lat, sr_param, use_ele_ix, synrad_mode)
call write_results(negative_x_wall, lat, sr_param, use_ele_ix, synrad_mode)

open (unit = 1, file = 'element_power.dat')

if (beam_direction == 0) then
  write (1, *) '  Ix  Name               |    S Position     |   Fwd_Power (W)    |   Back_Power (W)   |'
  write (1, *) '                         |   Start      End  | Radiated  Hit_Wall | Radiated  Hit_Wall |'
  do i = 1, lat%n_ele_max
    if (fwd_power(i)%radiated > 1 .or. back_power(i)%radiated > 1) then
      write (1, '(i4, 2x, a20, 2f10.3, 4f10.0)') i, lat%ele(i)%name, &
              lat%ele(i-1)%s, lat%ele(i)%s, &
              fwd_power(i)%radiated, fwd_power(i)%at_wall, &
              back_power(i)%radiated, back_power(i)%at_wall
    endif
  enddo

elseif (beam_direction == -1) then
  write (1, *) '  Ix  Name               |    S Position     |   Back_Power (W)   |'
  write (1, *) '                         |   Start      End  | Radiated  Hit_Wall |'
  do i = 1, lat%n_ele_max
    if (back_power(i)%radiated > 1) then
      write (1, '(i4, 2x, a20, 2f10.3, 2f10.0)') i, lat%ele(i)%name, &
              lat%ele(i-1)%s, lat%ele(i)%s, &
              back_power(i)%radiated, back_power(i)%at_wall
    endif
  enddo

elseif (beam_direction == 1) then
  write (1, *) '  Ix  Name               |    S Position     |    Fwd_Power (W)   |'
  write (1, *) '                         |   Start      End  | Radiated  Hit_Wall |'
  do i = 1, lat%n_ele_max
    if (fwd_power(i)%radiated > 1) then
      write (1, '(i4, 2x, a20, 2f10.3, 2f10.0)') i, lat%ele(i)%name, &
              lat%ele(i-1)%s, lat%ele(i)%s, &
              fwd_power(i)%radiated, fwd_power(i)%at_wall
    endif
  enddo
endif

close (unit = 1)
type *, 'Written: element_power.dat'

deallocate(orb)
deallocate(fwd_power, back_power)

!------------------------------------------------------------------------------
contains

subroutine synch_calc (direction, beam_type, power, synrad_mode)

type (synrad_mode_struct) :: synrad_mode
type (normal_modes_struct) :: mode
type (ele_power_struct), allocatable :: power(:)
integer direction
character(*) beam_type

!

if (beam_type == 'ELECTRON') then
  lat%param%particle = electron$
else if (beam_type == 'POSITRON') then
  lat%param%particle = positron$
endif

call twiss_and_track (lat, orb)
call calculate_synrad_power(lat, orb, direction, power, &
                              walls, sr_param, use_ele_ix)

call radiation_integrals (lat, orb, mode)

if (beam_type == 'ELECTRON') then
  synrad_mode%ele_mode = mode
  print *,'Electron horiz emittance:',mode%a%emittance
else if (beam_type == 'POSITRON') then
  synrad_mode%pos_mode = mode
  print *,'Positron horiz emittance:',mode%a%emittance
endif

end subroutine

end program
