program synrad

  use bmad
  use sr_mod

  implicit none

  type (lat_struct) :: lat
  type (coord_struct), allocatable :: orb(:)
  type (wall_struct) :: inside, outside
  type (synrad_param_struct) :: sr_param

  type (ele_power_struct), allocatable :: fwd_power(:), back_power(:)

  integer i, n, ix, n_arg, n_wall, ios, beam_direction

  character(100) this_lat, line, temp
  character(100) lat_file, in_file, wall_file
  character(16) forward_beam, backward_beam

  real(rp) end_s, wall_offset, s, x_in, x_out, seg_len

  logical err_flag

  namelist / synrad_params / sr_param, seg_len, wall_file, wall_offset, beam_direction, &
                             forward_beam, backward_beam

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
  forward_beam = "POSITRON"
  backward_beam = "ELECTRON"

! Read file

  print *, 'Opening: ', trim(in_file)
  open (1, file = in_file, status = "old")
  read (1, nml = synrad_params)
  close (1)

!

  call bmad_parser(sr_param%lat_file, lat)

  allocate(orb(0:lat%n_ele_max))
  allocate(back_power(lat%n_ele_max), fwd_power(lat%n_ele_max))

! create a wall outline and break into segments

  end_s = lat%ele(lat%n_ele_track)%s
  n = 2 * end_s / seg_len + 2
  allocate (outside%seg(n), inside%seg(n))

  inside%side  = inside$
  outside%side = outside$

  if (wall_file == 'NONE') then

    allocate (outside%pt(0:1), inside%pt(0:1))

    outside%pt(0)%s = 0.0
    outside%pt(0)%ix_pt = 0
    outside%pt(1)%s = end_s
    outside%pt(1)%ix_pt = 1

    outside%n_pt_tot = 1
    outside%pt(:)%type = no_alley$
    outside%pt(:)%name = 'OUTSIDE'
    outside%pt(:)%x = wall_offset
    outside%pt(:)%phantom = .false.

    inside%pt(0)%s = 0.0
    inside%pt(0)%ix_pt = 0
    inside%pt(1)%s = end_s
    inside%pt(1)%ix_pt = 1

    inside%n_pt_tot = 1
    inside%pt(:)%type = no_alley$
    inside%pt(:)%name = 'INSIDE'
    inside%pt(:)%x = -wall_offset
    inside%pt(:)%phantom = .false.

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
    allocate (outside%pt(0:n_wall), inside%pt(0:n_wall))
    rewind (1)
    call skip_header (1, err_flag)
    i = -1
    do 
      read (1, '(a)', iostat = ios) line
      if (ios < 0) exit
      if (line == '') cycle
      i = i + 1
      read (line, *) s, x_in, x_out
      outside%pt(i)%s = s
      outside%pt(i)%x = x_out
      outside%pt(i)%name = 'OUTSIDE'
      outside%pt(i)%type = no_alley$
      outside%pt(i)%phantom = .false.
      outside%pt(i)%ix_pt = i

      inside%pt(i)%s = s
      inside%pt(i)%x = x_in
      inside%pt(i)%name = 'INSIDE'
      inside%pt(i)%type = no_alley$
      inside%pt(i)%phantom = .false.
      inside%pt(i)%ix_pt = i
    enddo
    close (1)

    outside%n_pt_tot = i
    inside%n_pt_tot = i

    inside%pt(i)%s  = lat%ele(lat%n_ele_track)%s
    outside%pt(i)%s = lat%ele(lat%n_ele_track)%s
 
  endif

  do i = 0, inside%n_pt_tot

  enddo

!

  call delete_overlapping_wall_points(outside)
  call delete_overlapping_wall_points(inside)

  call break_wall_into_segments(inside, seg_len)
  call break_wall_into_segments(outside, seg_len)

! calculate power densities

  call init_wall(outside)
  call init_wall(inside)

! Synch calculation

  if (beam_direction == 0) then
    call synch_calc (1, forward_beam, fwd_power) 
    call synch_calc (-1, backward_beam, back_power) 
  elseif (beam_direction == -1) then
    call synch_calc (-1, backward_beam, back_power) 
  elseif (beam_direction == 1) then
    call synch_calc (1, forward_beam, fwd_power) 
  else 
    print *, 'INVALID BEAM DIRECITON:', beam_direction
  endif

! write out results

  call write_power_results(outside, lat, sr_param)
  call write_power_results(inside, lat, sr_param)

  open (unit = 1, file = 'element_power.dat', carriagecontrol = 'list')

  if (beam_direction == 0) then
    write (1, *) '  Ix  Name               |    S Position     |   Fwd_Power (W)    |   Back_Power (W)   |'
    write (1, *) '                         |   Start      End  | Radiated  Hit_Wall | Radiated  Hit_Wall |'
    do i = 1, lat%n_ele_max
      if (fwd_power(i)%radiated > 1 .or. back_power(i)%radiated > 1) then
        write (1, '(i4, 2x, a20, 2f10.3, 2f10.0)') i, lat%ele(i)%name, &
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

subroutine synch_calc (direction, beam_type, power)

  type (ele_power_struct), allocatable :: power(:)
  integer direction
  character(*) beam_type

!

  if (beam_type == 'ELECTRON') then
    lat%param%particle = electron$
  else if (beam_type == 'POSITRON') then
    lat%param%particle = positron$
  endif

  call twiss_at_start(lat)
  call closed_orbit_calc(lat, orb, 4)
  call lat_make_mat6(lat, -1, orb)
  call twiss_at_start(lat)
  call twiss_propagate_all(lat)

  call calculate_sr_power(lat, orb, direction, power, &
                                inside, outside, sr_param)

end subroutine

end program
