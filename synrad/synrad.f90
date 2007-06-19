program synrad

  use bmad
  use sr_mod

  implicit none

  type (lat_struct) :: lat
  type (coord_struct), allocatable :: orb(:)
  type (wall_struct) :: inside, outside
  type (synrad_param_struct) :: sr_param

  type (ele_power_struct), allocatable :: power(:)

  integer i, n, ix, n_arg, n_wall_pt_max

  character(80) this_lat, line, temp
  character(80) lat_file, in_file, wall_file

  real(rp) end_s, wall_offset, s, x_in, x_out, seg_len

  logical err_flag

  namelist / synrad_params / sr_param, seg_len, wall_file
  namelist / wall_params / n_wall_pt_max

! get parameters

  n_arg = cesr_iargc()
  if (n_arg > 1) then
    print *, 'Usage: synrad <input_file>'
    print *, 'Default: <input_file> = synrad.in'
    stop
  endif

  in_file = 'synrad.in'
  if (n_arg == 1) call cesr_getarg(1, in_file)

  print *, 'Opening: ', trim(in_file)
  open (1, file = in_file, status = "old")
  read (1, nml = synrad_params)
  close (1)

!

  call bmad_parser(sr_param%lat_file, lat)

  allocate(orb(0:lat%n_ele_max))
  allocate(power(lat%n_ele_max))

! calculate twiss, closed orbits

  call twiss_at_start(lat)
  call closed_orbit_calc(lat, orb, 4)
  call lat_make_mat6(lat, -1, orb)
  call twiss_at_start(lat)
  call twiss_propagate_all(lat)

! create a wall outline and break into segments

  end_s = lat%ele(lat%n_ele_track)%s
  n = 2 * end_s / seg_len + 2
  allocate (outside%seg(n), inside%seg(n))

  if (wall_file == 'NONE') then

    allocate (outside%pt(0:1), inside%pt(0:1))

    outside%pt(0)%s = 0.0
    outside%pt(0)%x = wall_offset
    outside%pt(0)%name = 'OUTSIDE'
    outside%pt(0)%ix_pt = 0
    outside%n_pt_tot = 1
    outside%pt(1)%s = end_s
    outside%pt(1)%x = wall_offset
    outside%pt(1)%name = 'OUTSIDE'
    outside%side = 1
    outside%pt(1)%ix_pt = 1

    inside%pt(0)%s = 0.0
    inside%pt(0)%x = -wall_offset
    inside%pt(0)%name = 'INSIDE'
    inside%pt(0)%ix_pt = 0
    inside%n_pt_tot = 1
    inside%pt(1)%s = end_s
    inside%pt(1)%x = -wall_offset
    inside%pt(1)%name = 'INSIDE'
    inside%side = -1
    inside%pt(1)%ix_pt = 1

  else
    open (1, file = wall_file, status = 'old')
    read (1, nml = wall_params)
    allocate (outside%pt(0:n_wall_pt_max), inside%pt(0:n_wall_pt_max))
    outside%n_pt_tot = n_wall_pt_max
    inside%n_pt_tot = n_wall_pt_max
    call skip_header (1, err_flag)
    do  i = 0, n_wall_pt_max
      read (1, '(a)') line
      read (line, *) ix, s, x_in, x_out
      if (ix /= i) then
        print *, 'ERROR: IN WALL FILE: ', trim(wall_file)
        print *, '       WALL INDEX NOT IN ORDER:', I
        call err_exit
      endif
      outside%pt(ix)%s = s
      outside%pt(ix)%x = x_out
      outside%pt(ix)%name = 'OUTSIDE'
      outside%pt(ix)%ix_pt = ix
      inside%pt(ix)%s = s
      inside%pt(ix)%x = x_in
      inside%pt(ix)%name = 'INSIDE'
      inside%pt(ix)%ix_pt = ix
    enddo
    close (1)
  endif

  call delete_overlapping_wall_points(outside)
  call delete_overlapping_wall_points(inside)

  call break_wall_into_segments(inside, seg_len)
  call break_wall_into_segments(outside, seg_len)

! calculate power densities

  call init_wall(outside)
  call init_wall(inside)

  call calculate_sr_power(lat, orb, +1, power, &
       inside, outside, sr_param)

! write out results

  call write_power_results(outside, lat, sr_param)
  call write_power_results(inside, lat, sr_param)

  open (unit = 1, file = 'element_power.dat', carriagecontrol = 'list')
  write (1, *) '  Ix  Name               |    S Position     |      Power (W)     |'
  write (1, *) '                         |   Start      End  | Radiated      Wall |'
  do i = 1, lat%n_ele_track
    if ((lat%ele(i)%key == sbend$) .or. (lat%ele(i)%key == rbend$) .or. &
          (lat%ele(i)%key == wiggler$)) then
      write (1, '(i4, 2x, a20, 2f10.3, 2f10.0)') i, lat%ele(i)%name, &
              lat%ele(i-1)%s, lat%ele(i)%s, &
              power(i)%radiated, power(i)%at_wall
    endif
  enddo

  close (unit = 1)
  type *, 'Written: element_power.dat'

  deallocate(orb)
  deallocate(power)

end program
