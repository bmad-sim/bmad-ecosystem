program sr
! Example program for getting synch radiation 
!   power results for a lattice
!
! mjf 2007.02.08
!

  use bmad
  use sr_mod

  implicit none

  type (lat_struct) :: lat
  type (coord_struct), allocatable :: orb(:)
  type (wall_struct) :: inside, outside
  type (general_lat_param_struct) :: gen_params

  type (ele_power_struct), allocatable :: power(:)

  integer i, ix
                                
  character this_lat*80, line*80, temp*80
  character*80 lat_file

  real(rp) wall_offset, seg_len, end_s


! get lattice
  type '(a, $)', ' What lattice file? <cerl/trunk/model/lat.bmad> : '
  accept '(a)', line
  call string_trim( line, line, ix )
  if (ix /= 0) then
    lat_file = line
  else
    lat_file = 'cerl/trunk/model/lat.bmad'
  endif
  gen_params%lattice = lat_file

  call bmad_parser( lat_file, lat )


  allocate( orb(0:lat%n_ele_max) )
  allocate( power(lat%n_ele_max) )

! query user for parameters

  type '(a, $)', ' Current per beam <CR = 0.1 Amps>: '
  accept '(a)', line
  call string_trim( line, line, ix )
  if (ix /= 0) then
    read (line, *) gen_params%i_beam
  else
    gen_params%i_beam = 0.1   ! 100 mA / beam
  endif

  type '(a, $)', ' Vertical Emittance <CR = 10.0e-12>: '
  accept '(a)', line
  call string_trim( line, line, ix )
  if (ix /= 0) then
    read (line, *) gen_params%epsilon_y
  else
    gen_params%epsilon_y = 0.01e-9 
  endif

  type '(a, $)', ' Wall offset <CR = 4.48e-2 m>: '
  accept '(a)', line
  call string_trim( line, line, ix )
  if (ix /= 0) then
    read (line, *) wall_offset
  else
    wall_offset = 4.48e-2   
  endif

  type '(a, $)', ' Wall segment length <CR = 0.1 m>: '
  accept '(a)', line
  call string_trim( line, line, ix )
  if (ix /= 0) then
    read (line, *) seg_len
  else
    seg_len = 0.1 
  endif

! calculate twiss, closed orbits

  call twiss_at_start( lat )
  call closed_orbit_calc( lat, orb, 4 )
  call lat_make_mat6( lat, -1, orb )
  call twiss_at_start( lat )
  call twiss_propagate_all( lat )


! create a wall outline and break into segments

  end_s = lat%ele(lat%n_ele_track)%s

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
  inside%pt(0)%x = wall_offset * -1
  inside%pt(0)%name = 'INSIDE'
  inside%pt(0)%ix_pt = 0
  inside%n_pt_tot = 1
  inside%pt(1)%s = end_s
  inside%pt(1)%x = wall_offset * -1
  inside%pt(1)%name = 'INSIDE'
  inside%side = -1
  inside%pt(1)%ix_pt = 1


  call delete_overlapping_wall_points( outside )
  call delete_overlapping_wall_points( inside )

  call break_wall_into_segments( inside, seg_len )
  call break_wall_into_segments( outside, seg_len )

! calculate power densities

  call init_wall( outside )
  call init_wall( inside )

  call calculate_sr_power( lat, orb, +1, power, &
       inside, outside, gen_params )

! write out results

  call write_power_results( outside, lat, gen_params )
  call write_power_results( inside, lat, gen_params )

  open (unit = 1, file = 'element_power.dat', carriagecontrol = 'list')
  write (1, *) '  Ix  Name               |    S Position     |      Power (W)     |'
  write (1, *) '                         |   Start      End  | Radiated      Wall |'
  do i = 1, lat%n_ele_track
    if ( (lat%ele(i)%key == sbend$) .or. (lat%ele(i)%key == rbend$) .or. &
          (lat%ele(i)%key == wiggler$) ) then
      write (1, '(i4, 2x, a20, 2f10.3, 2f10.0)') i, lat%ele(i)%name, &
              lat%ele(i-1)%s, lat%ele(i)%s, &
              power(i)%radiated, power(i)%at_wall
    endif
  enddo

  close (unit = 1)
  type *, 'Written: element_power.dat'

  deallocate(orb)
  deallocate(power)

end program sr
