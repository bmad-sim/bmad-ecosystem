subroutine do_synrad (walls, u, ring, gen_params, window)

  use sr_mod
  use cesrv_struct
  use cesrv_interface
  use synrad_window_mod

  implicit none

  type (walls_struct) walls
  type (universe_struct), target :: u
  type (lat_struct) ring, p_ring, e_ring
  type (coord_struct), allocatable :: e_orb(:), p_orb(:)
  type (synrad_param_struct) gen_params
  type (ele_power_struct), allocatable :: p_power(:), e_power(:)
  type (crotch_window_struct) window(:)
  type (normal_modes_struct) modes

  integer ix
  character*80 line

  allocate( e_orb(0:ring%n_ele_max) )
  allocate(  p_orb(0:ring%n_ele_max) )
  allocate(  e_power(ring%n_ele_max) )
  allocate(  p_power(ring%n_ele_max) )

  call get_input_string ('Current per beam <CR = 0.2 Amps>:', line)
  call string_trim (line, line, ix)
  if (ix /= 0) then
    read (line, *)gen_params%i_beam
  else
    gen_params%i_beam = 0.2   ! 100 mA / beam
  endif

  call radiation_integrals( ring, u%orb, modes )

  gen_params%epsilon_y = max( (modes%a%emittance * .02), modes%b%emittance )
  print *, "a%emit: ", modes%a%emittance, " b%emit: ", modes%b%emittance
  print *, 'Default emittance is: ',gen_params%epsilon_y
  call get_input_string ('Vertical Emittance <CR = DEFAULT>:', line)
  call string_trim (line, line, ix)
  if (ix /= 0) then
    read (line, *) gen_params%epsilon_y
    print *, ' Vertical Emittance set to: ',gen_params%epsilon_y
  endif

! initialize the window ray counters

  call window_ray_reset (window)

! calculate twiss, closed orbits

  if (ring%param%particle == electron$) then
    e_ring = ring
    e_ring%param%particle = electron$
    call closed_orbit_calc (e_ring, e_orb, 4)
    call lat_make_mat6 (e_ring, -1, e_orb)
    call twiss_at_start (e_ring)
    call twiss_propagate_all (e_ring)
    p_ring = ring
    p_ring%param%particle = positron$
    call twiss_at_start (p_ring)
    call closed_orbit_calc (p_ring, p_orb, 4)
    call lat_make_mat6 (p_ring, -1, p_orb)
    call twiss_at_start (p_ring)
    call twiss_propagate_all (p_ring)
  else
    p_ring = ring
    p_ring%param%particle = positron$
    call twiss_at_start (p_ring)
    call closed_orbit_calc (p_ring, p_orb, 4)
    call lat_make_mat6 (p_ring, -1, p_orb)
    call twiss_at_start (p_ring)
    call twiss_propagate_all (p_ring)
    e_ring = ring
    e_ring%param%particle = electron$
    call closed_orbit_calc (e_ring, e_orb, 4)
    call lat_make_mat6 (e_ring, -1, e_orb)
    call twiss_at_start (e_ring)
    call twiss_propagate_all (e_ring)
  endif

  call init_wall (walls%positive_x_wall)
  call init_wall (walls%negative_x_wall)

  call calculate_wind_power (p_ring, p_orb, +1, p_power, walls, gen_params, window)
  call calculate_wind_power (e_ring, e_orb, -1, e_power, walls, gen_params, window)

  call sigma_at_windows ( window, gen_params )

  deallocate( e_orb )
  deallocate( p_orb )
  deallocate( e_power )
  deallocate( p_power )

  print *,' Synchrotron radiation modeling and propagation is complete.'
end subroutine
