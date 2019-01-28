program long_term_tracking

use beam_mod
use twiss_and_track_mod
use ptc_map_with_radiation_mod

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele, ele_start
type (beam_init_struct) beam_init
type (bunch_struct), target :: bunch, bunch_init
type (coord_struct) orb_1_turn, orb_init
type (coord_struct), allocatable :: closed_orb(:), orb(:)
type (coord_struct), pointer :: p
type (branch_struct), pointer :: branch
type (ele_pointer_struct), allocatable :: eles(:)
type (tree_element_zhe) map_with_rad(3)
type (normal_modes_struct) modes
type (rad_int_all_ele_struct) rad_int_ele

real(rp) time0, time, timer_print_dtime, chrom_x, chrom_y, ring_length

integer ix_ele_start, n_turns, random_seed, map_order
integer i, ip, n, ie, output_every_n_turns, ix_branch, ix_ele_end
integer iu_out, iu_snap, n_loc, track_state, i_turn

logical err_flag, rfcavity_on, use_1_turn_map, one_tracking_data_file, add_closed_orbit_to_init_position
logical include_initial_position, err, is_lost

character(20) simulation_mode
character(40) ele_start_name, fmt
character(200) lat_file, init_file, tracking_data_file

namelist / params / lat_file, n_turns, ele_start_name, map_order, beam_init, one_tracking_data_file, random_seed, &
                  rfcavity_on, output_every_n_turns, bmad_com, add_closed_orbit_to_init_position, &
                  simulation_mode, use_1_turn_map, tracking_data_file, timer_print_dtime, include_initial_position

! Parse command line

init_file = 'long_term_tracking.init'

if (cesr_iargc() == 1) then
  call cesr_getarg(1, init_file)
elseif (cesr_iargc() > 1) then
  print *, 'Extra stuff on the command line!? Stopping here.'
  stop
endif

! Read parameters

random_seed = 0
include_initial_position = .false.
timer_print_dtime = 120
map_order = -1
ele_start_name = ''
rfcavity_on = .true.
output_every_n_turns = -1
random_seed = 0
add_closed_orbit_to_init_position = .true.
tracking_data_file = ''
one_tracking_data_file = .false.

print '(2a)', 'Initialization file: ', trim(init_file)

open (1, file = init_file, status = 'old')
read (1, nml = params)
close (1)

! Lattice init

call ran_seed_put (random_seed)
call bmad_parser (lat_file, lat)

! Read the master input file again so that bmad_com parameters set in the file
! take precedence over bmad_com parameters set in the lattice file.

open (1, file = init_file, status = 'old')
read (1, nml = params)  
close (1)

call upcase_string(simulation_mode)

!

if (ele_start_name == '') Then
  ele_start => lat%ele(0)
else
  call lat_ele_locator (ele_start_name, lat, eles, n_loc, err)
  if (err .or. n_loc == 0) then
    print *, 'Starting element not found: ', trim(ele_start_name)
    stop
  endif
  if (n_loc > 1) then
    print *, 'Multiple elements found with name: ', trim(ele_start_name)
    print *, 'Will stop here.'
    stop
  endif
  ele_start => eles(1)%ele
endif

ix_ele_start = ele_start%ix_ele
ix_branch = ele_start%ix_branch
branch => lat%branch(ix_branch)
ix_ele_end = ix_ele_start - 1
if (ix_ele_end == 0) ix_ele_end = branch%n_ele_track

! Init 1-turn map.

if ((use_1_turn_map .or. simulation_mode == 'CHECK') .and. simulation_mode /= 'STAT') then
  if (.not. rfcavity_on) print *, 'RF will not be turned OFF since 1-turn map is in use!'
  call run_timer ('START')
  call lat_to_ptc_layout (lat)
  call ptc_setup_map_with_radiation (map_with_rad, ele_start, ele_start, map_order = map_order)
  call run_timer ('READ', time)
  print '(a, f8.2)', '1-turn map setup time (min)', time/60

else
  if (.not. rfcavity_on) call set_on_off (rfcavity$, lat, off$)
endif

call twiss_and_track (lat, closed_orb, ix_branch = ix_branch)

! Print some info.

print *, '--------------------------------------'
print '(a, a)',  'lattice:         ', trim(lat_file)
print '(a, a)',  'simulation_mode: ', simulation_mode
print '(a, i8)', 'n_particle = ', beam_init%n_particle
print '(a, i8)', 'n_bunch    = ', beam_init%n_bunch
print '(a, i8)', 'n_turns    = ', n_turns
print *, 'bmad_com%max_aperture_limit: ', bmad_com%max_aperture_limit
print *, 'bmad_com%aperture_limit_on:  ', bmad_com%aperture_limit_on
print *, 'Radiation Damping:           ', bmad_com%radiation_damping_on
print *, 'Stochastic Fluctuations:     ', bmad_com%radiation_fluctuations_on
print *, 'Spin_tracking_on:            ', bmad_com%spin_tracking_on
print *, '--------------------------------------'

call run_timer ('START')

!-----------------------------------------
! A single turn tracking check

select case (simulation_mode)
case ('CHECK')
  call init_coord (orb_init, beam_init%center, ele_start, downstream_end$, spin = beam_init%spin)
  call reallocate_coord(orb, lat)

  if (add_closed_orbit_to_init_position) orb_init%vec = orb_init%vec + closed_orb(ix_ele_start)%vec

  orb_1_turn = orb_init
  orb(ix_ele_start) = orb_init

  call ptc_track_map_with_radiation (orb_1_turn, map_with_rad, .true., .true.)
  call track_many (lat, orb, ix_ele_start, ix_ele_start, 1, ix_branch, track_state)

  print *, 'Phase Space:'
  print '(a, 6f14.8)', '1-turn map:', orb_1_turn%vec
  print '(a, 6f14.8)', 'Ele-by-ele:', orb(ix_ele_start)%vec
  print '(a, 6f14.8)', 'Diff:      ', orb_1_turn%vec - orb(ix_ele_start)%vec

  print *
  print *, 'Spin:'
  print '(a, 6f14.8)', '1-turn map:', orb_1_turn%spin
  print '(a, 6f14.8)', 'Ele-by-ele:', orb(ix_ele_start)%spin
  print '(a, 6f14.8)', 'Diff:      ', orb_1_turn%spin - orb(ix_ele_start)%spin

!-----------------------------------------
! Single particle tracking

case ('SINGLE')

  call reallocate_coord (orb, lat)
  call init_coord (orb(ix_ele_start), beam_init%center, ele_start, downstream_end$, lat%param%particle, spin = beam_init%spin)

  if (add_closed_orbit_to_init_position) orb(ix_ele_start)%vec = orb(ix_ele_start)%vec + closed_orb(ix_ele_start)%vec

  fmt = '(i6, a, 3x, 3f10.6)'
  iu_out = lunget()
  if (tracking_data_file == '') tracking_data_file = 'single.dat'
  open(iu_out, file = tracking_data_file, recl = 200)
  write (iu_out, '(a)') ' #                                             * 10^3                                  |'           
  write (iu_out, '(a)') ' #Turn           x           px            y           py            z           pz    |   spin_x    spin_y    spin_z'
  write (iu_out, fmt) 0, reals_to_table_row(1d3*orb(ix_ele_start)%vec, 13, 7, 1), orb(ix_ele_start)%spin

  do i_turn = 1, n_turns
    if (use_1_turn_map) then
      orb_1_turn = orb(ix_ele_start)
      call ptc_track_map_with_radiation (orb_1_turn, map_with_rad, .true., .true.)
      is_lost = orbit_too_large(orb_1_turn)
      orb(ix_ele_start) = orb_1_turn
    else
      call track_many (lat, orb, ix_ele_start, ix_ele_start, 1, ix_branch, track_state)
      is_lost = (track_state /= moving_forward$)
    endif

    if (output_every_n_turns > 0 .and. modulo(i_turn, output_every_n_turns) == 0) then
      write (iu_out, fmt) i_turn, reals_to_table_row(1d3*orb(ix_ele_start)%vec, 13, 7, 1), orb(ix_ele_start)%spin
    endif

    if (is_lost) then
      ele => branch%ele(track_state)
      print '(a, i0, 8a)', 'Particle lost at turn: ', i_turn
      if (.not. use_1_turn_map) print '(5a)', 'Lost at element: ', trim(ele%name), ' (', ele_loc_to_string(ele), '), State: ', coord_state_name(orb(track_state)%state)
      exit
    endif
  enddo

  print *, 'Data file: ', trim(tracking_data_file)

!-----------------------------------------
! Beam tracking

case ('BUNCH')

  ! Init

  if (tracking_data_file == '') tracking_data_file = 'tracking.dat'
  iu_snap = 0

  call init_bunch_distribution (ele_start, lat%param, beam_init, ix_branch, bunch, err_flag)
  if (err_flag) stop

  do n = 1, size(bunch%particle)
    p => bunch%particle(n)
    if (add_closed_orbit_to_init_position) p%vec = p%vec + closed_orb(ix_ele_start)%vec
  enddo

  allocate(bunch_init%particle(size(bunch%particle)))
  bunch_init%particle = bunch%particle

  ! 

  call write_tracking_data (0)
  time0 = time

  do i_turn = 1, n_turns
    if (use_1_turn_map) then
      do ip = 1, size(bunch%particle)
        call ptc_track_map_with_radiation (bunch%particle(ip), map_with_rad, .true., .true.)
      enddo
    else
      call track_bunch (lat, bunch, ele_start, ele_start, err_flag)
    endif

    call run_timer('READ', time)

    if (time-time0 > timer_print_dtime) then
      print '(a, f10.2, a, i0)', 'Ellapsed time (min): ', time/60, ', At turn: ', i_turn
      time0 = time
    endif

    call write_tracking_data (i_turn)
  end do

  close(iu_snap)

!-----------------------------------------
! Lattice statistics (radiation integrals, etc.).

case ('STAT')

  open(unit=20,file="twiss.dat")
  open(unit=21,file="coupling.dat")
  open (unit=22, file = "closed_orbit.dat")

  write(20,'(52x,"s",11x,"b_x",9x,"b_y",9x,"a_x",9x,"a_y",8x,"mu_x",8x,"mu_y",9x,"D_x",9x,"D_y")')
  write(22, '(a)') "Closed orbit"
  write(22, '(a, 7x, a, 9x, a, 11x, a, 8x, a)') "pos", "x", "y", "dz", "dE/E"

  do i = 0, lat%n_ele_track
    ele => branch%ele(i)
    write(20,'(i4,2x,a16,2x,a,3f12.4,2x,6f12.6)') i, ele%name, key_name(ele%key), ele%s, &
                      ele%a%beta, ele%b%beta, ele%a%alpha, ele%b%alpha, ele%a%phi, ele%b%phi, ele%a%eta,  ele%b%eta
    write(21,'(i4,2x,a16,2x,a,1x,f10.3,4(1x,f12.6))') i, ele%name, key_name(ele%key), &
                                        ele%s,ele%c_mat(1,1),ele%c_mat(1,2),ele%c_mat(2,1),ele%c_mat(2,2)
    write(22,'(a,1x,f10.3,4(1x,f10.6))') ele%name, ele%s, closed_orb(i)%vec(1:5:2), closed_orb(i)%vec(6)
  enddo

  close(20)
  close(21)
  close(22)

  !

  ring_length = branch%param%total_length
  call chrom_calc(lat, 1.0d-6, chrom_x, chrom_y, err_flag, ix_branch = ix_branch)
  call calc_z_tune (lat)
  call radiation_integrals (lat, orb, modes, rad_int_by_ele = rad_int_ele, ix_branch = ix_branch)

  print *, 'Momentum Compaction:', modes%synch_int(1)/ring_length
  print *, 'dE/E=', modes%sigE_E
  print *, 'sig_z(m)=', modes%sig_z
  print *,'emit_I  (m)     : ',  modes%a%emittance
  print *,'emit_II (m)     : ',  modes%b%emittance
  print *,'emit_III(m)     : ',  modes%z%emittance
  print*, 'QI =',ele%a%phi/twopi
  print*, 'QII=',ele%b%phi/twopi
  print *,'QIII      : ', lat%z%tune / twopi
  print*, '# of elements: ', lat%n_ele_track
  print*, 'L=',ring_length
  print*, 'dQI =',chrom_x
  print*, 'dQII=',chrom_y

!-----------------------------------------
! Unknown

case default
  print *, 'BAD SIMULATION_MODE: ' // simulation_mode
end select

call run_timer ('READ', time)
print '(a, f8.2)', 'Tracking time (min)', time/60

!------------------------------------------------------------------------------------------
contains

subroutine write_tracking_data (nn)

type (coord_struct), pointer :: p, p0
integer nn, ix, ip, j
character(200) file_name
character(40) fmt

!

if (output_every_n_turns == -1 .and. nn /= n_turns) return
if (output_every_n_turns == 0 .and. nn /= 0 .and. nn /= n_turns) return
if (output_every_n_turns > 0 .and. modulo(nn, output_every_n_turns) /= 0) return

if (iu_snap == 0) then
  if (one_tracking_data_file) then
    file_name = tracking_data_file
  else
    j = int(log10(real(n_turns, rp)) + 1 + 1d-10)
    write (fmt, '(a, i0, a, i0, a)') '(a, i', j, '.', j, ', a)'
    ix = index(tracking_data_file, '#')
    if (ix == 0) then
      write (file_name, fmt) trim(tracking_data_file), nn
    else
      write (file_name, fmt) tracking_data_file(1:ix-1), nn, trim(tracking_data_file(ix+1:))
    endif
  endif

  iu_snap = lunget()
  open (iu_snap, file = file_name, recl = 300)

  if (include_initial_position) then
    write (iu_snap, '(2a)') '#                  |                                  Start  * 10^3                              |                Start            ', &
                                            '|                                       End  * 10^3                              |                 End'
    write (iu_snap, '(2a)') '#      Ix     Turn |        x           px            y           py            z           pz   |     spin_x    spin_y    spin_z  ', &
                                            '|           x           px            y           py            z           pz   |     spin_x    spin_y    spin_z    State'
  else
    write (iu_snap, '(a)')  '#                  |                                     * 10^3                                  |'
    write (iu_snap, '(a)')  '#      Ix     Turn |        x           px            y           py            z           pz   |     spin_x    spin_y    spin_z    State'
  endif
endif

do ip = 1, size(bunch%particle)
  p0 => bunch_init%particle(ip)
  p => bunch%particle(ip)
  if (include_initial_position) then
    write (iu_snap, '(i9, i9, 2(a, 3x, 3f10.6, 4x), a)') ip, nn, reals_to_table_row(1d3*p0%vec, 13, 7, 1), &
                              p0%spin, reals_to_table_row(1d3*p%vec, 13, 7, 1), p%spin, trim(coord_state_name(p%state))
  else
    write (iu_snap, '(i9, i9, a, 3x, 3f10.6, 4x, a)')  ip, nn, reals_to_table_row(1d3*p%vec, 13, 7, 1), p%spin, coord_state_name(p%state)
  endif
enddo

if (.not. one_tracking_data_file) then
  close(iu_snap)
  iu_snap = 0
endif

end subroutine write_tracking_data

end program
