program long_term_tracking

use beam_mod
use twiss_and_track_mod
use ptc_map_with_radiation_mod

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele, ele_start
type (beam_init_struct) beam_init
type (bunch_struct), target :: bunch, bunch_init
type (coord_struct) orb_end, orb_init
type (coord_struct), allocatable :: closed_orb(:), orb(:)
type (coord_struct), pointer :: p
type (branch_struct), pointer :: branch
type (ele_pointer_struct), allocatable :: eles(:)
type (ptc_map_with_rad_struct) map_with_rad
type (normal_modes_struct) modes
type (rad_int_all_ele_struct) rad_int_ele

real(rp) time0, time, timer_print_dtime, chrom_x, chrom_y, ring_length, dead_cutoff
real(rp) average(6), sigma(6,6), print_on_dead_loss

integer ix_ele_start, n_turns, random_seed, map_order, n_turn_sigma_average
integer i, ip, n, ie, output_every_n_turns, ix_branch, ix_ele_end, iu_ave
integer iu_out, n_loc, track_state, i_turn, n_sum, n_live, n_live_old, n_print_dead_loss

logical err_flag, rfcavity_on, use_1_turn_map, add_closed_orbit_to_init_position
logical output_initial_position, err, is_lost, merge_particle_output_files

character(20) simulation_mode, map_mode, one_tracking_data_file
character(40) ele_start_name, fmt
character(200) lat_file, init_file, tracking_data_file, common_master_input_file, particle_output_file
character(200) sigma_matrix_output_file, input_map_file, output_map_file, map_file, averages_output_file

namelist / params / lat_file, n_turns, ele_start_name, map_order, beam_init, one_tracking_data_file, random_seed, &
                rfcavity_on, output_every_n_turns, bmad_com, add_closed_orbit_to_init_position, dead_cutoff, &
                simulation_mode, use_1_turn_map, tracking_data_file, timer_print_dtime, output_initial_position, &
                common_master_input_file, input_map_file, output_map_file, sigma_matrix_output_file, map_file, &
                print_on_dead_loss, averages_output_file, merge_particle_output_files, n_turn_sigma_average, &
                particle_output_file

! Parse command line

init_file = 'long_term_tracking.init'

if (cesr_iargc() == 1) then
  call cesr_getarg(1, init_file)
elseif (cesr_iargc() > 1) then
  print '(a)', 'Extra stuff on the command line!? Stopping here.'
  stop
endif

! Read parameters

dead_cutoff = 0
print_on_dead_loss = -1
random_seed = 0
output_initial_position = .false.
timer_print_dtime = 120
map_order = -1
ele_start_name = ''
rfcavity_on = .true.
output_every_n_turns = -1
n_turn_sigma_average = -1
random_seed = 0
add_closed_orbit_to_init_position = .true.
tracking_data_file = ''
one_tracking_data_file = ''
common_master_input_file = ''
input_map_file = ''
output_map_file = ''
sigma_matrix_output_file = ''
map_file = ''
particle_output_file = ''
merge_particle_output_files = .false.

print '(2a)', 'Initialization file: ', trim(init_file)

open (1, file = init_file, status = 'old')
read (1, nml = params)
close (1)

if (common_master_input_file /= '') then
  print '(2a)', 'Using common_master_input_file: ', trim(common_master_input_file)

  open (1, file = common_master_input_file, status = 'old')
  read (1, nml = params)
  close (1)

  open (1, file = init_file, status = 'old')
  read (1, nml = params)
  close (1)
endif

if (tracking_data_file /= '') then
  print '(a)', 'The "tracking_data_file" file has been renamed to "particle_output_file".'
  print '(a)', 'Please make the change in your input file.'
  stop
endif

if (one_tracking_data_file /= '') then
  print '(a)', 'The "one_tracking_data_file" has been renamed to "merge_particle_output_files".'
  print '(a)', 'Please make the change in your input file.'
  stop
endif

! Lattice init

bmad_com%auto_bookkeeper = .false.
call ran_seed_put (random_seed)

call bmad_parser (lat_file, lat)

! Read the master input file again so that bmad_com parameters set in the file
! take precedence over bmad_com parameters set in the lattice file.

open (1, file = init_file, status = 'old')
read (1, nml = params)  
close (1)

!

if (input_map_file /= '' .or. output_map_file /= '') then
  print '(a)', 'Note: "input_map_file and "output_map_file" parameters no longer exist.'
  print '(a)', 'Use the "map_file" parameter instead.'
  print '(a)', 'Stopping here.'
  stop
endif

if (map_file(1:6) == 'WRITE:') then
  map_file = map_file(7:)
  map_mode = 'write'
elseif (map_file == '') then
  map_mode = 'create'
else
  map_mode = 'read'
endif

call upcase_string(simulation_mode)

!

if (ele_start_name == '') Then
  ele_start => lat%ele(0)
else
  call lat_ele_locator (ele_start_name, lat, eles, n_loc, err)
  if (err .or. n_loc == 0) then
    print '(2a)', 'Starting element not found: ', trim(ele_start_name)
    stop
  endif
  if (n_loc > 1) then
    print '(2a)', 'Multiple elements found with name: ', trim(ele_start_name)
    print '(a)', 'Will stop here.'
    stop
  endif
  ele_start => eles(1)%ele
endif

ix_ele_start = ele_start%ix_ele
ix_branch = ele_start%ix_branch
branch => lat%branch(ix_branch)
ix_ele_end = ix_ele_start - 1
if (ix_ele_end == 0) ix_ele_end = branch%n_ele_track

call twiss_and_track (lat, closed_orb, ix_branch = ix_branch)

! Init 1-turn map.

if (map_mode == 'write' .or. ((use_1_turn_map .or. simulation_mode == 'CHECK') .and. simulation_mode /= 'STAT')) then
  use_1_turn_map = .true.
  if (map_mode == 'read') then
    call ptc_read_map_with_radiation(map_file, map_with_rad)
    print '(2a)',    'Map read in from file: ', trim(map_file)
    print '(2a)',    'Lattice file used for map:         ', trim(map_with_rad%lattice_file)
    print '(a, l1)', 'Map saved with radiation damping:  ', map_with_rad%radiation_damping_on
  else
    if (.not. rfcavity_on) print '(a)', 'RF will not be turned OFF since 1-turn map is in use!'
    print '(a)', 'Creating map file...'
    call run_timer ('START')
    call lat_to_ptc_layout (lat)
    call ptc_setup_map_with_radiation (map_with_rad, ele_start, ele_start, map_order, bmad_com%radiation_damping_on)
    call run_timer ('READ', time)
    print '(a, f8.2)', 'Map setup time (min)', time/60
    if (map_mode == 'write') then
      call ptc_write_map_with_radiation(map_file, map_with_rad)
      print '(2a)', 'Created map file: ', trim(map_file)
    endif
endif

else
  use_1_turn_map = .false.
  if (.not. rfcavity_on) call set_on_off (rfcavity$, lat, off$)
endif

! Print some info.

print '(a)', '--------------------------------------'
print '(a, a)',  'lattice:         ', trim(lat_file)
print '(a, a)',  'simulation_mode: ', simulation_mode
print '(a, l1)', 'Radiation Damping:           ', bmad_com%radiation_damping_on
print '(a, l1)', 'Stochastic Fluctuations:     ', bmad_com%radiation_fluctuations_on
print '(a, l1)', 'Spin_tracking_on:            ', bmad_com%spin_tracking_on
if (use_1_turn_map) then
  print '(a, i0)', 'Map order:          ', map_with_rad%map_order
endif
print '(a)', '--------------------------------------'

call run_timer ('START')

average = 0
sigma = 0
n_sum = 0

!-----------------------------------------
! A single turn tracking check

select case (simulation_mode)
case ('CHECK')
  call init_coord (orb_init, beam_init%center, ele_start, downstream_end$, spin = beam_init%spin)
  call reallocate_coord(orb, lat)

  if (add_closed_orbit_to_init_position) orb_init%vec = orb_init%vec + closed_orb(ix_ele_start)%vec

  orb_end = orb_init
  orb(ix_ele_start) = orb_init

  call ptc_track_map_with_radiation (orb_end, map_with_rad)
  call track_many (lat, orb, ix_ele_start, ix_ele_start, 1, ix_branch, track_state)

  print '(a)', 'Phase Space at Track End:'
  print '(a, 6f14.8)', '1-turn map tracking:', orb_end%vec
  print '(a, 6f14.8)', 'Ele-by-ele tracking:', orb(ix_ele_start)%vec
  print '(a, 6f14.8)', 'Diff:               ', orb_end%vec - orb(ix_ele_start)%vec

  print *
  print '(a)', 'Spin at Track End:'
  print '(a, 6f14.8)', '1-turn map tracking:', orb_end%spin
  print '(a, 6f14.8)', 'Ele-by-ele tracking:', orb(ix_ele_start)%spin
  print '(a, 6f14.8)', 'Diff:               ', orb_end%spin - orb(ix_ele_start)%spin

!-----------------------------------------
! Single particle tracking

case ('SINGLE')

  call reallocate_coord (orb, lat)
  call init_coord (orb(ix_ele_start), beam_init%center, ele_start, downstream_end$, lat%param%particle, spin = beam_init%spin)

  if (add_closed_orbit_to_init_position) orb(ix_ele_start)%vec = orb(ix_ele_start)%vec + closed_orb(ix_ele_start)%vec

  fmt = '(i6, 6es16.8, 3x, 3f10.6)'
  iu_out = lunget()
  if (particle_output_file == '') particle_output_file = 'single.dat'
  open(iu_out, file = particle_output_file, recl = 200)
  call write_this_header(iu_out, beam_init, map_with_rad, n_turns, simulation_mode)
  write (iu_out, '(a)') '# Turn |            x              px               y              py               z              pz    |   spin_x    spin_y    spin_z'
  write (iu_out, fmt) 0, orb(ix_ele_start)%vec, orb(ix_ele_start)%spin

  do i_turn = 1, n_turns
    if (use_1_turn_map) then
      orb_end = orb(ix_ele_start)
      call ptc_track_map_with_radiation (orb_end, map_with_rad)
      is_lost = orbit_too_large(orb_end)
    else
      call track_many (lat, orb, ix_ele_start, ix_ele_start, 1, ix_branch, track_state)
      is_lost = (track_state /= moving_forward$)
      orb_end = orb(ix_ele_start)
    endif

    if (output_every_n_turns > 0 .and. modulo(i_turn, output_every_n_turns) == 0) then
      write (iu_out, fmt) i_turn, orb_end%vec, orb_end%spin
    endif

    if (is_lost) then
      ele => branch%ele(track_state)
      print '(a, i0, 8a)', 'Particle lost at turn: ', i_turn
      if (.not. use_1_turn_map) print '(5a)', 'Lost at element: ', trim(ele%name), ' (', ele_location(ele), '), State: ', coord_state_name(orb(track_state)%state)
      exit
    endif

    if (sigma_matrix_output_file /= '') then
      average = average + orb_end%vec
      sigma = sigma + outer_product(orb_end%vec, orb_end%vec)
      n_sum = n_sum + 1
    endif
  enddo

  print '(2a)', 'Particle output file: ', trim(particle_output_file)
  close (iu_out)

  if (sigma_matrix_output_file /= '') then
    call write_sigma_file (sigma_matrix_output_file, n_sum, average, sigma)
  endif

!-----------------------------------------
! Beam tracking

case ('BUNCH')

  ! Init

  call init_bunch_distribution (ele_start, lat%param, beam_init, ix_branch, bunch, err_flag)
  if (err_flag) stop

  print '(a, i8)',   'n_particle             = ', size(bunch%particle)

  iu_ave = 0
  if (averages_output_file /= '') then
    print '(a)', 'Averages_output_file: ', trim(averages_output_file)
    iu_ave = lunget()
    open (iu_ave, file = averages_output_file, recl = 200)
    write (iu_ave, '(a1, a8, a9, a14, 3a14, 6a14)') '#', 'Turn', 'N_live', 'Polarization', &
                       '<Sx>', '<Sy>', '<Sz>', '<x>', '<px>', '<y>', '<py>', '<z>', '<pz>'
  endif

  do n = 1, size(bunch%particle)
    p => bunch%particle(n)
    if (add_closed_orbit_to_init_position) p%vec = p%vec + closed_orb(ix_ele_start)%vec
  enddo

  n_print_dead_loss = max(1, nint(print_on_dead_loss * size(bunch%particle)))
  allocate(bunch_init%particle(size(bunch%particle)))
  bunch_init%particle = bunch%particle
  n_live_old = size(bunch%particle)

  ! 

  call write_particle_data (0, bunch, map_with_rad)
  call write_averaged_data (0, bunch)

  time0 = 0

  do i_turn = 1, n_turns
    if (use_1_turn_map) then
      do ip = 1, size(bunch%particle)
        call ptc_track_map_with_radiation (bunch%particle(ip), map_with_rad)
      enddo
    else
      call track_bunch (lat, bunch, ele_start, ele_start, err_flag)
    endif

    n_live = count(bunch%particle%state == alive$)
    if ((n_live - n_live_old) >= n_print_dead_loss) then
      print '(a, i0, a, i0)', 'Number dead on turn ', i_turn, ': ', size(bunch%particle) - n_live
      n_live_old = n_live
    endif

    if (size(bunch%particle) - n_live > dead_cutoff * size(bunch%particle)) then
      print '(a)', 'Particle loss greater than set by dead_cutoff. Stopping now.'
      exit
    endif

    if (sigma_matrix_output_file /= '') then
      if (n_turn_sigma_average < 0 .or. n_turns - i_turn < n_turn_sigma_average) then
        do ip = 1, size(bunch%particle)
          p => bunch%particle(ip)
          if (p%state /= alive$) cycle
          average = average + p%vec
          sigma = sigma + outer_product(p%vec, p%vec)
          n_sum = n_sum + 1
        enddo
      endif
    endif

    call run_timer('READ', time)

    if (time-time0 > timer_print_dtime) then
      print '(a, f10.2, a, i0)', 'Ellapsed time (min): ', time/60, ', At turn: ', i_turn
      time0 = time
    endif

    call write_particle_data (i_turn, bunch, map_with_rad)
    call write_averaged_data (i_turn, bunch)
  end do

  print '(2a)', 'Tracking data file: ', trim(particle_output_file)

  if (sigma_matrix_output_file /= '') then
    call write_sigma_file (sigma_matrix_output_file, n_sum, average, sigma)
  endif

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
  call radiation_integrals (lat, closed_orb, modes, rad_int_by_ele = rad_int_ele, ix_branch = ix_branch)

  print *, 'Momentum Compaction:', modes%synch_int(1)/ring_length
  print *, 'dE/E=', modes%sigE_E
  print *, 'sig_z(m)=', modes%sig_z
  print *, 'emit_I  (m): ',  modes%a%emittance
  print *, 'emit_II (m): ',  modes%b%emittance
  print *, 'emit_III(m): ',  modes%z%emittance
  print *, 'QI =',ele%a%phi/twopi
  print *, 'QII=',ele%b%phi/twopi
  print *, 'QIII: ', lat%z%tune / twopi
  print *, '# of elements: ', lat%n_ele_track
  print *, 'L=',ring_length
  print *, 'dQI =',chrom_x
  print *, 'dQII=',chrom_y

!-----------------------------------------
! Unknown

case default
  print *, 'BAD SIMULATION_MODE: ' // simulation_mode
end select

call run_timer ('READ', time)
print '(a, f8.2)', 'Tracking time (min)', time/60

!------------------------------------------------------------------------------------------
contains

subroutine write_particle_data (i_turn, bunch, map_with_rad)

type (bunch_struct), target :: bunch
type (ptc_map_with_rad_struct) map_with_rad
type (coord_struct), pointer :: p, p0
integer i_turn, ix, ip, j
integer, save :: iu_snap = 0
character(200) file_name
character(40) fmt

!

if (particle_output_file == '') return
if (output_every_n_turns == -1 .and. i_turn /= n_turns) return
if (output_every_n_turns == 0 .and. i_turn /= 0 .and. i_turn /= n_turns) return
if (output_every_n_turns > 0 .and. modulo(i_turn, output_every_n_turns) /= 0) return

if (iu_snap == 0) then
  if (merge_particle_output_files) then
    file_name = particle_output_file
  else
    j = int(log10(real(n_turns, rp)) + 1 + 1d-10)
    write (fmt, '(a, i0, a, i0, a)') '(a, i', j, '.', j, ', a)'
    ix = index(particle_output_file, '#')
    if (ix == 0) then
      write (file_name, fmt) trim(particle_output_file), i_turn
    else
      write (file_name, fmt) particle_output_file(1:ix-1), i_turn, trim(particle_output_file(ix+1:))
    endif
  endif

  iu_snap = lunget()
  open (iu_snap, file = file_name, recl = 300)
  call write_this_header(iu_snap, beam_init, map_with_rad, n_turns, simulation_mode)

  if (output_initial_position) then
    write (iu_snap, '(2a)') '#      Ix     Turn | Start:    x              px               y              py               z              pz         spin_x    spin_y    spin_z  ', &
                                            '| End:         x              px               y              py               z              pz         spin_x    spin_y    spin_z    State'
  else
    write (iu_snap, '(a)')  '#      Ix     Turn |           x              px               y              py               z              pz   |     spin_x    spin_y    spin_z    State'
  endif
endif

do ip = 1, size(bunch%particle)
  p0 => bunch_init%particle(ip)
  p => bunch%particle(ip)
  if (output_initial_position) then
    write (iu_snap, '(i9, i9, 2(6es16.8, 3x, 3f10.6, 4x), a)') ip, i_turn, p0%vec, p0%spin, p%vec, p%spin, trim(coord_state_name(p%state))
  else
    write (iu_snap, '(i9, i9, 6es16.8, 3x, 3f10.6, 4x, a)')  ip, i_turn, p%vec, p%spin, trim(coord_state_name(p%state))
  endif
enddo

if (.not. merge_particle_output_files) then
  close(iu_snap)
  iu_snap = 0
endif

end subroutine write_particle_data

!------------------------------------------------------------------------------------------
! contains

subroutine write_this_header(iu, beam_init, map_with_rad, n_turns, simulation_mode)

type (beam_init_struct) beam_init
type (ptc_map_with_rad_struct) map_with_rad

integer iu, n_turns

character(20) simulation_mode

!

write (iu,  '(3a)')      '# lattice = "', trim(lat_file), '"'
write (iu,  '(3a)')      '# simulation_mode = "', trim(simulation_mode), '"'
write (iu,  '(a, i8)')   '# beam_init%n_particle   = ', beam_init%n_particle
write (iu,  '(a, i8)')   '# n_particle             = ', size(bunch%particle)
write (iu,  '(a, i8)')   '# n_turns                = ', n_turns
write (iu,  '(a, l1)')   '# Radiation_Damping      = ', bmad_com%radiation_damping_on
write (iu,  '(a, l1)')   '# Radiation_Fluctuations = ', bmad_com%radiation_fluctuations_on
write (iu,  '(a, l1)')   '# Spin_tracking_on       = ', bmad_com%spin_tracking_on
if (use_1_turn_map) then
  write (iu, '(3a)')     '# Map_file               = "', trim(map_file), '"'
  write (iu, '(a, i0)')  '# Map_order              = ', map_with_rad%map_order
else
  write (iu, '(a)') '#'
  write (iu, '(a)') '#'
endif
write (iu, '(a)') '#'

end subroutine write_this_header

!------------------------------------------------------------------------------------------
! contains

subroutine write_averaged_data (i_turn, bunch)

type (bunch_struct) bunch
real(rp) average(6), s_ave(3)
integer i, i_turn, n_live

!

if (iu_ave == 0) return
if (output_every_n_turns == -1 .and. i_turn /= n_turns) return
if (output_every_n_turns == 0 .and. i_turn /= 0 .and. i_turn /= n_turns) return
if (output_every_n_turns > 0 .and. modulo(i_turn, output_every_n_turns) /= 0) return

n_live = count(bunch%particle%state == alive$)

do i = 1, 6
  average(i) = sum(bunch%particle%vec(i), bunch%particle%state == alive$) / n_live
enddo

do i = 1, 3
  s_ave(i) = sum(bunch%particle%spin(i), bunch%particle%state == alive$) / n_live
enddo

write (iu_ave, '(i9, i9, f14.9, 2x, 3f14.9, 2x, 6es14.6)') i_turn, n_live, norm2(s_ave), s_ave, average

end subroutine write_averaged_data

!------------------------------------------------------------------------------------------
! contains

subroutine write_sigma_file (output_file, n_sum, average, sigma)

real(rp) average(6), sigma(6,6)
integer i, n_sum
character(*) output_file

!

if (n_sum == 0) then
  write (1, '(a, i0)') '# n_turn_sigma_average = ', n_turn_sigma_average
  write (1, '(a)') '# NO DATA TO AVERAGE OVER!'
  return
endif

!

average = average / n_sum
sigma = sigma / n_sum - outer_product(average, average)


open (1, file = sigma_matrix_output_file)

write (1, '(a, i0)') '# n_turn_sigma_average = ', n_turn_sigma_average
write (1, '(a)') '# Average:'
write (1, '(5x, 6es16.8)') average
write (1, *)
write (1, '(a)') '# Sigma:'
do i = 1, 6
  write (1, '(5x, 6es16.8)') sigma(i,:)
enddo

close (1)

print '(2a)', 'Sigma matrix data file: ', trim(sigma_matrix_output_file)

end subroutine write_sigma_file

end program
