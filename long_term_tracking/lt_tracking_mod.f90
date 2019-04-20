!+
! Module lt_tracking_mod
!
! Routines used by the long_term_tracking program.
!-

module lt_tracking_mod

use beam_mod
use twiss_and_track_mod
use ptc_map_with_radiation_mod

implicit none

integer, parameter :: master_rank$ = 0
integer, parameter :: results_tag$ = 1000
integer, parameter :: is_done_tag$ = 1001

type ltt_params_struct
  type (ele_struct), pointer :: ele_start
  character(20) :: simulation_mode = ''
  character(40) :: ele_start_name = ''
  character(200) :: lat_file = ''
  character(200) :: common_master_input_file = ''
  character(200) :: particle_output_file = ''
  character(200) :: sigma_matrix_output_file = ''
  character(200) :: map_file = ''
  character(200) :: averages_output_file = ''
  integer :: n_turns = -1
  integer :: random_seed = -1
  integer :: map_order = -1
  integer :: n_turn_sigma_average = -1
  integer :: output_every_n_turns = -1
  real(rp) :: print_on_dead_loss = -1
  real(rp) :: timer_print_dtime = 120
  real(rp) :: dead_cutoff = 0
  logical :: rfcavity_on = .true.
  logical :: use_1_turn_map = .false.
  logical :: add_closed_orbit_to_init_position = .true.
  logical :: output_initial_position = .false.
  logical :: merge_particle_output_files = .false.
  integer :: mpi_rank  = master_rank$
  integer :: mpi_n_proc = 1                    ! Number of processeses including master
  integer :: mpi_num_runs = 10                 ! Number of runs a slave process will take on average.
  integer :: mpi_n_particles_per_run = 0        ! Number of particles per run.
  logical :: using_mpi = .false.
  logical :: debug = .false.
end type

type ltt_sum_data_struct
  integer :: i_turn = 0
  integer :: n_live = 0
  real(rp) :: orb_sum(6) = 0    ! Orbit average
  real(rp) :: spin_sum(3) = 0   ! Spin
end type

contains

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------


subroutine ltt_init_params(lttp, lat, beam_init)

type (ltt_params_struct) lttp
type (lat_struct) lat
type (beam_init_struct) beam_init

real(rp) timer_print_dtime, dead_cutoff, print_on_dead_loss

integer n_turns, random_seed, map_order, n_turn_sigma_average, output_every_n_turns, ir

logical rfcavity_on, use_1_turn_map, add_closed_orbit_to_init_position
logical output_initial_position, merge_particle_output_files

character(20) simulation_mode, map_mode, one_tracking_data_file
character(40) ele_start_name
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

use_1_turn_map = .false.
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

if (.not. lttp%using_mpi .or. lttp%mpi_rank == master_rank$) then
  print '(2a)', 'Initialization file: ', trim(init_file)
endif

open (1, file = init_file, status = 'old', action = 'read')
read (1, nml = params)
close (1)

if (common_master_input_file /= '') then
  print '(2a)', 'Using common_master_input_file: ', trim(common_master_input_file)

  open (1, file = common_master_input_file, status = 'old', action = 'read')
  read (1, nml = params)
  close (1)

  open (1, file = init_file, status = 'old', action = 'read')
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
call ptc_set_map_with_radiation_ran_seed(random_seed)

if (lttp%using_mpi) then
  call ran_seed_get (ir)
  call ran_seed_put (ir + 10 * lttp%mpi_rank)
  call ptc_set_map_with_radiation_ran_seed(ir + 10 * lttp%mpi_rank)
endif

call bmad_parser (lat_file, lat)

! Read the master input file again so that bmad_com parameters set in the file
! take precedence over bmad_com parameters set in the lattice file.

open (1, file = init_file, status = 'old', action = 'read')
read (1, nml = params)  
close (1)

!

if (input_map_file /= '' .or. output_map_file /= '') then
  print '(a)', 'Note: "input_map_file and "output_map_file" parameters no longer exist.'
  print '(a)', 'Use the "map_file" parameter instead.'
  print '(a)', 'Stopping here.'
  stop
endif

!

lttp%simulation_mode                     = simulation_mode
lttp%ele_start_name                      = ele_start_name
lttp%lat_file                            = lat_file
lttp%common_master_input_file            = common_master_input_file
lttp%particle_output_file                = particle_output_file
lttp%sigma_matrix_output_file            = sigma_matrix_output_file
lttp%map_file                            = map_file
lttp%averages_output_file                = averages_output_file
lttp%n_turns                             = n_turns
lttp%random_seed                         = random_seed
lttp%map_order                           = map_order
lttp%n_turn_sigma_average                = n_turn_sigma_average
lttp%output_every_n_turns                = output_every_n_turns
lttp%timer_print_dtime                   = timer_print_dtime
lttp%dead_cutoff                         = dead_cutoff
lttp%print_on_dead_loss                  = print_on_dead_loss
lttp%rfcavity_on                         = rfcavity_on
lttp%use_1_turn_map                      = use_1_turn_map
lttp%add_closed_orbit_to_init_position   = add_closed_orbit_to_init_position
lttp%output_initial_position             = output_initial_position
lttp%merge_particle_output_files         = merge_particle_output_files

call upcase_string(lttp%simulation_mode)

end subroutine ltt_init_params

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_init_tracking(lttp, lat, closed_orb, map_with_rad)

type (ltt_params_struct) lttp
type (lat_struct), target :: lat
type (ptc_map_with_rad_struct) map_with_rad
type (ele_pointer_struct), allocatable :: eles(:)
type (coord_struct), allocatable :: closed_orb(:), orb(:)
type (branch_struct), pointer :: branch

real(rp) time

integer n_loc, ix_ele_start, ix_ele_end, ix_branch

logical err, ok

character(20) map_mode

!

if (lttp%map_file(1:6) == 'WRITE:') then
  lttp%map_file = lttp%map_file(7:)
  map_mode = 'write'
  if (lttp%mpi_rank /= master_rank$) map_mode = 'read'  ! Master has created the file

elseif (lttp%map_file == '') then
  if (lttp%using_mpi) then
    lttp%map_file = 'map.dat'
    map_mode = 'write'
    if (lttp%mpi_rank /= master_rank$) map_mode = 'read'  ! Master has created the file
  else
    map_mode = 'create'
  endif
  
else
  inquire (file = lttp%map_file, exist = ok)
  if (ok) then
    map_mode = 'read'
  else
    map_mode = 'write'
    print *, 'Map_file does not exist: ', trim(lttp%map_file)
    print *, 'A map_file will be created.'
  endif
endif

!

if (lttp%ele_start_name == '') Then
  lttp%ele_start => lat%ele(0)
else
  call lat_ele_locator (lttp%ele_start_name, lat, eles, n_loc, err)
  if (err .or. n_loc == 0) then
    print '(2a)', 'Starting element not found: ', trim(lttp%ele_start_name)
    stop
  endif
  if (n_loc > 1) then
    print '(2a)', 'Multiple elements found with name: ', trim(lttp%ele_start_name)
    print '(a)', 'Will stop here.'
    stop
  endif
  lttp%ele_start => eles(1)%ele
endif

ix_ele_start = lttp%ele_start%ix_ele
ix_branch = lttp%ele_start%ix_branch
branch => lat%branch(ix_branch)
ix_ele_end = ix_ele_start - 1
if (ix_ele_end == 0) ix_ele_end = branch%n_ele_track

call twiss_and_track (lat, closed_orb, ix_branch = ix_branch)

! Init 1-turn map.

if (map_mode == 'write' .or. ((lttp%use_1_turn_map .or. lttp%simulation_mode == 'CHECK') .and. lttp%simulation_mode /= 'STAT')) then
  lttp%use_1_turn_map = .true.
  if (map_mode == 'read') then
    call ptc_read_map_with_radiation(lttp%map_file, map_with_rad)
    if (lttp%mpi_rank == master_rank$) then
      print '(2a)',    'Map read in from file: ', trim(lttp%map_file)
      print '(2a)',    'Lattice file used for map:         ', trim(map_with_rad%lattice_file)
      print '(a, l1)', 'Map saved with radiation damping:  ', map_with_rad%radiation_damping_on
    endif
  else
    if (.not. lttp%rfcavity_on) print '(a)', 'RF will not be turned OFF since 1-turn map is in use!'
    print '(a)', 'Creating map file...'
    call run_timer ('START')
    call lat_to_ptc_layout (lat)
    call ptc_setup_map_with_radiation (map_with_rad, lttp%ele_start, lttp%ele_start, lttp%map_order, bmad_com%radiation_damping_on)
    call run_timer ('READ', time)
    print '(a, f8.2)', 'Map setup time (min)', time/60
    if (map_mode == 'write') then
      call ptc_write_map_with_radiation(lttp%map_file, map_with_rad)
      print '(2a)', 'Created map file: ', trim(lttp%map_file)
    endif
endif

else
  lttp%use_1_turn_map = .false.
  if (.not. lttp%rfcavity_on) call set_on_off (rfcavity$, lat, off$)
endif

end subroutine ltt_init_tracking

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_print_inital_info (lttp, map_with_rad)

type (ltt_params_struct) lttp
type (ptc_map_with_rad_struct) map_with_rad

! Print some info.

print '(a)', '--------------------------------------'
print '(a, a)',  'lattice:         ', trim(lttp%lat_file)
print '(a, a)',  'simulation_mode: ', lttp%simulation_mode
print '(a, l1)', 'Radiation Damping:           ', bmad_com%radiation_damping_on
print '(a, l1)', 'Stochastic Fluctuations:     ', bmad_com%radiation_fluctuations_on
print '(a, l1)', 'Spin_tracking_on:            ', bmad_com%spin_tracking_on
if (lttp%use_1_turn_map) then
  print '(a, i0)', 'Map order:          ', map_with_rad%map_order
endif
print '(a)', '--------------------------------------'

end subroutine ltt_print_inital_info

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_run_check_mode (lttp, lat, map_with_rad, beam_init, closed_orb)

type (ltt_params_struct) lttp
type (lat_struct), target :: lat
type (ptc_map_with_rad_struct) map_with_rad
type (beam_init_struct) beam_init
type (coord_struct), allocatable :: closed_orb(:), orb(:)
type (coord_struct) orb_end, orb_init

integer track_state, ix_ele_start

! Run serial in check mode.

ix_ele_start = lttp%ele_start%ix_ele
call init_coord (orb_init, beam_init%center, lttp%ele_start, downstream_end$, spin = beam_init%spin)
call reallocate_coord(orb, lat)

if (lttp%add_closed_orbit_to_init_position) orb_init%vec = orb_init%vec + closed_orb(ix_ele_start)%vec

orb_end = orb_init
orb(lttp%ele_start%ix_ele) = orb_init

call ptc_track_map_with_radiation (orb_end, map_with_rad)
call track_many (lat, orb, lttp%ele_start%ix_ele, lttp%ele_start%ix_ele, 1, lttp%ele_start%ix_branch, track_state)

print '(a)', 'Phase Space at Track End:'
print '(a, 6f14.8)', '1-turn map tracking:', orb_end%vec
print '(a, 6f14.8)', 'Ele-by-ele tracking:', orb(ix_ele_start)%vec
print '(a, 6f14.8)', 'Diff:               ', orb_end%vec - orb(ix_ele_start)%vec

print *
print '(a)', 'Spin at Track End:'
print '(a, 6f14.8)', '1-turn map tracking:', orb_end%spin
print '(a, 6f14.8)', 'Ele-by-ele tracking:', orb(ix_ele_start)%spin
print '(a, 6f14.8)', 'Diff:               ', orb_end%spin - orb(ix_ele_start)%spin

end subroutine ltt_run_check_mode

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_run_single_mode (lttp, lat, beam_init, closed_orb, map_with_rad)

type (ltt_params_struct) lttp
type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (beam_init_struct) beam_init
type (coord_struct), allocatable :: closed_orb(:), orb(:)
type (coord_struct) orb_end
type (ele_struct), pointer :: ele
type (ptc_map_with_rad_struct) map_with_rad

real(rp) average(6), sigma(6,6)
integer n_sum, ix_ele_start, iu_out, i_turn, track_state, ix_branch
logical is_lost

character(40) fmt

! Run serial in single mode.

n_sum = 0
sigma = 0
average = 0
ix_ele_start = lttp%ele_start%ix_ele
ix_branch = lttp%ele_start%ix_branch
branch => lat%branch(ix_branch)

call reallocate_coord (orb, lat)
call init_coord (orb(ix_ele_start), beam_init%center, lttp%ele_start, downstream_end$, lat%param%particle, spin = beam_init%spin)

if (lttp%add_closed_orbit_to_init_position) orb(ix_ele_start)%vec = orb(ix_ele_start)%vec + closed_orb(ix_ele_start)%vec

fmt = '(i6, 6es16.8, 3x, 3f10.6)'
iu_out = lunget()
if (lttp%particle_output_file == '') lttp%particle_output_file = 'single.dat'
open(iu_out, file = lttp%particle_output_file, recl = 200)
call ltt_write_this_header(lttp, iu_out, map_with_rad, 1)
write (iu_out, '(a)') '# Turn |            x              px               y              py               z              pz    |   spin_x    spin_y    spin_z'
write (iu_out, fmt) 0, orb(ix_ele_start)%vec, orb(ix_ele_start)%spin

do i_turn = 1, lttp%n_turns
  if (lttp%use_1_turn_map) then
    orb_end = orb(ix_ele_start)
    call ptc_track_map_with_radiation (orb_end, map_with_rad)
    is_lost = orbit_too_large(orb_end)
  else
    call track_many (lat, orb, ix_ele_start, ix_ele_start, 1, ix_branch, track_state)
    is_lost = (track_state /= moving_forward$)
    orb_end = orb(ix_ele_start)
  endif

  if (lttp%output_every_n_turns > 0 .and. modulo(i_turn, lttp%output_every_n_turns) == 0) then
    write (iu_out, fmt) i_turn, orb_end%vec, orb_end%spin
  endif

  if (is_lost) then
    ele => branch%ele(track_state)
    print '(a, i0, 8a)', 'Particle lost at turn: ', i_turn
    if (.not. lttp%use_1_turn_map) print '(5a)', 'Lost at element: ', trim(ele%name), ' (', ele_location(ele), '), State: ', coord_state_name(orb(track_state)%state)
    exit
  endif

  if (lttp%sigma_matrix_output_file /= '') then
    average = average + orb_end%vec
    sigma = sigma + outer_product(orb_end%vec, orb_end%vec)
    n_sum = n_sum + 1
  endif
enddo

print '(2a)', 'Particle output file: ', trim(lttp%particle_output_file)
close (iu_out)

if (lttp%sigma_matrix_output_file /= '') then
  call ltt_write_sigma_file (lttp, n_sum, average, sigma)
endif

end subroutine ltt_run_single_mode

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_run_bunch_mode (lttp, lat, beam_init, closed_orb, map_with_rad, sum_data_array)

type (ltt_params_struct) lttp
type (lat_struct) lat
type (beam_init_struct) beam_init
type (coord_struct), allocatable :: closed_orb(:), orb(:)
type (ptc_map_with_rad_struct) map_with_rad
type (bunch_struct), target :: bunch, bunch_init
type (coord_struct), pointer :: p
type (ltt_sum_data_struct), allocatable, optional :: sum_data_array(:)
type (ltt_sum_data_struct), allocatable, target :: sum_data_arr(:)
type (ltt_sum_data_struct), pointer :: sd

real(rp) time0, time
real(rp) average(6), sigma(6,6)

integer n, n_sum, ix_ele_start, n_print_dead_loss, i_turn
integer ip, n_live_old, n_live

logical err_flag

! Init

n_sum = 0
sigma = 0
average = 0
ix_ele_start = lttp%ele_start%ix_ele

if (lttp%using_mpi) beam_init%n_particle = lttp%mpi_n_particles_per_run
call init_bunch_distribution (lttp%ele_start, lat%param, beam_init, lttp%ele_start%ix_branch, bunch, err_flag)
if (err_flag) stop

if (lttp%mpi_rank == master_rank$) print '(a, i8)',   'n_particle             = ', size(bunch%particle)

do n = 1, size(bunch%particle)
  p => bunch%particle(n)
  if (lttp%add_closed_orbit_to_init_position) p%vec = p%vec + closed_orb(ix_ele_start)%vec
enddo

n_print_dead_loss = max(1, nint(lttp%print_on_dead_loss * size(bunch%particle)))
allocate(bunch_init%particle(size(bunch%particle)))
bunch_init%particle = bunch%particle
n_live_old = size(bunch%particle)

! 

if (.not. lttp%using_mpi) call ltt_write_particle_data (lttp, 0, bunch, map_with_rad)
call ltt_calc_bunch_sums (lttp, 0, bunch, sum_data_arr)

time0 = 0

do i_turn = 1, lttp%n_turns
  if (lttp%use_1_turn_map) then
    do ip = 1, size(bunch%particle)
      call ptc_track_map_with_radiation (bunch%particle(ip), map_with_rad)
    enddo
  else
    call track_bunch (lat, bunch, lttp%ele_start, lttp%ele_start, err_flag)
  endif

  n_live = count(bunch%particle%state == alive$)
  if ((n_live - n_live_old) >= n_print_dead_loss) then
    print '(a, i0, a, i0)', 'Number dead on turn ', i_turn, ': ', size(bunch%particle) - n_live
    n_live_old = n_live
  endif

  if (size(bunch%particle) - n_live > lttp%dead_cutoff * size(bunch%particle) .and. .not. lttp%using_mpi) then
    print '(a)', 'Particle loss greater than set by dead_cutoff. Stopping now.'
    exit
  endif

  if (lttp%sigma_matrix_output_file /= '') then
    if (lttp%n_turn_sigma_average < 0 .or. lttp%n_turns - i_turn < lttp%n_turn_sigma_average) then
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

  if (time-time0 > lttp%timer_print_dtime) then
    print '(a, f10.2, a, i0)', 'Ellapsed time (min): ', time/60, ', At turn: ', i_turn
    time0 = time
  endif

  if (.not. lttp%using_mpi) call ltt_write_particle_data (lttp, i_turn, bunch, map_with_rad)
  call ltt_calc_bunch_sums (lttp, i_turn, bunch, sum_data_arr)
end do

if (lttp%mpi_rank == master_rank$) print '(2a)', 'Tracking data file: ', trim(lttp%particle_output_file)

if (lttp%sigma_matrix_output_file /= '' .and. .not. lttp%using_mpi) then
  call ltt_write_sigma_file (lttp, n_sum, average, sigma)
endif

if (lttp%using_mpi) then
  if (allocated(sum_data_array)) deallocate (sum_data_array)
  call move_alloc (sum_data_arr, sum_data_array)
else
  call ltt_write_bunch_averages(lttp, sum_data_arr)
endif

end subroutine ltt_run_bunch_mode

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_run_stat_mode (lttp, lat, closed_orb)

type (ltt_params_struct) lttp
type (lat_struct), target :: lat
type (coord_struct), allocatable :: closed_orb(:)
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele
type (normal_modes_struct) modes
type (rad_int_all_ele_struct) rad_int_ele

real(rp) chrom_x, chrom_y, ring_length
integer i
logical err_flag

! Run serial in stat mode.

branch => lat%branch(lttp%ele_start%ix_branch)

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
call chrom_calc(lat, 1.0d-6, chrom_x, chrom_y, err_flag, ix_branch = branch%ix_branch)
call calc_z_tune (lat)
call radiation_integrals (lat, closed_orb, modes, rad_int_by_ele = rad_int_ele, ix_branch = branch%ix_branch)

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

end subroutine ltt_run_stat_mode

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_write_particle_data (lttp, i_turn, bunch, map_with_rad)

type (ltt_params_struct) lttp
type (bunch_struct), target :: bunch, bunch_init
type (ptc_map_with_rad_struct) map_with_rad
type (coord_struct), pointer :: p, p0

integer i_turn, ix, ip, j
integer, save :: iu_snap = 0

character(200) file_name
character(40) fmt

!

if (lttp%particle_output_file == '') return
if (lttp%output_every_n_turns == -1 .and. i_turn /= lttp%n_turns) return
if (lttp%output_every_n_turns == 0 .and. i_turn /= 0 .and. i_turn /= lttp%n_turns) return
if (lttp%output_every_n_turns > 0 .and. modulo(i_turn, lttp%output_every_n_turns) /= 0) return

if (iu_snap == 0) then
  if (lttp%merge_particle_output_files) then
    file_name = lttp%particle_output_file
  else
    j = int(log10(real(lttp%n_turns, rp)) + 1 + 1d-10)
    write (fmt, '(a, i0, a, i0, a)') '(a, i', j, '.', j, ', a)'
    ix = index(lttp%particle_output_file, '#')
    if (ix == 0) then
      write (file_name, fmt) trim(lttp%particle_output_file), i_turn
    else
      write (file_name, fmt) lttp%particle_output_file(1:ix-1), i_turn, trim(lttp%particle_output_file(ix+1:))
    endif
  endif

  iu_snap = lunget()
  open (iu_snap, file = file_name, recl = 300)
  call ltt_write_this_header(lttp, iu_snap, map_with_rad, size(bunch%particle))

  if (lttp%output_initial_position) then
    write (iu_snap, '(2a)') '#      Ix     Turn | Start:    x              px               y              py               z              pz         spin_x    spin_y    spin_z  ', &
                                            '| End:         x              px               y              py               z              pz         spin_x    spin_y    spin_z    State'
  else
    write (iu_snap, '(a)')  '#      Ix     Turn |           x              px               y              py               z              pz   |     spin_x    spin_y    spin_z    State'
  endif
endif

do ip = 1, size(bunch%particle)
  p0 => bunch_init%particle(ip)
  p => bunch%particle(ip)
  if (lttp%output_initial_position) then
    write (iu_snap, '(i9, i9, 2(6es16.8, 3x, 3f10.6, 4x), a)') ip, i_turn, p0%vec, p0%spin, p%vec, p%spin, trim(coord_state_name(p%state))
  else
    write (iu_snap, '(i9, i9, 6es16.8, 3x, 3f10.6, 4x, a)')  ip, i_turn, p%vec, p%spin, trim(coord_state_name(p%state))
  endif
enddo

if (.not. lttp%merge_particle_output_files) then
  close(iu_snap)
  iu_snap = 0
endif

end subroutine ltt_write_particle_data

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_write_this_header(lttp, iu, map_with_rad, n_particle)

type (ltt_params_struct) lttp
type (ptc_map_with_rad_struct) map_with_rad

integer iu, n_particle

!

write (iu,  '(3a)')      '# lattice = "', trim(lttp%lat_file), '"'
write (iu,  '(3a)')      '# simulation_mode = "', trim(lttp%simulation_mode), '"'
write (iu,  '(a, i8)')   '# n_particle             = ', n_particle
write (iu,  '(a, i8)')   '# n_turns                = ', lttp%n_turns
write (iu,  '(a, l1)')   '# Radiation_Damping      = ', bmad_com%radiation_damping_on
write (iu,  '(a, l1)')   '# Radiation_Fluctuations = ', bmad_com%radiation_fluctuations_on
write (iu,  '(a, l1)')   '# Spin_tracking_on       = ', bmad_com%spin_tracking_on
if (lttp%use_1_turn_map) then
  write (iu, '(3a)')     '# Map_file               = "', trim(lttp%map_file), '"'
  write (iu, '(a, i0)')  '# Map_order              = ', map_with_rad%map_order
else
  write (iu, '(a)') '#'
  write (iu, '(a)') '#'
endif
write (iu, '(a)') '#'

end subroutine ltt_write_this_header


!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_calc_bunch_sums (lttp, i_turn, bunch, sum_data_arr)

type (ltt_params_struct) lttp
type (bunch_struct) bunch
type (ltt_sum_data_struct), allocatable, target :: sum_data_arr(:)
type (ltt_sum_data_struct), pointer :: sd

integer i, ix, i_turn

!

if (lttp%averages_output_file == '') return

select case (lttp%output_every_n_turns)
case (-1)
  if (i_turn /= lttp%n_turns) return
  if (.not. allocated(sum_data_arr)) allocate (sum_data_arr(1))
  sd => sum_data_arr(1)

case (0)
  if (i_turn /= 0 .and. i_turn /= lttp%n_turns) return
  if (.not. allocated(sum_data_arr)) allocate (sum_data_arr(0:1))
  ix = i_turn/lttp%output_every_n_turns
  sd => sum_data_arr(ix)

case default
  if (modulo(i_turn, lttp%output_every_n_turns) /= 0) return
  if (.not. allocated(sum_data_arr)) allocate (sum_data_arr(0:lttp%n_turns / lttp%output_every_n_turns))
  ix = i_turn/lttp%output_every_n_turns
  sd => sum_data_arr(ix)
end select

!

sd%n_live = count(bunch%particle%state == alive$)
sd%i_turn = i_turn

do i = 1, 6
  sd%orb_sum(i) = sum(bunch%particle%vec(i), bunch%particle%state == alive$) 
enddo

do i = 1, 3
  sd%spin_sum(i) = sum(bunch%particle%spin(i), bunch%particle%state == alive$)
enddo

end subroutine ltt_calc_bunch_sums

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_write_bunch_averages (lttp, sum_data_arr)

type (ltt_params_struct) lttp
type (ltt_sum_data_struct), target :: sum_data_arr(:)
type (ltt_sum_data_struct), pointer :: sd

integer iu, ix

!

if (lttp%averages_output_file == '') return

print '(a)', 'Averages_output_file: ', trim(lttp%averages_output_file)
iu = lunget()
open (iu, file = lttp%averages_output_file, recl = 200)
write (iu, '(a1, a8, a9, a14, 3a14, 6a14)') '#', 'Turn', 'N_live', 'Polarization', &
                     '<Sx>', '<Sy>', '<Sz>', '<x>', '<px>', '<y>', '<py>', '<z>', '<pz>'

!

do ix = lbound(sum_data_arr, 1) , ubound(sum_data_arr, 1)
  sd => sum_data_arr(ix)
  if (sd%n_live == 0) exit
  write (iu, '(i9, i9, f14.9, 2x, 3f14.9, 2x, 6es14.6)') sd%i_turn, sd%n_live, &
                        norm2(sd%spin_sum/sd%n_live), sd%spin_sum/sd%n_live, sd%orb_sum/sd%n_live
enddo

close (iu)

end subroutine ltt_write_bunch_averages

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_write_sigma_file (lttp, n_sum, average, sigma)

type (ltt_params_struct) lttp

real(rp) average(6), sigma(6,6)
integer i, n_sum

!

if (n_sum == 0) then
  write (1, '(a, i0)') '# n_turn_sigma_average = ', lttp%n_turn_sigma_average
  write (1, '(a)') '# NO DATA TO AVERAGE OVER!'
  return
endif

!

average = average / n_sum
sigma = sigma / n_sum - outer_product(average, average)

open (1, file = lttp%sigma_matrix_output_file)

write (1, '(a, i0)') '# n_turn_sigma_average = ', lttp%n_turn_sigma_average
write (1, '(a)') '# Average:'
write (1, '(5x, 6es16.8)') average
write (1, *)
write (1, '(a)') '# Sigma:'
do i = 1, 6
  write (1, '(5x, 6es16.8)') sigma(i,:)
enddo

close (1)

print '(2a)', 'Sigma matrix data file: ', trim(lttp%sigma_matrix_output_file)

end subroutine ltt_write_sigma_file

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine print_mpi_info (lttp, line, do_print)

type (ltt_params_struct) lttp

character(*) line
character(20) dtime
logical, optional :: do_print

!

if (.not. logic_option(lttp%debug, do_print)) return

call date_and_time_stamp (dtime)
print '(a, 2x, i0, 2a)', dtime, lttp%mpi_rank, ': ', trim(line)

end subroutine print_mpi_info

end module
