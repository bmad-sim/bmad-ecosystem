!+
! Module lt_tracking_mod
!
! Routines used by the long_term_tracking program.
!-

module lt_tracking_mod

use beam_mod
use twiss_and_track_mod
use ptc_map_with_radiation_mod
use high_energy_space_charge_mod
use s_fitting_new, only: probe, internal_state, track_probe, assignment(=)

implicit none

integer, parameter :: master_rank$ = 0
integer, parameter :: results_tag$ = 1000
integer, parameter :: is_done_tag$ = 1001

type ltt_params_struct
  type (lat_ele_loc_struct) :: start
  character(20) :: simulation_mode = ''
  character(20) :: tracking_method = 'BMAD'
  character(40) :: ele_start_name = ''
  character(200) :: lat_file = ''
  character(200) :: common_master_input_file = ''
  character(200) :: particle_output_file = ''
  character(200) :: sigma_matrix_output_file = ''
  character(200) :: map_file = ''
  character(200) :: averages_output_file = ''
  integer :: n_turns = -1
  integer :: random_seed = 0
  integer :: map_order = -1
  integer :: n_turn_sigma_average = -1
  integer :: output_every_n_turns = -1
  real(rp) :: ptc_aperture(2) = 0.1
  real(rp) :: print_on_dead_loss = -1
  real(rp) :: timer_print_dtime = 120
  real(rp) :: dead_cutoff = 1
  real(rp) :: a_emittance = 0   ! Used for space charge calculation.
  real(rp) :: b_emittance = 0   ! Used for space charge calculation.
  logical :: rfcavity_on = .true.
  logical :: add_closed_orbit_to_init_position = .true.
  logical :: output_initial_position = .false.
  logical :: merge_particle_output_files = .false.
  integer :: mpi_rank  = master_rank$
  integer :: mpi_n_proc = 1                    ! Number of processeses including master
  integer :: mpi_num_runs = 10                 ! Number of runs a slave process will take on average.
  integer :: mpi_n_particles_per_run = 0       ! Number of particles per run.
  logical :: using_mpi = .false.
  logical :: debug = .false.
  logical :: need_map = .false.     ! Internal var
  logical :: write_map = .false.
end type

type ltt_internal_struct
  ! Internal vars
  type (internal_state) ptc_state
  type (coord_struct), allocatable :: bmad_closed_orb(:)
  real(rp) ptc_closed_orb(6)
end type

type ltt_sum_data_struct
  integer :: i_turn = 0
  integer :: n_live = 0
  real(rp) :: orb_sum(6) = 0    ! Orbit average
  real(rp) :: orb2_sum(6,6) = 0
  real(rp) :: spin_sum(3) = 0   ! Spin
end type

contains

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_init_params(ltt, lat, beam_init)

type (ltt_params_struct) ltt
type (lat_struct) lat
type (beam_init_struct) beam_init

integer ir, i
character(200) init_file, arg

namelist / params / bmad_com, beam_init, ltt

! Parse command line

init_file = ''

i = 0
do while (i < cesr_iargc())
  i = i + 1
  call cesr_getarg(i, arg)
  select case (arg)
  case ('-write')
    ltt%write_map = .true.
  case default
    if (init_file /= '') then
      print '(2a)', 'Extra stuff on the command line: ', quote(arg)
      print '(a)',  'Stopping here.'
      stop
    endif
    init_file = arg
  end select
end do

if (init_file == '') init_file = 'long_term_tracking.init'

! Read parameters

if (.not. ltt%using_mpi .or. ltt%mpi_rank == master_rank$) then
  print '(2a)', 'Initialization file: ', trim(init_file)
endif

open (1, file = init_file, status = 'old', action = 'read')
read (1, nml = params)
close (1)

call upcase_string(ltt%simulation_mode)

if (ltt%common_master_input_file /= '') then
  print '(2a)', 'Using common_master_input_file: ', trim(ltt%common_master_input_file)

  open (1, file = ltt%common_master_input_file, status = 'old', action = 'read')
  read (1, nml = params)
  close (1)

  open (1, file = init_file, status = 'old', action = 'read')
  read (1, nml = params)
  close (1)
endif

! Lattice init

bmad_com%auto_bookkeeper = .false.

call ran_seed_put (ltt%random_seed)
call ptc_ran_seed_put (ltt%random_seed)

if (ltt%using_mpi) then
  call ran_seed_get (ir)
  call ran_seed_put (ir + 10 * ltt%mpi_rank)
  call ptc_ran_seed_put (ir + 10 * ltt%mpi_rank)
endif

call bmad_parser (ltt%lat_file, lat)

! Read the master input file again so that bmad_com parameters set in the file
! take precedence over bmad_com parameters set in the lattice file.

open (1, file = init_file, status = 'old', action = 'read')
read (1, nml = params)  
close (1)

! PTC has an internal aperture of 1.0 meter. To be safe use an aperture of 0.9 meter

ltt%ptc_aperture = min([0.9_rp, 0.9_rp], ltt%ptc_aperture)

!

if (ltt%tracking_method == 'STANDARD') then
  ltt%tracking_method = 'BMAD'
  print *, 'Note: LTT%TRACKING_METHOD = "STANDARD" has been renamed "BMAD". Will run normally...'
endif

! map file logic

if (ltt%map_file(1:6) == 'WRITE:') then
  ltt%map_file = ltt%map_file(7:)
  ltt%write_map = .true.
endif

if (ltt%write_map) ltt%need_map = .true.
if (ltt%tracking_method == 'MAP') ltt%need_map = .true.
if (ltt%simulation_mode == 'CHECK') ltt%need_map = .true.
if (ltt%simulation_mode == 'STAT') ltt%need_map = .false.

if (ltt%need_map .and. ltt%map_file == '') write (ltt%map_file, '(a, i0, a)') 'Order', ltt%map_order, '.map'

end subroutine ltt_init_params

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_init_tracking(lttp, lat, ltt_internal, rad_map)

type (ltt_params_struct) lttp
type (lat_struct), target :: lat
type (ltt_internal_struct) ltt_internal
type (ptc_map_with_rad_struct) rad_map
type (ele_pointer_struct), allocatable :: eles(:)
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele_start

real(rp) time, closed_orb(6)

integer n_loc, ix_ele_end, ix_ele_start, ix_branch

logical err, map_file_exists

!

if (lttp%using_mpi .and. lttp%tracking_method == 'PTC') then
  print '(a)', 'TRACKING_METHOD = PTC NOT ALLOWED WITH MPI VERSION OF THIS PROGRAM.'
  stop
endif

!

if (lttp%ele_start_name == '') Then
  lttp%start = lat_ele_loc_struct(0, 0)
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
  lttp%start = lat_ele_loc_struct(eles(1)%ele%ix_ele, eles(1)%ele%ix_branch)
endif

ix_ele_start = lttp%start%ix_ele
ix_branch = lttp%start%ix_branch
ele_start => pointer_to_ele(lat, lttp%start)
branch => lat%branch(ix_branch)
ix_ele_end = ix_ele_start - 1
if (ix_ele_end == 0) ix_ele_end = branch%n_ele_track

call twiss_and_track (lat, ltt_internal%bmad_closed_orb, ix_branch = ix_branch)

! Init 1-turn map.

inquire (file = lttp%map_file, exist = map_file_exists)

if (.not. lttp%rfcavity_on) then
  if (lttp%need_map) then
    print '(a)', 'RF will not be turned OFF since 1-turn map is in use!'
  else
    if (.not. lttp%rfcavity_on) call set_on_off (rfcavity$, lat, off$)
  endif
endif

if (lttp%need_map .and. map_file_exists .and. .not. lttp%write_map) then
  call ptc_read_map_with_radiation(lttp%map_file, rad_map)

  if (lttp%mpi_rank == master_rank$) then
    print '(2a)',   'Map read in from file: ', trim(lttp%map_file)
    print '(2a)',   'Lattice file used for map:         ', trim(rad_map%lattice_file)
    if ((rad_map%radiation_damping_on .neqv. bmad_com%radiation_damping_on) .or. &
                                                           (rad_map%map_order /= lttp%map_order)) then
      print '(a)',  'Map in file does not have the same map order or radiation_damping_on setting as in input files.'
      print '(a)',  'Will make a new map...'
    else
      lttp%need_map = .false.
    endif

  else
    ! Wait until master has created proper map if needed
    do
      if ((rad_map%radiation_damping_on .eqv. bmad_com%radiation_damping_on) .and. &
                                                            (rad_map%map_order == lttp%map_order)) exit
      call milli_sleep(1000)
      call ptc_read_map_with_radiation(lttp%map_file, rad_map)
    enddo
  endif
endif

if (lttp%need_map) then
  print '(a)', 'Creating map...'
  call run_timer ('START')
  call lat_to_ptc_layout (lat)
  call ptc_setup_map_with_radiation (rad_map, ele_start, ele_start, lttp%map_order, bmad_com%radiation_damping_on)
  call run_timer ('READ', time)
  print '(a, f8.2)', 'Map setup time (min)', time/60
endif

if (lttp%write_map) then
  call ptc_write_map_with_radiation(lttp%map_file, rad_map)
  print '(2a)', 'Created map file: ', trim(lttp%map_file)
endif

! Init PTC layout if needed

if (lttp%tracking_method == 'PTC' .or. lttp%simulation_mode == 'CHECK') then
  if (.not. associated(lat%branch(0)%ptc%m_t_layout)) call lat_to_ptc_layout(lat)
  call ptc_setup_tracking_with_damping_and_excitation(lat%branch(0), bmad_com%radiation_damping_on, &
                                          bmad_com%radiation_fluctuations_on, ltt_internal%ptc_state, ltt_internal%ptc_closed_orb)
endif

end subroutine ltt_init_tracking

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_print_inital_info (lttp, rad_map)

type (ltt_params_struct) lttp
type (ptc_map_with_rad_struct) rad_map

! Print some info.

print '(a)', '--------------------------------------'
print '(a, a)',  'lattice:                  ', quote(lttp%lat_file)
print '(a, a)',  'Averages_output_file:     ', quote(lttp%averages_output_file)
print '(a, a)',  'sigma_matrix_output_file: ', quote(lttp%sigma_matrix_output_file)
print '(a, a)',  'particle_output_file:     ', quote(lttp%particle_output_file)
print '(a, a)',  'averages_output_file:     ', quote(lttp%averages_output_file)
print '(a, a)',  'map_file:                 ', quote(lttp%map_file)
print '(a, a)',  'simulation_mode: ', trim(lttp%simulation_mode)
print '(a, a)',  'tracking_method: ', trim(lttp%tracking_method)
if (lttp%tracking_method == 'MAP') print '(a, i0)', 'map_order: ', rad_map%map_order
print '(a, l1)', 'Radiation Damping:           ', bmad_com%radiation_damping_on
print '(a, l1)', 'Stochastic Fluctuations:     ', bmad_com%radiation_fluctuations_on
print '(a, l1)', 'Spin_tracking_on:            ', bmad_com%spin_tracking_on

print '(a)', '--------------------------------------'

end subroutine ltt_print_inital_info

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_run_check_mode (lttp, lat, rad_map, beam_init, ltt_internal)

type (ltt_params_struct) lttp
type (lat_struct), target :: lat
type (ptc_map_with_rad_struct) rad_map
type (beam_init_struct) beam_init
type (ltt_internal_struct) ltt_internal
type (coord_struct), allocatable :: orb(:)
type (coord_struct) orb_end, orb_init
type (ele_struct), pointer :: ele_start
type (probe) prb

real(rp) ptc_track, ptc_spin(3)
integer track_state, ix_ele_start, ix_branch

! Run serial in check mode.

ix_ele_start = lttp%start%ix_ele
ele_start => pointer_to_ele(lat, lttp%start)
ix_branch = ele_start%ix_branch

call init_coord (orb_init, beam_init%center, ele_start, downstream_end$, spin = beam_init%spin)
call reallocate_coord(orb, lat)

if (lttp%add_closed_orbit_to_init_position) orb_init%vec = orb_init%vec + ltt_internal%bmad_closed_orb(ix_ele_start)%vec

orb_end = orb_init
orb(ix_ele_start) = orb_init

call ptc_track_map_with_radiation (orb_end, rad_map)
call track_many (lat, orb, ix_ele_start, ix_ele_start, 1, ix_branch, track_state)

prb = orb_init%vec
prb%q%x = [1, 0, 0, 0]   ! Unit quaternion
call track_probe (prb, ltt_internal%ptc_state, fibre1 = lat%branch(ix_branch)%ele(1)%ptc_fibre)
ptc_spin = rotate_vec_given_quat(prb%q%x, orb_init%spin)

print '(a, 6f14.8)', 'Map closed orbit at start:  ', rad_map%sub_map(1)%fix0
print '(a, 6f14.8)', 'Bmad closed orbit at start: ', ltt_internal%bmad_closed_orb(ix_ele_start)%vec

print *
print '(a, 6f14.8)', 'Starting orbit for tracking:', orb_init%vec

print *
print '(a)', 'Phase Space at Track End (Bmad closed orbit used with map tracking:'
print '(a, 6f14.8)', 'Map tracking:   ', orb_end%vec
print '(a, 6f14.8)', 'Bmad tracking:  ', orb(ix_ele_start)%vec
print '(a, 6f14.8)', 'PTC tracking:   ', prb%x
print '(a, 6f14.8)', 'Diff Map - PTC:   ', orb_end%vec - prb%x
print '(a, 6f14.8)', 'Diff Bmad - PTC:  ', orb(ix_ele_start)%vec - prb%x
print '(a, 6f14.8)', 'Diff Bmad - Map:  ', orb(ix_ele_start)%vec - orb_end%vec

print *
print '(a, 6f14.8)', 'Initial spin: ', orb_init%spin

print *
print '(a)', 'Spin at Track End:'
print '(a, 6f14.8)', 'Map tracking: ', orb_end%spin
print '(a, 6f14.8)', 'Bmad tracking:', orb(ix_ele_start)%spin
print '(a, 6f14.8)', 'PTC tracking: ', ptc_spin
print '(a, 6f14.8)', 'Diff Map - PTC: ', orb_end%spin - ptc_spin
print '(a, 6f14.8)', 'Diff Bmad - PTC:', orb(ix_ele_start)%spin - ptc_spin
print '(a, 6f14.8)', 'Diff Bmad - Map:', orb(ix_ele_start)%spin - orb_end%spin

end subroutine ltt_run_check_mode

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_run_single_mode (lttp, lat, beam_init, ltt_internal,rad_map)

type (ltt_params_struct) lttp
type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (beam_init_struct) beam_init
type (ltt_internal_struct) ltt_internal
type (coord_struct), allocatable :: orb(:)
type (coord_struct) :: orbit, orbit_old
type (ele_struct), pointer :: ele
type (ptc_map_with_rad_struct) rad_map
type (ele_struct), pointer :: ele_start
type (probe) prb

real(rp) average(6), sigma(6,6)
integer n_sum, ix_ele_start, iu_out, i_turn, track_state, ix_branch
logical is_lost

character(40) fmt

! Run serial in single mode.

n_sum = 0
sigma = 0
average = 0
ix_ele_start = lttp%start%ix_ele
ix_branch = lttp%start%ix_branch
ele_start => pointer_to_ele(lat, lttp%start)
branch => lat%branch(ix_branch)
orbit_old%vec = real_garbage$

call ltt_setup_high_energy_space_charge(lttp, branch, beam_init, ltt_internal)

if (lttp%tracking_method == 'BMAD') call reallocate_coord (orb, lat)
call init_coord (orbit, beam_init%center, ele_start, downstream_end$, lat%param%particle, spin = beam_init%spin)

if (lttp%add_closed_orbit_to_init_position) then
  select case (lttp%tracking_method)
  case ('MAP');       orbit%vec = orbit%vec + rad_map%sub_map(1)%fix0
  case ('BMAD');      orbit%vec = orbit%vec + ltt_internal%bmad_closed_orb(ix_ele_start)%vec
  case ('PTC');       orbit%vec = orbit%vec + ltt_internal%ptc_closed_orb
  end select
endif

fmt = '(i6, 6es16.8, 3x, 3f10.6)'
iu_out = lunget()
if (lttp%particle_output_file == '') lttp%particle_output_file = 'single.dat'
open(iu_out, file = lttp%particle_output_file, recl = 200)
call ltt_write_this_header(lttp, iu_out, rad_map, 1)
write (iu_out, '(a)') '# Turn |            x              px               y              py               z              pz    |   spin_x    spin_y    spin_z'
write (iu_out, fmt) 0, orbit%vec, orbit%spin

do i_turn = 1, lttp%n_turns
  select case (lttp%tracking_method)
  case ('MAP')
    call ptc_track_map_with_radiation (orbit, rad_map)
    is_lost = orbit_too_large(orbit)
  case ('BMAD')
    orb(ix_ele_start) = orbit
    call track_many (lat, orb, ix_ele_start, ix_ele_start, 1, ix_branch, track_state)
    is_lost = (track_state /= moving_forward$)
    orbit = orb(ix_ele_start)
  case ('PTC')
    prb = orbit%vec
    prb%q%x = [1, 0, 0, 0]  ! Unit quaternion
    call track_probe (prb, ltt_internal%ptc_state, fibre1 = lat%branch(ix_branch)%ele(1)%ptc_fibre)
    orbit%vec = prb%x
    orbit%spin = rotate_vec_given_quat(prb%q%x, orbit%spin)
    is_lost = (abs(orbit%vec(1)) > lttp%ptc_aperture(1) .or. abs(orbit%vec(3)) > lttp%ptc_aperture(2))
    if (all(orbit%vec == orbit_old%vec)) is_lost = .true.
    orbit_old%vec = orbit%vec
  case default
    print *, 'Unknown tracking_method: ' // lttp%tracking_method
    stop
  end select

  if (lttp%output_every_n_turns > 0 .and. modulo(i_turn, lttp%output_every_n_turns) == 0) then
    write (iu_out, fmt) i_turn, orbit%vec, orbit%spin
  endif

  if (is_lost) then
    ele => branch%ele(track_state)
    print '(a, i0, 8a)', 'Particle lost at turn: ', i_turn
    if (lttp%tracking_method == 'BMAD') print '(5a)', 'Lost at element: ', trim(ele%name), ' (', ele_loc_name(ele), '), State: ', coord_state_name(orb(track_state)%state)
    exit
  endif

  if (lttp%sigma_matrix_output_file /= '') then
    average = average + orbit%vec
    sigma = sigma + outer_product(orbit%vec, orbit%vec)
    n_sum = n_sum + 1
  endif
enddo

print '(2a)', 'Particle output file: ', trim(lttp%particle_output_file)
close (iu_out)

if (lttp%sigma_matrix_output_file /= '') then
  call ltt_write_single_mode_sigma_file (lttp, n_sum, average, sigma)
endif

end subroutine ltt_run_single_mode

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_run_bunch_mode (lttp, lat, beam_init, ltt_internal,rad_map, sum_data_array)

type (ltt_params_struct) lttp
type (lat_struct) lat
type (beam_init_struct) beam_init
type (ltt_internal_struct) ltt_internal
type (coord_struct), allocatable :: orb(:)
type (ptc_map_with_rad_struct) rad_map
type (bunch_struct), target :: bunch, bunch_init
type (coord_struct), pointer :: p
type (ltt_sum_data_struct), allocatable, optional :: sum_data_array(:)
type (ltt_sum_data_struct), allocatable, target :: sum_data_arr(:)
type (ltt_sum_data_struct), pointer :: sd
type (ele_struct), pointer :: ele_start
type (probe) prb

real(rp) time0, time

integer n, ix_ele_start, n_print_dead_loss, i_turn
integer ip, n_live_old, n_live, ix_branch

logical err_flag

! Init

ix_ele_start = lttp%start%ix_ele
ele_start => pointer_to_ele (lat, lttp%start)
ix_branch = ele_start%ix_branch

if (lttp%using_mpi) beam_init%n_particle = lttp%mpi_n_particles_per_run
call init_bunch_distribution (ele_start, lat%param, beam_init, lttp%start%ix_branch, bunch, err_flag)
if (err_flag) stop

call ltt_setup_high_energy_space_charge(lttp, lat%branch(ix_branch), beam_init, ltt_internal)

if (lttp%mpi_rank == master_rank$) print '(a, i8)',   'n_particle             = ', size(bunch%particle)

do n = 1, size(bunch%particle)
  p => bunch%particle(n)
  if (lttp%add_closed_orbit_to_init_position) then
    select case (lttp%tracking_method)
    case ('MAP');       p%vec = p%vec + rad_map%sub_map(1)%fix0
    case ('BMAD');      p%vec = p%vec + ltt_internal%bmad_closed_orb(ix_ele_start)%vec
    case ('PTC');       p%vec = p%vec + ltt_internal%ptc_closed_orb
    end select
  endif
enddo

n_print_dead_loss = max(1, nint(lttp%print_on_dead_loss * size(bunch%particle)))
allocate(bunch_init%particle(size(bunch%particle)))
bunch_init%particle = bunch%particle
n_live_old = size(bunch%particle)

! 

if (.not. lttp%using_mpi) call ltt_write_particle_data (lttp, 0, bunch, bunch_init, rad_map)
call ltt_calc_bunch_sums (lttp, 0, bunch, sum_data_arr)

time0 = 0

do i_turn = 1, lttp%n_turns
  select case (lttp%tracking_method)
  case ('MAP')
    do ip = 1, size(bunch%particle)
      p => bunch%particle(ip)
      if (p%state /= alive$) cycle
      call ptc_track_map_with_radiation (p, rad_map)
      if (abs(p%vec(1)) > lttp%ptc_aperture(1) .or. abs(p%vec(3)) > lttp%ptc_aperture(2)) p%state = lost$
    enddo

  case ('BMAD')
    call track_bunch (lat, bunch, ele_start, ele_start, err_flag)

  case ('PTC')
    do ip = 1, size(bunch%particle)
      p => bunch%particle(ip)
      if (p%state /= alive$) cycle
      prb = p%vec
      prb%q%x = [1, 0, 0, 0]  ! Unit quaternion
      call track_probe (prb, ltt_internal%ptc_state, fibre1 = lat%branch(ix_branch)%ele(1)%ptc_fibre)
      p%vec = prb%x
      p%spin = rotate_vec_given_quat(prb%q%x, p%spin)
      if (abs(p%vec(1)) > lttp%ptc_aperture(1) .or. abs(p%vec(3)) > lttp%ptc_aperture(2)) p%state = lost$
    enddo

  case default
    print *, 'Unknown tracking_method: ' // lttp%tracking_method
    stop
  end select

  n_live = count(bunch%particle%state == alive$)
  if (n_live_old - n_live >= n_print_dead_loss) then
    print '(a, i0, a, i0)', 'Total number dead on turn ', i_turn, ': ', size(bunch%particle) - n_live
    n_live_old = n_live
  endif

  if (size(bunch%particle) - n_live > lttp%dead_cutoff * size(bunch%particle) .and. .not. lttp%using_mpi) then
    print '(a)', 'Particle loss greater than set by dead_cutoff. Stopping now.'
    exit
  endif

  call run_timer('READ', time)

  if (time-time0 > lttp%timer_print_dtime) then
    print '(a, f10.2, a, i0)', 'Ellapsed time (min): ', time/60, ', At turn: ', i_turn
    time0 = time
  endif

  if (.not. lttp%using_mpi) call ltt_write_particle_data (lttp, i_turn, bunch, bunch_init, rad_map)
  call ltt_calc_bunch_sums (lttp, i_turn, bunch, sum_data_arr)

  if (.not. lttp%using_mpi) then
    call ltt_write_bunch_averages(lttp, i_turn, sum_data_arr)
    call ltt_write_sigma_matrix (lttp, i_turn, sum_data_arr)
  endif

end do

if (lttp%mpi_rank == master_rank$) print '(2a)', 'Tracking data file: ', trim(lttp%particle_output_file)

if (lttp%using_mpi) then
  if (allocated(sum_data_array)) deallocate (sum_data_array)
  call move_alloc (sum_data_arr, sum_data_array)
endif

end subroutine ltt_run_bunch_mode

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_run_stat_mode (lttp, lat, ltt_internal)

type (ltt_params_struct) lttp
type (lat_struct), target :: lat
type (ltt_internal_struct) ltt_internal
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele
type (normal_modes_struct) modes
type (rad_int_all_ele_struct) rad_int_ele

real(rp) chrom_x, chrom_y, ring_length
integer i
logical err_flag

! Run serial in stat mode.

branch => lat%branch(lttp%start%ix_branch)

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
  write(22,'(a,1x,f10.3,4(1x,f10.6))') ele%name, ele%s, ltt_internal%bmad_closed_orb(i)%vec(1:5:2), ltt_internal%bmad_closed_orb(i)%vec(6)
enddo

close(20)
close(21)
close(22)

!

ring_length = branch%param%total_length
call chrom_calc(lat, 1.0d-6, chrom_x, chrom_y, err_flag, ix_branch = branch%ix_branch)
call calc_z_tune (lat)
call radiation_integrals (lat, ltt_internal%bmad_closed_orb, modes, rad_int_by_ele = rad_int_ele, ix_branch = branch%ix_branch)

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

subroutine ltt_setup_high_energy_space_charge(lttp, branch, beam_init, ltt_internal)

type (ltt_params_struct) lttp
type (branch_struct) branch
type (beam_init_struct) beam_init
type (ltt_internal_struct) ltt_internal
type (normal_modes_struct) mode
real(rp) n_particle

!

if (.not. branch%param%high_energy_space_charge_on) return

if (lttp%tracking_method == 'MAP') then
  print *, 'NOTE: Space effects are not present when using a 1-turn map for tracking!'
  return
endif

if (.not. lttp%rfcavity_on) then
  print *, 'NOTE: RF is not on. Cannot calculate a longitudinal bunch length.'
  print *, '      Therefore no space charge kick will be applied.'
  return
endif

call radiation_integrals(branch%lat, ltt_internal%bmad_closed_orb, mode, ix_branch = branch%ix_branch)
if (lttp%a_emittance /= 0) mode%a%emittance = lttp%a_emittance
if (lttp%b_emittance /= 0) mode%b%emittance = lttp%b_emittance
n_particle = abs(beam_init%bunch_charge / (e_charge * charge_of(ltt_internal%bmad_closed_orb(0)%species)))
call setup_high_energy_space_charge_calc (.true., branch, n_particle, mode)

end subroutine ltt_setup_high_energy_space_charge

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_write_particle_data (lttp, i_turn, bunch, bunch_init, rad_map)

type (ltt_params_struct) lttp
type (bunch_struct), target :: bunch, bunch_init
type (ptc_map_with_rad_struct) rad_map
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
  call ltt_write_this_header(lttp, iu_snap, rad_map, size(bunch%particle))

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

subroutine ltt_write_this_header(lttp, iu, rad_map, n_particle)

type (ltt_params_struct) lttp
type (ptc_map_with_rad_struct) rad_map

integer iu, n_particle

!

write (iu,  '(3a)')      '# lattice = "', trim(lttp%lat_file), '"'
write (iu,  '(3a)')      '# simulation_mode = "', trim(lttp%simulation_mode), '"'
write (iu,  '(a, i8)')   '# n_particle             = ', n_particle
write (iu,  '(a, i8)')   '# n_turns                = ', lttp%n_turns
write (iu,  '(a, l1)')   '# Radiation_Damping      = ', bmad_com%radiation_damping_on
write (iu,  '(a, l1)')   '# Radiation_Fluctuations = ', bmad_com%radiation_fluctuations_on
write (iu,  '(a, l1)')   '# Spin_tracking_on       = ', bmad_com%spin_tracking_on
if (lttp%tracking_method == 'MAP') then
  write (iu, '(3a)')     '# Map_file               = "', trim(lttp%map_file), '"'
  write (iu, '(a, i0)')  '# Map_order              = ', rad_map%map_order
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

integer i, j, ix, i_turn

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
  do j = i, 6
    sd%orb2_sum(i,j) = sum(bunch%particle%vec(i) * bunch%particle%vec(j), bunch%particle%state == alive$) 
  enddo
enddo

do i = 1, 3
  sd%spin_sum(i) = sum(bunch%particle%spin(i), bunch%particle%state == alive$)
enddo

end subroutine ltt_calc_bunch_sums

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_write_bunch_averages (lttp, i_turn, sum_data_arr)

type (ltt_params_struct) lttp
type (ltt_sum_data_struct), target :: sum_data_arr(:)
type (ltt_sum_data_struct), pointer :: sd

real(rp) sigma(6)
integer i_turn, i, iu, ix, ix_last

logical new

! i_turn = -1 is used with mpi version.

if (lttp%averages_output_file == '') return

if (i_turn == -1) then
  new = .true.
  ix_last = -1
else
  select case (lttp%output_every_n_turns)
  case (-1)
    if (i_turn /= lttp%n_turns) return
    new = .true.
    ix_last = 1
  case (0)
    if (i_turn /= 0 .and. i_turn /= lttp%n_turns) return
    ix_last = i_turn/lttp%output_every_n_turns
    new = (i_turn == 0)
  case default
    if (modulo(i_turn, lttp%output_every_n_turns) /= 0) return
    ix_last = i_turn/lttp%output_every_n_turns
    new = (i_turn == lttp%output_every_n_turns)
  end select
endif 

!

iu = lunget()

if (new) then
  open (iu, file = lttp%averages_output_file, recl = 300)
  write (iu, '(a1, a8, a9, a14, 3a14, 12a14)') '#', 'Turn', 'N_live', 'Polarization', &
                     '<Sx>', '<Sy>', '<Sz>', '<x>', '<px>', '<y>', '<py>', '<z>', '<pz>', &
                     'Sig_x', 'Sig_px', 'Sig_y', 'Sig_py', 'Sig_z', 'Sig_pz'
else
  open (iu, file = lttp%averages_output_file, recl = 300, access = 'append')
endif

! If i_turn = -1 (mpi) then dump all the averages.

do ix = lbound(sum_data_arr, 1), ubound(sum_data_arr, 1)
  if (i_turn /= -1 .and. ix /= ix_last) cycle
  sd => sum_data_arr(ix)
  if (sd%n_live == 0) exit
  forall (i = 1:6) sigma(i) = sd%orb2_sum(i,i)/sd%n_live - (sd%orb_sum(i)/sd%n_live)**2
  sigma = sqrt(max(0.0_rp, sigma))
  write (iu, '(i9, i9, f14.9, 2x, 3f14.9, 2x, 12es14.6)') sd%i_turn, sd%n_live, &
          norm2(sd%spin_sum/sd%n_live), sd%spin_sum/sd%n_live, sd%orb_sum/sd%n_live, sigma
enddo

!

close (iu)

end subroutine ltt_write_bunch_averages

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_write_sigma_matrix (lttp, i_turn, sum_data_arr)

type (ltt_params_struct) lttp
type (ltt_sum_data_struct), target :: sum_data_arr(:)
type (ltt_sum_data_struct), pointer :: sd

real(rp) sigma(21)
integer i_turn, ix, i, j, k, iu, ix_last
logical new

! i_turn = -1 is used with mpi version.

if (lttp%sigma_matrix_output_file == '') return

if (i_turn == -1) then
  new = .true.
  ix_last = -1
else
  select case (lttp%output_every_n_turns)
  case (-1)
    if (i_turn /= lttp%n_turns) return
    new = .true.
    ix_last = 1
  case (0)
    if (i_turn /= 0 .and. i_turn /= lttp%n_turns) return
    ix_last = i_turn/lttp%output_every_n_turns
    new = (i_turn == 0)
  case default
    if (modulo(i_turn, lttp%output_every_n_turns) /= 0) return
    ix_last = i_turn/lttp%output_every_n_turns
    new = (i_turn == lttp%output_every_n_turns)
  end select
endif 

!

iu = lunget()

if (new) then
  open (iu, file = lttp%sigma_matrix_output_file, recl = 400)
  write (iu, '(a1, a8, a9, 21a14)') '#', 'Turn', 'N_live', &
    '<x.x>', '<x.px>', '<x.y>', '<x.py>', '<x.z>', '<x.pz>', '<px.px>', '<px.y>', '<px.py>', '<px.z>', '<px.pz>', &
    '<y.y>', '<y.py>', '<y.z>', '<y.pz>', '<py.py>', '<py.z>', '<py.pz>', '<z.z>', '<z.pz>', '<pz.pz>'
else
  open (iu, file = lttp%sigma_matrix_output_file, recl = 400, access = 'append')
endif

! If i_turn = -1 (mpi) then dump all the averages.

do ix = lbound(sum_data_arr, 1), ubound(sum_data_arr, 1)
  if (i_turn /= -1 .and. ix /= ix_last) cycle
  sd => sum_data_arr(ix)
  if (sd%n_live == 0) exit
  k = 0
  do i = 1, 6
  do j = i, 6
    k = k + 1
    sigma(k) = sd%orb2_sum(i,j) / sd%n_live - sd%orb_sum(i) * sd%orb_sum(j) / sd%n_live**2
  enddo
  enddo
  write (iu, '(i9, i9, 2x, 21es14.6)') sd%i_turn, sd%n_live, (sigma(k), k = 1, 21)
enddo

!

close (iu)

end subroutine ltt_write_sigma_matrix

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_write_single_mode_sigma_file (lttp, n_sum, average, sigma)

type (ltt_params_struct) lttp

real(rp) average(6), sigma(6,6)
integer i, n_sum

!

if (lttp%sigma_matrix_output_file == '') return

open (1, file = lttp%sigma_matrix_output_file)

if (n_sum == 0) then
  write (1, '(a, i0)') '# n_turn_sigma_average = ', lttp%n_turn_sigma_average
  write (1, '(a)') '# NO DATA TO AVERAGE OVER!'
  return
endif

!

average = average / n_sum
sigma = sigma / n_sum - outer_product(average, average)

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

end subroutine ltt_write_single_mode_sigma_file

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
