!+
! Module lt_tracking_mod
!
! Routines used by the long_term_tracking program.
!-

module lt_tracking_mod

use beam_mod
use twiss_and_track_mod
use ptc_map_with_radiation_mod
use space_charge_mod

implicit none

integer, parameter :: master_rank$ = 0
integer, parameter :: results_tag$ = 1000
integer, parameter :: is_done_tag$ = 1001

type ltt_params_struct
  type (lat_ele_loc_struct) :: start
  character(20) :: simulation_mode = ''
  character(20) :: tracking_method = 'STANDARD'
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
  real(rp) :: print_on_dead_loss = -1
  real(rp) :: timer_print_dtime = 120
  real(rp) :: dead_cutoff = 0
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
call ptc_set_map_with_radiation_ran_seed(ltt%random_seed)

if (ltt%using_mpi) then
  call ran_seed_get (ir)
  call ran_seed_put (ir + 10 * ltt%mpi_rank)
  call ptc_set_map_with_radiation_ran_seed(ir + 10 * ltt%mpi_rank)
endif

call bmad_parser (ltt%lat_file, lat)

! Read the master input file again so that bmad_com parameters set in the file
! take precedence over bmad_com parameters set in the lattice file.

open (1, file = init_file, status = 'old', action = 'read')
read (1, nml = params)  
close (1)

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

subroutine ltt_init_tracking(lttp, lat, closed_orb, rad_map)

type (ltt_params_struct) lttp
type (lat_struct), target :: lat
type (ptc_map_with_rad_struct) rad_map
type (ele_pointer_struct), allocatable :: eles(:)
type (coord_struct), allocatable :: closed_orb(:), orb(:)
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele_start

real(rp) time

integer n_loc, ix_ele_end, ix_ele_start, ix_branch

logical err, map_file_exists

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

call twiss_and_track (lat, closed_orb, ix_branch = ix_branch)

! Init 1-turn map.

inquire (file = lttp%map_file, exist = map_file_exists)

if (.not. lttp%rfcavity_on) then
  if (lttp%need_map) then
    print '(a)', 'RF will not be turned OFF since 1-turn map is in use!'
  else
    if (.not. lttp%rfcavity_on) call set_on_off (rfcavity$, lat, off$)
  endif
endif

if (lttp%need_map .and. map_file_exists) then
  call ptc_read_map_with_radiation(lttp%map_file, rad_map)

  if (lttp%mpi_rank == master_rank$) then
    print '(2a)',   'Map read in from file: ', trim(lttp%map_file)
    print '(2a)',   'Lattice file used for map:         ', trim(rad_map%lattice_file)
    if (rad_map%radiation_damping_on /= bmad_com%radiation_damping_on .or. rad_map%map_order /= lttp%map_order) then
      print '(a)',  'Map in file does not have the same map order or radiation_damping_on setting as in input files.'
      print '(a)',  'Will make a new map...'
    else
      lttp%need_map = .false.
    endif

  else
    ! Wait until master has created proper map if needed
    do
      if (rad_map%radiation_damping_on == bmad_com%radiation_damping_on .and. rad_map%map_order == lttp%map_order) exit
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

end subroutine ltt_init_tracking

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_print_inital_info (lttp, rad_map)

type (ltt_params_struct) lttp
type (ptc_map_with_rad_struct) rad_map

! Print some info.

print '(a)', '--------------------------------------'
print '(a, a)',  'lattice:         ', trim(lttp%lat_file)
print '(a, a)',  'simulation_mode: ', lttp%simulation_mode
print '(a, l1)', 'Radiation Damping:           ', bmad_com%radiation_damping_on
print '(a, l1)', 'Stochastic Fluctuations:     ', bmad_com%radiation_fluctuations_on
print '(a, l1)', 'Spin_tracking_on:            ', bmad_com%spin_tracking_on
if (lttp%need_map) then
  print '(a, i0)', 'Map order:          ', rad_map%map_order
endif
print '(a)', '--------------------------------------'

end subroutine ltt_print_inital_info

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_run_check_mode (lttp, lat, rad_map, beam_init, closed_orb)

type (ltt_params_struct) lttp
type (lat_struct), target :: lat
type (ptc_map_with_rad_struct) rad_map
type (beam_init_struct) beam_init
type (coord_struct), allocatable :: closed_orb(:), orb(:)
type (coord_struct) orb_end, orb_init
type (ele_struct), pointer :: ele_start

integer track_state, ix_ele_start

! Run serial in check mode.

ix_ele_start = lttp%start%ix_ele
ele_start => pointer_to_ele(lat, lttp%start)

call init_coord (orb_init, beam_init%center, ele_start, downstream_end$, spin = beam_init%spin)
call reallocate_coord(orb, lat)

if (lttp%add_closed_orbit_to_init_position) orb_init%vec = orb_init%vec + closed_orb(ix_ele_start)%vec

orb_end = orb_init
orb(ix_ele_start) = orb_init

call ptc_track_map_with_radiation (orb_end, rad_map)
call track_many (lat, orb, ix_ele_start, ix_ele_start, 1, lttp%start%ix_branch, track_state)

print '(a, 6f14.8)', 'Map closed orbit at start: ', rad_map%sub_map(1)%fix0
print '(a, 6f14.8)', 'Bmad closed orbit at start:', closed_orb(ix_ele_start)%vec

print *
print '(a)', 'Phase Space at Track End (Bmad closed orbit used with map tracking:'
print '(a, 6f14.8)', 'Map tracking:       ', orb_end%vec
print '(a, 6f14.8)', 'Ele-by-ele tracking:', orb(ix_ele_start)%vec
print '(a, 6f14.8)', 'Diff:               ', orb_end%vec - orb(ix_ele_start)%vec

print *
print '(a)', 'Spin at Track End:'
print '(a, 6f14.8)', 'Map tracking:       ', orb_end%spin
print '(a, 6f14.8)', 'Ele-by-ele tracking:', orb(ix_ele_start)%spin
print '(a, 6f14.8)', 'Diff:               ', orb_end%spin - orb(ix_ele_start)%spin

end subroutine ltt_run_check_mode

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_run_single_mode (lttp, lat, beam_init, closed_orb, rad_map)

type (ltt_params_struct) lttp
type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (beam_init_struct) beam_init
type (coord_struct), allocatable :: closed_orb(:), orb(:)
type (coord_struct) orb_end
type (ele_struct), pointer :: ele
type (ptc_map_with_rad_struct) rad_map
type (ele_struct), pointer :: ele_start

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

call ltt_setup_space_charge(lttp, lat, beam_init, closed_orb)

call reallocate_coord (orb, lat)
call init_coord (orb(ix_ele_start), beam_init%center, ele_start, downstream_end$, lat%param%particle, spin = beam_init%spin)

if (lttp%add_closed_orbit_to_init_position) then
  select case (lttp%tracking_method)
  case ('MAP');       orb(ix_ele_start)%vec = orb(ix_ele_start)%vec + rad_map%sub_map(1)%fix0
  case ('STANDARD');  orb(ix_ele_start)%vec = orb(ix_ele_start)%vec + closed_orb(ix_ele_start)%vec
  case ('PTC');       call err_exit  ! Not yet implemented 
  end select
endif

fmt = '(i6, 6es16.8, 3x, 3f10.6)'
iu_out = lunget()
if (lttp%particle_output_file == '') lttp%particle_output_file = 'single.dat'
open(iu_out, file = lttp%particle_output_file, recl = 200)
call ltt_write_this_header(lttp, iu_out, rad_map, 1)
write (iu_out, '(a)') '# Turn |            x              px               y              py               z              pz    |   spin_x    spin_y    spin_z'
write (iu_out, fmt) 0, orb(ix_ele_start)%vec, orb(ix_ele_start)%spin

do i_turn = 1, lttp%n_turns
  select case (lttp%tracking_method)
  case ('MAP')
    orb_end = orb(ix_ele_start)
    call ptc_track_map_with_radiation (orb_end, rad_map)
    is_lost = orbit_too_large(orb_end)
  case ('STANDARD')
    call track_many (lat, orb, ix_ele_start, ix_ele_start, 1, ix_branch, track_state)
    is_lost = (track_state /= moving_forward$)
    orb_end = orb(ix_ele_start)
  case ('PTC')
    print *, 'Not yet implemented...'
    stop
  case default
    print *, 'Unknown tracking_method: ' // lttp%tracking_method
    stop
  end select

  if (lttp%output_every_n_turns > 0 .and. modulo(i_turn, lttp%output_every_n_turns) == 0) then
    write (iu_out, fmt) i_turn, orb_end%vec, orb_end%spin
  endif

  if (is_lost) then
    ele => branch%ele(track_state)
    print '(a, i0, 8a)', 'Particle lost at turn: ', i_turn
    if (lttp%tracking_method /= 'MAP') print '(5a)', 'Lost at element: ', trim(ele%name), ' (', ele_location(ele), '), State: ', coord_state_name(orb(track_state)%state)
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

subroutine ltt_run_bunch_mode (lttp, lat, beam_init, closed_orb, rad_map, sum_data_array)

type (ltt_params_struct) lttp
type (lat_struct) lat
type (beam_init_struct) beam_init
type (coord_struct), allocatable :: closed_orb(:), orb(:)
type (ptc_map_with_rad_struct) rad_map
type (bunch_struct), target :: bunch, bunch_init
type (coord_struct), pointer :: p
type (ltt_sum_data_struct), allocatable, optional :: sum_data_array(:)
type (ltt_sum_data_struct), allocatable, target :: sum_data_arr(:)
type (ltt_sum_data_struct), pointer :: sd
type (ele_struct), pointer :: ele_start

real(rp) time0, time
real(rp) average(6), sigma(6,6)

integer n, n_sum, ix_ele_start, n_print_dead_loss, i_turn
integer ip, n_live_old, n_live

logical err_flag

! Init

n_sum = 0
sigma = 0
average = 0
ix_ele_start = lttp%start%ix_ele
ele_start => pointer_to_ele (lat, lttp%start)

if (lttp%using_mpi) beam_init%n_particle = lttp%mpi_n_particles_per_run
call init_bunch_distribution (ele_start, lat%param, beam_init, lttp%start%ix_branch, bunch, err_flag)
if (err_flag) stop

call ltt_setup_space_charge(lttp, lat, beam_init, closed_orb)

if (lttp%mpi_rank == master_rank$) print '(a, i8)',   'n_particle             = ', size(bunch%particle)

do n = 1, size(bunch%particle)
  p => bunch%particle(n)
  if (lttp%add_closed_orbit_to_init_position) then
    select case (lttp%tracking_method)
    case ('MAP');       p%vec = p%vec + rad_map%sub_map(1)%fix0
    case ('STANDARD');  p%vec = p%vec + closed_orb(ix_ele_start)%vec
    case ('PTC');       call err_exit  ! Not yet implemented 
    end select
  endif
enddo

n_print_dead_loss = max(1, nint(lttp%print_on_dead_loss * size(bunch%particle)))
allocate(bunch_init%particle(size(bunch%particle)))
bunch_init%particle = bunch%particle
n_live_old = size(bunch%particle)

! 

if (.not. lttp%using_mpi) call ltt_write_particle_data (lttp, 0, bunch, rad_map)
call ltt_calc_bunch_sums (lttp, 0, bunch, sum_data_arr)

time0 = 0

do i_turn = 1, lttp%n_turns
  select case (lttp%tracking_method)
  case ('MAP')
    do ip = 1, size(bunch%particle)
      call ptc_track_map_with_radiation (bunch%particle(ip), rad_map)
      if (lttp%debug .and. ip == 1) print '(i6, 6f15.9)', i_turn, bunch%particle(ip)%vec
    enddo
  case ('STANDARD')
    call track_bunch (lat, bunch, ele_start, ele_start, err_flag)
  case ('PTC')
    print *, 'Not yet implemented...'
    stop
  case default
    print *, 'Unknown tracking_method: ' // lttp%tracking_method
    stop
  end select

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

  if (.not. lttp%using_mpi) call ltt_write_particle_data (lttp, i_turn, bunch, rad_map)
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

subroutine ltt_setup_space_charge(lttp, lat, beam_init, closed_orb)

type (ltt_params_struct) lttp
type (lat_struct) lat
type (beam_init_struct) beam_init
type (coord_struct) :: closed_orb(0:)
type (normal_modes_struct) mode
real(rp) n_particle

!

if (.not. bmad_com%space_charge_on) return

if (lttp%tracking_method == 'MAP') then
  print *, 'NOTE: Space effects are not present when using a 1-turn map for tracking!'
  return
endif

if (.not. lttp%rfcavity_on) then
  print *, 'NOTE: RF is not on. Cannot calculate a longitudinal bunch length.'
  print *, '      Therefore no space charge kick will be applied.'
  return
endif

call radiation_integrals(lat, closed_orb, mode, ix_branch = lttp%start%ix_branch)
if (lttp%a_emittance /= 0) mode%a%emittance = lttp%a_emittance
if (lttp%b_emittance /= 0) mode%b%emittance = lttp%b_emittance
n_particle = abs(beam_init%bunch_charge / (e_charge * charge_of(closed_orb(0)%species)))
call setup_ultra_rel_space_charge_calc (.true., lat, n_particle, mode, closed_orb)

end subroutine ltt_setup_space_charge

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_write_particle_data (lttp, i_turn, bunch, rad_map)

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
