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
use s_fitting_new, only: probe, internal_state, track_probe, assignment(=), operator(+), default, spin0
use superimpose_mod, only: add_superimpose
use radiation_mod
use expression_mod

implicit none

integer, parameter :: master_rank$  = 0
integer, parameter :: results_tag$  = 1000
integer, parameter :: is_done_tag$  = 1001
integer, parameter :: particle_tag$ = 1002
integer, parameter :: num_tag$      = 1003

! Essentially: The ltt_params_struct holds user setable parameters while the ltt_com_struct holds
! parameters that are not setable.

type ltt_column_struct
  character(120) :: param = ''
  character(40) :: header_str = ''
  character(20) :: format = ''
end type

! User settable parameters

type ltt_params_struct
  character(20) :: simulation_mode = ''       ! CHECK, SINGLE, BEAM, STAT
  character(20) :: tracking_method = 'BMAD'   ! MAP, PTC, BMAD
  character(100) :: exclude_from_maps = 'beambeam::*'
  character(40) :: ele_start = ''
  character(40) :: ele_stop = ''
  character(200) :: lat_file = ''
  character(200) :: particle_output_file = ''
  character(200) :: beam_binary_output_file = ''
  character(200) :: sigma_matrix_output_file = ''
  character(200) :: custom_output_file = ''
  character(200) :: map_file_prefix = ''
  character(200) :: averages_output_file = ''
  type (ltt_column_struct) column(100)
  integer :: n_turns = -1
  integer :: random_seed = 0
  integer :: map_order = -1
  integer :: averaging_window = 1
  integer :: particle_output_every_n_turns = -1
  integer :: averages_output_every_n_turns = -1
  real(rp) :: ramping_start_time
  real(rp) :: ptc_aperture(2) = 0.1
  real(rp) :: print_on_dead_loss = -1
  real(rp) :: timer_print_dtime = 120
  real(rp) :: dead_cutoff = 1
  real(rp) :: a_emittance = 0   ! Used for space charge calculation.
  real(rp) :: b_emittance = 0   ! Used for space charge calculation.
  logical :: only_live_particles_out = .true.
  logical :: ramping_on = .false.
  logical :: rfcavity_on = .true.
  logical :: add_closed_orbit_to_init_position = .true.
  logical :: symplectic_map_tracking = .false.
  logical :: split_bends_for_radiation = .false.
  integer :: mpi_runs_per_subprocess = 4        ! Number of runs a slave process will take on average.
  logical :: debug = .false.
  logical :: regression_test = .false.          ! Only used for regression testing. Not of general interest.
end type

! A section structure is either:
!   1) A map between two points.
!   2) A pointer to an element where:
!         a) Radiation is to be put in. Or
!         b) The element is to be tracked.

integer, parameter :: map$ = 1, ele$ = 2

type ltt_section_struct
  integer :: type = 0
  type (ptc_map_with_rad_struct), allocatable :: map
  type (ele_struct), pointer :: ele => null()
end type

! Common vars

type ltt_com_struct
  type (lat_struct) :: lat
  type (lat_struct) :: tracking_lat      ! Used for tracking with PTC and maps. Can contain radiation markers for SLICK sectioning.
  type (internal_state) ptc_state
  type (coord_struct), allocatable :: bmad_closed_orb(:)
  type (normal_modes_struct) modes
  type (ltt_section_struct), allocatable :: sec(:)   ! Array of sections indexed from 0. The first one marks the beginning.
  type (ele_pointer_struct), allocatable :: ramper(:)    ! Ramper element locations.
  type (bunch_struct) bunch                    ! Used when regression testing.
  integer :: n_ramper_loc = 0
  integer, allocatable :: ix_wake_ele(:)     ! List of element indexes where a wake is applied.
  integer :: ix_branch = 0                   ! Lattice branch being tracked.
  real(rp) :: ptc_closed_orb(6) = 0
  real(rp) :: time_start = 0
  logical :: wrote_particle_file_header = .false.
  logical :: ramp_in_track1_preprocess = .false.
  logical :: debug = .false.
  integer :: n_particle      ! Num particles per bunch. Needed with MPI.
  integer :: mpi_rank = master_rank$
  integer :: mpi_run_index = 0                 ! Run index
  integer :: mpi_ix0_particle = 0              ! Index of first particle
  logical :: using_mpi = .false.
  character(200) :: mpi_data_dir = ''          ! Temporary directory For storing particle position data when running with mpi.
  character(200) :: master_input_file = ''
end type

integer, parameter :: new$ = 0,  valid$ = 1, written$ = 2

type ltt_bunch_data_struct
  integer :: n_live = 0     ! Number alive on last turn used for averaging.
  integer :: n_count = 0    ! Number of particles counted over all turns used for averaging.
  real(rp) :: orb_sum(6) = 0    ! Orbit average
  real(rp) :: orb2_sum(6,6) = 0
  real(rp) :: spin_sum(3) = 0   ! Spin
  real(rp) :: p0c_sum = 0
  real(rp) :: time_sum = 0
  integer :: species = 0
end type

type ltt_turn_data_struct
  type (ltt_bunch_data_struct), allocatable :: bunch(:)
  integer :: status = new$  ! Has all data been gathered (valid$)? Has data been written (written$)?
  integer :: i_turn = 0
end type

type ltt_beam_data_struct
  type (ltt_turn_data_struct), allocatable :: turn(:)
end type

type (ltt_params_struct), pointer, save :: ltt_params_global   ! Needed for track1_preprocess and track1_bunch_hook
type (ltt_com_struct),    pointer, save :: ltt_com_global      ! Needed for track1_preprocess and track1_bunch_hook

contains

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_init_params(ltt, ltt_com, beam_init)

type (ltt_params_struct), target :: ltt
type (ltt_com_struct), target :: ltt_com
type (beam_init_struct) beam_init
type (lat_struct), pointer :: lat
type (branch_struct), pointer :: branch
type (ele_pointer_struct), allocatable :: eles(:)

real(rp) dummy
integer ir, i, n, ix, n_loc
character(200) arg
character(40) m_name
character(*), parameter :: r_name = 'ltt_init_params'
logical err

namelist / params / bmad_com, beam_init, ltt

!

ltt_params_global => ltt
ltt_com_global => ltt_com

! Parse command line

lat => ltt_com%lat
ltt_com%master_input_file = ''

i = 0
do while (i < cesr_iargc())
  i = i + 1
  call cesr_getarg(i, arg)
  call match_word (arg, ['-debug'], ix, .false., .true., m_name)
  select case (m_name)
  case ('-debug')
    ltt_com%debug = .true.
  case default
    if (ltt_com%master_input_file /= '') then
      print '(2a)', 'Extra stuff on the command line: ', quote(arg)
      print '(a)',  'Stopping here.'
      stop
    endif
    ltt_com%master_input_file = arg
  end select
end do

if (ltt_com%master_input_file == '') ltt_com%master_input_file = 'long_term_tracking.init'

! Read parameters

if (.not. ltt_com%using_mpi .or. ltt_com%mpi_rank == master_rank$) then
  print '(2a)', 'Initialization file: ', trim(ltt_com%master_input_file)
endif

open (1, file = ltt_com%master_input_file, status = 'old', action = 'read')
read (1, nml = params)
close (1)

call upcase_string(ltt%simulation_mode)

! Lattice init

bmad_com%auto_bookkeeper = .false.

call ran_seed_put (ltt%random_seed)
call ptc_ran_seed_put (ltt%random_seed)

if (ltt_com%using_mpi) then
  call ran_seed_get (ir)
  call ran_seed_put (ir + 10 * ltt_com%mpi_rank)
  call ptc_ran_seed_put (ir + 10 * ltt_com%mpi_rank)
endif

call bmad_parser (ltt%lat_file, lat)

! Read the master input file again so that bmad_com parameters set in the file
! take precedence over bmad_com parameters set in the lattice file.

open (1, file = ltt_com%master_input_file, status = 'old', action = 'read')
read (1, nml = params)  
close (1)

! Sanity checks

select case (ltt%simulation_mode)
case ('CHECK', 'SINGLE', 'BEAM', 'STAT')
case ('BUNCH')
  print '(a)', '"BUNCH" SETTING FOR LTT%SIMULATION_MODE HAS BEEN CHANGED TO "BEAM"'
case default
  print '(a)', 'UNKNOWN LTT%SIMULATION_MODE: ' // ltt%simulation_mode
  stop
end select

if (ltt%simulation_mode /= 'CHECK') then
  select case (ltt%tracking_method)
  case ('', 'MAP', 'PTC', 'BMAD')
  case default
    print '(a)', 'UNKNOWN LTT%TRACKING_METHOD: ' // ltt%tracking_method
    stop
  end select
endif

! There is a possible problem with splitting the bends in a lattice if the 
! ele_start or ele_stop names are the index of the element.
! In this case the element index in the tracking_lat will not be the same. 
! This is solved by setting the type string to something unique.

if (ltt%ele_start /= '') then
  call lat_ele_locator (ltt%ele_start, lat, eles, n_loc, err)
  if (err .or. n_loc == 0) then
    print '(2a)', 'Starting element not found: ', trim(ltt%ele_start)
    stop
  endif
  if (n_loc > 1) then
    print '(2a)', 'Multiple elements found with ele_start name: ', trim(ltt%ele_start)
    print '(a)', 'Will stop here.'
    stop
  endif
  branch => pointer_to_branch(eles(1)%ele)
  ltt_com%ix_branch = branch%ix_branch
  eles(1)%ele%type = '@START_ELE'
else
  ltt_com%ix_branch = 0
endif

branch => lat%branch(ltt_com%ix_branch)

if (ltt%simulation_mode == 'CHECK') then
  if (ltt%ele_stop /= '' .and. ltt%ele_stop /= ltt%ele_start) then
    call lat_ele_locator (ltt%ele_stop, lat, eles, n_loc, err)
    if (err .or. n_loc == 0) then
      print '(2a)', 'Stopping element not found: ', trim(ltt%ele_stop)
      stop
    endif
    if (n_loc > 1) then
      print '(2a)', 'Multiple elements found with ele_stop name: ', trim(ltt%ele_stop)
      print '(a)', 'Will stop here.'
      stop
    endif
    eles(1)%ele%type = '@STOP_ELE'
  endif
endif

if (ltt_com%using_mpi .and. ltt%tracking_method == 'PTC') then
  print '(a)', 'TRACKING_METHOD = PTC NOT ALLOWED WITH MPI VERSION OF THIS PROGRAM.'
  stop
endif

if (beam_init%use_particle_start_for_center) beam_init%center = lat%particle_start%vec

!

if (ltt%ramping_on) then
  if (ltt%tracking_method /= 'BMAD') THEN
    print *, 'NOTE: ltt%tracking_method MUST BE SET TO "BMAD" IF LTT%RAMPING_ON IS SET TO TRUE.'
  endif

  call lat_ele_locator ('RAMPER::*', lat, ltt_com%ramper, ltt_com%n_ramper_loc, err)
  if (ltt_com%n_ramper_loc == 0) then
    print '(2a)', 'Warning! NO RAMPER ELEMENTS FOUND IN LATTICE.'
    stop
  endif
endif

! Get list of wake elements.

n = 0
do i = 1, branch%n_ele_track
  if (.not. associated(pointer_to_wake_ele(branch%ele(i), dummy))) cycle
  n = n + 1
  call re_allocate(ltt_com%ix_wake_ele, n)
  ltt_com%ix_wake_ele(n) = i
enddo

!

call fullfilename (ltt%particle_output_file, ltt%particle_output_file)
call fullfilename (ltt%sigma_matrix_output_file, ltt%sigma_matrix_output_file)
call fullfilename (ltt%custom_output_file, ltt%custom_output_file)

end subroutine ltt_init_params

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_init_tracking (lttp, ltt_com)

type (ltt_params_struct) lttp
type (ltt_com_struct), target :: ltt_com
type (lat_struct), pointer :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele, ele_start, ele_stop

real(rp) closed_orb(6)

integer i, ix_branch, ib, n_slice, ie, ir

logical err, map_file_exists

character(40) start_name, stop_name

! PTC has an internal aperture of 1.0 meter. To be safe, use an aperture of 0.9 meter

call ltt_make_tracking_lat(lttp, ltt_com)
lat => ltt_com%tracking_lat

lttp%ptc_aperture = min([0.9_rp, 0.9_rp], lttp%ptc_aperture)

if (lttp%tracking_method == 'PTC' .or. lttp%simulation_mode == 'CHECK') then
  if (.not. associated(lat%branch(0)%ptc%m_t_layout)) call lat_to_ptc_layout(lat)
  call ptc_setup_tracking_with_damping_and_excitation(lat%branch(0), bmad_com%radiation_damping_on, &
                                          bmad_com%radiation_fluctuations_on, ltt_com%ptc_state, ltt_com%ptc_closed_orb)
  if (bmad_com%spin_tracking_on) ltt_com%ptc_state = ltt_com%ptc_state + SPIN0
endif

! If using ramping elements, setup the lattice using lttp%ramping_start_time

if (lttp%ramping_on) then
  do i = 1, ltt_com%n_ramper_loc
    ltt_com%ramper(i)%ele%control%var(1)%value = lttp%ramping_start_time
  enddo
  branch => lat%branch(ltt_com%ix_branch)
  do ie = 0, branch%n_ele_max
    call apply_ramper (branch%ele(ie), ltt_com%ramper(1:ltt_com%n_ramper_loc), err)
  enddo
endif  

!

if (.not. lttp%rfcavity_on) call set_on_off (rfcavity$, lat, off$)
call twiss_and_track (ltt_com%tracking_lat, ltt_com%bmad_closed_orb, ix_branch = ltt_com%ix_branch)
call radiation_integrals (lat, ltt_com%bmad_closed_orb, ltt_com%modes, ix_branch = ltt_com%ix_branch)

if (lttp%simulation_mode == 'CHECK') bmad_com%radiation_fluctuations_on = .false.

if (lttp%simulation_mode == 'CHECK' .or. lttp%tracking_method == 'MAP') then
  call ltt_read_map (lttp, ltt_com, err)
  if (err) then
    if (ltt_com%mpi_rank == master_rank$) then
      call ltt_make_map(lttp, ltt_com)
      call ltt_write_map(lttp, ltt_com)
    else
      do 
        call milli_sleep(1000)
        call ltt_read_map(lttp, ltt_com, err)
        if (.not. err) exit
      enddo
    endif
  endif
endif

if (lttp%simulation_mode == 'CHECK' .or. lttp%tracking_method == 'PTC') call lat_to_ptc_layout(ltt_com%tracking_lat)

end subroutine ltt_init_tracking

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_print_inital_info (lttp, ltt_com)

type (ltt_params_struct) lttp
type (ltt_com_struct), target :: ltt_com
integer i, nm

! Print some info.

nm = 0
if (allocated(ltt_com%sec)) then
  do i = 0, ubound(ltt_com%sec, 1)
    if (allocated(ltt_com%sec(i)%map)) nm = nm + 1
  enddo
endif

print *
print '(a)', '--------------------------------------'
print '(a, a)',    'ltt%lat_file:                       ', quote(lttp%lat_file)
print '(a, a)',    'ltt%sigma_matrix_output_file:       ', quote(lttp%sigma_matrix_output_file)
print '(a, a)',    'ltt%particle_output_file:           ', quote(lttp%particle_output_file)
print '(a, a)',    'ltt%averages_output_file:           ', quote(lttp%averages_output_file)
print '(a, a)',    'ltt%custom_output_file:             ', quote(lttp%custom_output_file)
print '(a, a)',    'ltt%map_file_prefix:                ', quote(lttp%map_file_prefix)
print '(a, a)',    'ltt%simulation_mode:                ', trim(lttp%simulation_mode)
print '(a, a)',    'ltt%tracking_method:                ', trim(lttp%tracking_method)
if (lttp%tracking_method == 'MAP' .or. lttp%simulation_mode == 'CHECK') then
  print '(a, i0)', 'ltt%map_order:                      ', lttp%map_order
  print '(a, a)',  'ltt%ele_start:                      ', quote(lttp%ele_start)
  print '(a, a)',  'ltt%ele_stop:                       ', quote(lttp%ele_stop)
  print '(a, a)',  'ltt%exclude_from_maps:              ', quote(lttp%exclude_from_maps)
  print '(a, l1)', 'ltt%split_bends_for_radiation:      ', lttp%split_bends_for_radiation
  print '(a, l1)', 'ltt%symplectic_map_tracking:        ', lttp%symplectic_map_tracking
  print '(a, i0)', 'Number of maps:                     ', nm
endif
print '(a, l1)',   'bmad_com%radiation_damping_on:      ', bmad_com%radiation_damping_on
print '(a, l1)',   'bmad_com%radiation_fluctuations_on: ', bmad_com%radiation_fluctuations_on
print '(a, l1)',   'bmad_com%spin_tracking_on:          ', bmad_com%spin_tracking_on
print '(a, l1)',   'bmad_com%sr_wakes_on:               ', bmad_com%sr_wakes_on
if (bmad_com%sr_wakes_on) then
  print '(a, i8)',   'Number of wake elements:            ', size(ltt_com%ix_wake_ele)
endif
print '(a, i8)',   'ltt%n_turns:                        ', lttp%n_turns
print '(a, i8)',   'ltt%particle_output_every_n_turns:  ', lttp%particle_output_every_n_turns
print '(a, i8)',   'ltt%averages_output_every_n_turns:  ', lttp%averages_output_every_n_turns
print '(a, l1)',   'ltt%ramping_on:                     ', lttp%ramping_on
print '(2a)',      'ltt%ramping_start_time:             ', real_str(lttp%ramping_start_time, 6)
print '(a, i8)',   'ltt%averaging_window:               ', lttp%averaging_window
print '(a)', '--------------------------------------'
print *

end subroutine ltt_print_inital_info

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_run_check_mode (lttp, ltt_com, beam_init)

type (ltt_params_struct) lttp
type (beam_init_struct) beam_init
type (ltt_com_struct), target :: ltt_com
type (lat_struct), pointer :: lat
type (coord_struct), allocatable :: orb(:)
type (coord_struct) orb_map, orb_map1, orb_map2, orb_bmad, orb_bmad1, orb_bmad2, orb_init
type (ele_struct), pointer :: ele_start, ele_stop
type (probe) prb_ptc, prb_ptc1, prb_ptc2

real(rp) mat_map(6,6), mat_bmad(6,6), mat_ptc(6,6), spin_ptc(3)
integer i
integer track_state, ix_branch, ix_start, ix_stop

! Run serial in check mode.

bmad_com%radiation_fluctuations_on = .false.
lat => ltt_com%tracking_lat

call ltt_pointer_to_map_ends(lttp, lat, ele_start, ele_stop)
ix_branch = ltt_com%ix_branch
ix_start = ele_start%ix_ele
ix_stop  = ele_stop%ix_ele

call init_coord (orb_init, beam_init%center, ele_start, downstream_end$, spin = beam_init%spin)
call reallocate_coord(orb, lat)

if (lttp%add_closed_orbit_to_init_position) orb_init%vec = orb_init%vec + ltt_com%bmad_closed_orb(ix_start)%vec

call track_check_all (orb_init, orb, 0, 0, orb_map, orb_bmad, prb_ptc, spin_ptc, ele_start, ele_stop)

if (allocated(ltt_com%sec(1)%map)) print '(a, 6f14.8)', 'Map closed orbit at start:  ', ltt_com%sec(1)%map%sub_map(1)%fix0
print '(a, 6f14.8)', 'Bmad closed orbit at start: ', ltt_com%bmad_closed_orb(ix_start)%vec

print '(a)'
print '(a, 6f14.8)', 'Starting orbit for tracking:', orb_init%vec
print '(a, 6f14.8)', 'Initial spin: ', orb_init%spin

print '(a)'
print '(a)', 'Phase Space at Track End (Bmad closed orbit used with map tracking:'
print '(a, 6f14.8)', 'Bmad tracking:    ', orb_bmad%vec
print '(a, 6f14.8)', 'PTC tracking:     ', prb_ptc%x
print '(a, 6f14.8)', 'Map tracking:     ', orb_map%vec
print *
print '(a, 6f14.8)', 'Max Diff Tracking:', max(orb_map%vec, prb_ptc%x, orb_bmad%vec) - &
                                           min(orb_map%vec, prb_ptc%x, orb_bmad%vec)

print '(a)'
print '(a)', 'Spin at Track End:'
print '(a, 6f14.8)', 'Bmad tracking:   ', orb_bmad%spin
print '(a, 6f14.8)', 'PTC tracking:    ', spin_ptc
print '(a, 6f14.8)', 'Map tracking:    ', orb_map%spin
print *
print '(a, 6f14.8)', 'Max Diff Spin:   ', max(orb_map%spin, spin_ptc, orb_bmad%spin) - &
                                          min(orb_map%spin, spin_ptc, orb_bmad%spin) 
!

do i = 1, 6
  call track_check_all (orb_init, orb, i, -1, orb_map1, orb_bmad1, prb_ptc1, spin_ptc, ele_start, ele_stop)
  call track_check_all (orb_init, orb, i, +1, orb_map2, orb_bmad2, prb_ptc2, spin_ptc, ele_start, ele_stop)
  mat_ptc(1:6, i)   = (prb_ptc2%x - prb_ptc1%x) / (2 * bmad_com%d_orb(i))
  mat_bmad(1:6, i)  = (orb_bmad2%vec - orb_bmad1%vec) / (2 * bmad_com%d_orb(i))
  mat_map(1:6, i)   = (orb_map2%vec - orb_map1%vec) / (2 * bmad_com%d_orb(i))
enddo

print '(a)'
print '(a)', 'Transfer Matrix for Each Tracking Method Formed Using Finite Differences:'

do i = 1, 6
  print '(a, i0, 4x, 6es14.6)', ' Bmad:MatrixRow:', i, mat_bmad(i,:)
  print '(a, i0, 4x, 6es14.6)', '  PTC:MatrixRow:', i, mat_ptc(i,:)
  print '(a, i0, 4x, 6es14.6)', '  Map:MatrixRow:', i, mat_map(i,:)
  print '(a)'
enddo

!------------------------------------------------------------
contains

subroutine track_check_all (orb_init, orb, i_ps, sgn, orb_map, orb_bmad, prb_ptc, spin_ptc, ele_start, ele_stop)

type (coord_struct) orb(0:), orb_init, orb_map, orb_bmad, orb_start
type (probe) prb_ptc
type (lat_struct), pointer :: lat
type (ele_struct), pointer :: ele_start, ele_stop

real(rp) spin_ptc(3)
integer i_ps, sgn
logical is_lost

!

lat => ltt_com%tracking_lat
orb_start = orb_init
if (i_ps > 0) orb_start%vec(i_ps) = orb_start%vec(i_ps) + sgn * bmad_com%d_orb(i)

! Map

orb_map = orb_start
call ltt_track_map (lttp, ltt_com, orb_map, is_lost)

! Bmad

orb(ix_start) = orb_start
call track_many (lat, orb, ix_start, ix_stop, 1, ix_branch, track_state)
orb_bmad = orb(ix_stop)

! PTC

prb_ptc = orb_start%vec
prb_ptc%q%x = [1, 0, 0, 0]   ! Unit quaternion
if (ix_stop == ix_start) then
  call track_probe (prb_ptc, ltt_com%ptc_state, fibre1 = pointer_to_fibre(ele_start))
else
  call track_probe (prb_ptc, ltt_com%ptc_state, fibre1 = pointer_to_fibre(ele_start), &
                                                fibre2 = pointer_to_fibre(ele_stop))
endif
spin_ptc = rotate_vec_given_quat(prb_ptc%q%x, orb_start%spin)

end subroutine track_check_all

end subroutine ltt_run_check_mode

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_run_single_mode (lttp, ltt_com, beam_init)

use radiation_mod, only: track1_radiation_center

type (ltt_params_struct) lttp
type (ltt_com_struct), target :: ltt_com
type (beam_init_struct) beam_init
type (branch_struct), pointer :: branch
type (lat_struct), pointer :: lat
type (coord_struct), allocatable :: orb(:)
type (coord_struct) :: orbit, orbit_old
type (ele_struct), pointer :: ele
type (ele_struct), pointer :: ele_start
type (probe) prb

real(rp) average(6), sigma(6,6), dt
integer i, n_sum, iu_part, i_turn, track_state, ix_branch
logical is_lost

character(40) fmt

! Run serial in single mode.

lat => ltt_com%tracking_lat

n_sum = 0
sigma = 0
average = 0
ix_branch = ltt_com%ix_branch
branch => lat%branch(ix_branch)
call ltt_pointer_to_map_ends(lttp, lat, ele_start)

orbit_old%vec = real_garbage$
ltt_com%ramp_in_track1_preprocess = .true.

call ltt_setup_high_energy_space_charge(lttp, ltt_com, branch, beam_init)

if (lttp%tracking_method == 'BMAD') call reallocate_coord (orb, lat)
call init_coord (orbit, beam_init%center, ele_start, downstream_end$, lat%param%particle, spin = beam_init%spin)

if (lttp%add_closed_orbit_to_init_position) then
  select case (lttp%tracking_method)
  case ('MAP');     orbit%vec = orbit%vec + ltt_com%bmad_closed_orb(ele_start%ix_ele)%vec
  case ('BMAD');    orbit%vec = orbit%vec + ltt_com%bmad_closed_orb(ele_start%ix_ele)%vec
  case ('PTC');     orbit%vec = orbit%vec + ltt_com%ptc_closed_orb
  end select
endif

!

fmt = '(i6, 6es16.8, 3x, 3f10.6)'
iu_part = lunget()
if (lttp%particle_output_file == '') lttp%particle_output_file = 'single.dat'
open(iu_part, file = lttp%particle_output_file, recl = 200)
call ltt_write_params_header(lttp, ltt_com, iu_part, 1)
write (iu_part, '(a)') '# Turn |            x              px               y              py               z              pz    |   spin_x    spin_y    spin_z'
write (iu_part, fmt) 0, orbit%vec, orbit%spin

if (lttp%custom_output_file /= '') call ltt_write_custom (0, lttp, ltt_com, orbit = orbit)

do i_turn = 1, lttp%n_turns
  select case (lttp%tracking_method)

  case ('BMAD')
    orb(ele_start%ix_ele) = orbit
    call track_many (lat, orb, ele_start%ix_ele, ele_start%ix_ele, 1, ix_branch, track_state)
    is_lost = (track_state /= moving_forward$)
    orbit = orb(ele_start%ix_ele)

  case ('PTC')
    prb = orbit%vec
    prb%q%x = [1, 0, 0, 0]  ! Unit quaternion
    call track_probe (prb, ltt_com%ptc_state, fibre1 = pointer_to_fibre(ele_start))
    orbit%vec = prb%x
    orbit%spin = rotate_vec_given_quat(prb%q%x, orbit%spin)
    is_lost = (abs(orbit%vec(1)) > lttp%ptc_aperture(1) .or. abs(orbit%vec(3)) > lttp%ptc_aperture(2))
    if (all(orbit%vec == orbit_old%vec)) is_lost = .true.
    orbit_old%vec = orbit%vec

  case ('MAP')
    call ltt_track_map(lttp, ltt_com, orbit, is_lost)

  case default
    print '(a)', 'Unknown tracking_method: ' // lttp%tracking_method
    stop
  end select

  if (lttp%particle_output_every_n_turns > 0 .and. modulo(i_turn, lttp%particle_output_every_n_turns) == 0) then
    write (iu_part, fmt) i_turn, orbit%vec, orbit%spin
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

  if (lttp%custom_output_file /= '') call ltt_write_custom (i_turn, lttp, ltt_com, orbit = orbit)
enddo

print '(2a)', 'Particle output file: ', trim(lttp%particle_output_file)
close (iu_part)

if (lttp%sigma_matrix_output_file /= '') then
  call ltt_write_single_mode_sigma_file (lttp, n_sum, average, sigma)
endif

end subroutine ltt_run_single_mode

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_run_beam_mode (lttp, ltt_com, beam_init, beam_data_in, beam_in)

type (ltt_params_struct) lttp
type (beam_init_struct) beam_init
type (beam_struct), target :: beam0
type (beam_struct), pointer :: beam
type (beam_struct), optional, target :: beam_in  ! Used with mpi
type (bunch_struct), pointer :: bunch
type (ltt_com_struct), target :: ltt_com
type (lat_struct), pointer :: lat
type (coord_struct), allocatable :: orb(:)
type (coord_struct), pointer :: p
type (ltt_beam_data_struct), optional :: beam_data_in
type (ltt_beam_data_struct), target :: beam_data
type (ele_struct), pointer :: ele_start
type (probe) prb

real(rp) time0, time1, time_now

integer n, n_print_dead_loss, i_turn, ie, n_part_tot
integer ip, ib, n_live_old, n_live, ix_branch

logical err_flag, is_lost

character(16) prefix_str
character(200) file_name

! Init

ltt_com%mpi_run_index = ltt_com%mpi_run_index + 1
lat => ltt_com%tracking_lat
call run_timer('ABS', time0)
time1 = time0

ix_branch = ltt_com%ix_branch
call ltt_pointer_to_map_ends(lttp, lat, ele_start)
  
if (ltt_com%using_mpi) then
  beam => beam_in
  prefix_str = 'Slave ' // int_str(ltt_com%mpi_rank) // ':'
else
  call init_beam_distribution (ele_start, lat%param, beam_init, beam0, err_flag, modes = ltt_com%modes)
  if (err_flag) stop
  beam => beam0
  print '(a, i8)',   'n_particle:                    ', size(beam%bunch(1)%particle)
  prefix_str = ''
  ltt_com%n_particle = size(beam%bunch(1)%particle)
endif

if (bmad_com%spin_tracking_on .and. all(beam_init%spin == 0)) then
  ie = ele_start%ix_ele
  do ib = 1, size(beam%bunch)
    forall (n = 1:3) beam%bunch(ib)%particle%spin(n) = ltt_com%bmad_closed_orb(ie)%spin(n)
  enddo
endif

call ltt_setup_high_energy_space_charge(lttp, ltt_com, lat%branch(ix_branch), beam_init)

do ib = 1, size(beam%bunch)
  do n = 1, size(beam%bunch(ib)%particle)
    p => beam%bunch(ib)%particle(n)
    if (lttp%add_closed_orbit_to_init_position) then
      select case (lttp%tracking_method)
      case ('MAP');    p%vec = p%vec + ltt_com%bmad_closed_orb(ele_start%ix_ele)%vec
      case ('BMAD');   p%vec = p%vec + ltt_com%bmad_closed_orb(ele_start%ix_ele)%vec
      case ('PTC');    p%vec = p%vec + ltt_com%ptc_closed_orb
      end select
    endif
  enddo
enddo

n_part_tot = 0
do ib = 1, size(beam%bunch)
  n_part_tot = n_part_tot + size(beam%bunch(ib)%particle)
enddo
n_print_dead_loss = max(1, nint(lttp%print_on_dead_loss * n_part_tot))
n_live_old = n_part_tot

! 

call ltt_write_particle_data (lttp, ltt_com, 0, beam)

! If using mpi then this data will all be written out at the end.
! Here partial writes are used so the user can monitor progress if they want.
call ltt_calc_beam_sums (lttp, 0, beam, beam_data)
if (.not. ltt_com%using_mpi) then
  call ltt_write_beam_averages(lttp, beam_data)
  call ltt_write_sigma_matrix (lttp, beam_data)
  where (beam_data%turn%status == valid$) beam_data%turn%status = written$
endif

if (lttp%custom_output_file /= '' .and. (ltt_com%mpi_rank == master_rank$ .or. ltt_com%mpi_rank == 1) .and. &
                          ltt_com%mpi_run_index == 1) call ltt_write_custom (0, lttp, ltt_com, beam = beam)

do i_turn = 1, lttp%n_turns
  n_live = 0
  do ib = 1, size(beam%bunch)
    bunch => beam%bunch(ib)

    select case (lttp%tracking_method)
    case ('MAP')
      do ip = 1, size(bunch%particle)
        p => bunch%particle(ip)
        if (p%state /= alive$) cycle
        call ltt_track_map (lttp, ltt_com, p, is_lost)
        if (is_lost) p%state = lost$
      enddo

    case ('BMAD')
      call track_bunch (lat, bunch, ele_start, ele_start, err_flag)

    case ('PTC')
      do ip = 1, size(bunch%particle)
        p => bunch%particle(ip)
        if (p%state /= alive$) cycle
        prb = p%vec
        prb%q%x = [1, 0, 0, 0]  ! Unit quaternion
        call track_probe (prb, ltt_com%ptc_state, fibre1 = lat%branch(ix_branch)%ele(1)%ptc_fibre)
        p%vec = prb%x
        p%spin = rotate_vec_given_quat(prb%q%x, p%spin)
        if (abs(p%vec(1)) > lttp%ptc_aperture(1) .or. abs(p%vec(3)) > lttp%ptc_aperture(2)) p%state = lost$
      enddo

    case default
      print '(a)', 'Unknown tracking_method: ' // lttp%tracking_method
      stop
    end select

    n_live = n_live + count(bunch%particle%state == alive$)
  enddo

  if (n_live_old - n_live >= n_print_dead_loss) then
    print '(2a, i0, a, i0)', trim(prefix_str), ' Cumulative number dead on turn ', i_turn, ': ', n_part_tot - n_live
    n_live_old = n_live
  endif

  if (n_live  == 0 .and. .not. ltt_com%using_mpi) then
    print '(2a)', trim(prefix_str), ' NO LIVE PARTICLES! STOPPING NOW.'
    exit
  endif

  if (n_part_tot - n_live > lttp%dead_cutoff * n_part_tot .and. .not. ltt_com%using_mpi) then
    print '(2a)', trim(prefix_str), ' Particle loss greater than set by dead_cutoff. Stopping now.'
    exit
  endif

  call run_timer('ABS', time_now)
  if (time_now-time1 > lttp%timer_print_dtime .and. .not. ltt_com%using_mpi) then
    print '(2a, f10.2, a, i0)', trim(prefix_str), ' Ellapsed time (min): ', (time_now-time0)/60, ', At turn: ', i_turn
    time1 = time_now
  endif

  call ltt_write_particle_data (lttp, ltt_com, i_turn, beam)

  ! If using mpi then this data will all be written out at the end.
  ! Here partial writes are used so the user can monitor progress if they want.

  call ltt_calc_beam_sums (lttp, i_turn, beam, beam_data)

  if (.not. ltt_com%using_mpi) then
    call ltt_write_beam_averages(lttp, beam_data)
    call ltt_write_sigma_matrix (lttp, beam_data)
    where (beam_data%turn%status == valid$) beam_data%turn%status = written$
  endif

  if (lttp%custom_output_file /= '' .and. (ltt_com%mpi_rank == master_rank$ .or. ltt_com%mpi_rank == 1) .and. &
                          ltt_com%mpi_run_index == 1) call ltt_write_custom (i_turn, lttp, ltt_com, beam = beam)
end do

if (lttp%beam_binary_output_file /= '') then
  if (ltt_com%using_mpi) then
    file_name = trim(ltt_com%mpi_data_dir) // 'b' // int_str(ltt_com%mpi_ix0_particle)
  else
    file_name = lttp%beam_binary_output_file
  endif
  call hdf5_write_beam (file_name, beam%bunch, .false., err_flag)
endif

if (present(beam_data_in)) beam_data_in = beam_data
if (lttp%regression_test) ltt_com%bunch = beam%bunch(1)

end subroutine ltt_run_beam_mode

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_run_stat_mode (lttp, ltt_com)

type (ltt_params_struct) lttp
type (ltt_com_struct), target :: ltt_com
type (lat_struct), pointer :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele

real(rp) chrom_x, chrom_y, ring_length
integer i
logical err_flag

! Run serial in stat mode.

lat => ltt_com%tracking_lat
branch => lat%branch(ltt_com%ix_branch)

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
  write(22,'(a,1x,f10.3,4(1x,f10.6))') ele%name, ele%s, ltt_com%bmad_closed_orb(i)%vec(1:5:2), ltt_com%bmad_closed_orb(i)%vec(6)
enddo

close(20)
close(21)
close(22)

!

ring_length = branch%param%total_length
call chrom_calc(lat, 1.0d-6, chrom_x, chrom_y, err_flag, ix_branch = branch%ix_branch)
call calc_z_tune (lat)

print *, 'Momentum Compaction:', ltt_com%modes%synch_int(1)/ring_length
print *, 'dE/E=', ltt_com%modes%sigE_E
print *, 'sig_z(m)=', ltt_com%modes%sig_z
print *, 'emit_I  (m): ',  ltt_com%modes%a%emittance
print *, 'emit_II (m): ',  ltt_com%modes%b%emittance
print *, 'emit_III(m): ',  ltt_com%modes%z%emittance
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

subroutine ltt_setup_high_energy_space_charge(lttp, ltt_com, branch, beam_init)

type (ltt_params_struct) lttp
type (branch_struct) branch
type (beam_init_struct) beam_init
type (ltt_com_struct), target :: ltt_com
type (normal_modes_struct) modes
real(rp) n_particle

!

if (.not. branch%param%high_energy_space_charge_on) return

if (lttp%tracking_method == 'MAP') then
  print '(a)', 'NOTE: Space effects are not present when using a map tracking!'
  return
endif

if (.not. lttp%rfcavity_on) then
  print '(a)', 'NOTE: RF is not on. Cannot calculate a longitudinal bunch length.'
  print '(a)', '      Therefore no space charge kick will be applied.'
  return
endif

modes = ltt_com%modes
if (beam_init%a_emit > 0) modes%a%emittance = beam_init%a_emit
if (beam_init%b_emit > 0) modes%a%emittance = beam_init%b_emit
n_particle = abs(beam_init%bunch_charge / (e_charge * charge_of(ltt_com%bmad_closed_orb(0)%species)))
call setup_high_energy_space_charge_calc (.true., branch, n_particle, modes)

end subroutine ltt_setup_high_energy_space_charge

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_write_particle_data (lttp, ltt_com, i_turn, beam, use_mpi)

type (ltt_params_struct) lttp
type (ltt_com_struct), target :: ltt_com
type (beam_struct), target :: beam
type (bunch_struct), pointer :: bunch
type (coord_struct), pointer :: p

integer i_turn, ix, ip, ib, j
integer iu_part

character(200) file_name
character(40) fmt, str

logical, optional :: use_mpi
logical using_mpi
!

if (lttp%particle_output_file == '') return
if (lttp%particle_output_every_n_turns == -1 .and. i_turn /= lttp%n_turns) return
if (lttp%particle_output_every_n_turns == 0 .and. i_turn /= 0 .and. i_turn /= lttp%n_turns) return
if (lttp%particle_output_every_n_turns > 0 .and. modulo(i_turn, lttp%particle_output_every_n_turns) /= 0) return

!

using_mpi = logic_option(ltt_com%using_mpi, use_mpi)

do ib = 1, size(beam%bunch)
  bunch => beam%bunch(ib)

  if (using_mpi) then
    file_name = trim(ltt_com%mpi_data_dir) // int_str(i_turn) // '_' // int_str(ib) // '_' // &
                                                       int_str(ltt_com%mpi_rank) // '_' // int_str(ltt_com%mpi_run_index)
  else
    j = int(log10(real(lttp%n_turns, rp)) + 1 + 1d-10)
    str = int_str(i_turn, j)
    if (size(beam%bunch) > 1) str = int_str(ib) // '-' // str

    ix = index(lttp%particle_output_file, '#')
    if (ix == 0) then
      file_name = trim(lttp%particle_output_file) // str
    else
      file_name = lttp%particle_output_file(1:ix-1) // trim(str) // trim(lttp%particle_output_file(ix+1:))
    endif
  endif

  !

  iu_part = lunget()

  if (using_mpi) then
    open (iu_part, file = file_name, recl = 300)
  else
    open (iu_part, file = file_name, recl = 300)
    call ltt_write_params_header(lttp, ltt_com, iu_part, ltt_com%n_particle, size(beam%bunch))
    write (iu_part, '(a)')  '#      Ix     Turn |           x              px               y              py               z              pz   |     spin_x    spin_y    spin_z    State'
  endif

  ltt_com%wrote_particle_file_header = .true.

  !

  do ip = 1, size(bunch%particle)
    p => bunch%particle(ip)
    ix = ip + ltt_com%mpi_ix0_particle
    if (using_mpi) then
      write (iu_part, '(i9, i9, 6es16.8, 3x, 3f10.6, i5)')  ix, i_turn, p%vec, p%spin, p%state
    else
      if (lttp%only_live_particles_out .and. p%state /= alive$) cycle
      write (iu_part, '(i9, i9, 6es16.8, 3x, 3f10.6, 4x, a)')  ix, i_turn, p%vec, p%spin, trim(coord_state_name(p%state))
    endif
  enddo

  close(iu_part)
enddo

end subroutine ltt_write_particle_data

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
! The optional argument "beam" is passed when doing beam tracking. Otherwise orbit is passed. Never both.


subroutine ltt_write_custom (i_turn, lttp, ltt_com, orbit, beam)

type (ltt_params_struct), target :: lttp
type (ltt_com_struct), target :: ltt_com
type (coord_struct), optional :: orbit
type (beam_struct), optional, target :: beam
type (bunch_struct), pointer :: bunch
type (coord_struct) orb
type (coord_struct), allocatable :: orb_b(:)
type (ltt_column_struct), pointer :: col
type (ele_pointer_struct), allocatable :: eles(:)
type (all_pointer_struct) a_ptr
type (expression_atom_struct), allocatable, target :: stack(:)
type (expression_atom_struct), pointer :: st

integer i_turn
integer i, n, ib, iu, ix, ix1, ix2, is, multi, power, width, digits, n_loc, n_stack

logical err

character(1000) line
character(100) err_str
character(40) fmt, str, ele_name, param
character(4) code

!

if (lttp%custom_output_file == '') return

select case (lttp%averages_output_every_n_turns)
! Stats only for end.
case (-1)
  if (i_turn /= lttp%n_turns) return

! Stats for beginning and end
case (0)
  if (i_turn /= 0 .and. i_turn /= lttp%n_turns) return

! Stats for every %averages_output_every_n_turns 
case default
  if (mod(i_turn, lttp%averages_output_every_n_turns) /= 0) return
end select

!

iu = lunget()

if (i_turn == 0 .or. lttp%averages_output_every_n_turns == -1) then
  open (iu, file = lttp%custom_output_file, recl = 1000)
  line = '#'
  do i = 1, size(lttp%column)
    col => lttp%column(i)
    if (col%param == '') exit
    call parse_fortran_format (col%format, multi, power, code, width, digits)
    line = trim(line) // adjustr(col%header_str(1:width))
  enddo
  write (iu, '(a)') trim(line)

else
  open (iu, file = lttp%custom_output_file, access = 'append', recl = 1000)
endif

!

if (present(orbit)) then
  allocate (orb_b(1))
  orb_b(1) = orbit

else  ! beam is passed
  allocate (orb_b(size(beam%bunch)))

  do ib = 1, size(beam%bunch)
    bunch => beam%bunch(ib)
    n = count(bunch%particle%state == alive$)
    orb_b(ib)%t = sum(bunch%particle%t, bunch%particle%state == alive$) / n
    do i = 1, 6
      orb_b(ib)%vec(i) = sum(bunch%particle%vec(i), bunch%particle%state == alive$) / n
    enddo
    do i = 1, 3
      orb_b(ib)%spin(i) = sum(bunch%particle%spin(i), bunch%particle%state == alive$) / n
    enddo
  enddo
endif

line = ''

do i = 1, size(lttp%column)
  col => lttp%column(i)
  if (col%param == '') exit

  fmt = '(' // trim(col%format) // ')'

  call expression_string_to_stack (col%param, stack, n_stack, err, err_str)
  if (err) then
    print *, err_str
    exit
  endif

  do is = 1, n_stack
    st => stack(is)
    if (st%type /= variable$) cycle

    ix1 = index(st%name, '[')
    if (ix1 /= 0) then
      ix2 = index(st%name, ']')
      if (ix2 == 0) then
        print '(a)', 'NO "]" FOUND IN CUSTOM COLUMN PARAMETER NAME: ' // quote(st%name)
        exit
      endif
      if (st%name(ix2+1:) /= '') then
        print '(a)', 'MALFORMED CUSTOM COLUMN PARAMETER NAME: ' // quote(st%name)
        exit
      endif
      ele_name = st%name(1:ix1-1)
      param = st%name(ix1+1:ix2-1)
      call lat_ele_locator(ele_name, ltt_com%tracking_lat, eles, n_loc, err, .true.)
      if (n_loc == 0) then
        print '(a)', 'CANNOT FIND ELEMENT REFERENCED IN CUSTOM COLUMN PARAMETER: ' // quote(st%name)
        exit
      endif
      call pointer_to_attribute (eles(1)%ele, param, .true., a_ptr, err, .true.)
      if (err) exit
      if (associated(a_ptr%r)) then
        st%value = a_ptr%r
      elseif (associated(a_ptr%i)) then
        st%value = a_ptr%i
      endif
      cycle
    endif

    ix = index(st%name, '@')
    if (ix == 0) then
      orb = orb_b(1)
      str = st%name
    elseif (is_integer(st%name(1:ix-1), ib)) then
      orb = orb_b(ib)
      str = st%name(ix+1:)
    endif

    select case (str)
    case ('n_turn');      st%value = i_turn
    case ('x');           st%value = orb%vec(1)
    case ('px');          st%value = orb%vec(2)
    case ('y');           st%value = orb%vec(3)
    case ('py');          st%value = orb%vec(4)
    case ('z');           st%value = orb%vec(5)
    case ('pz');          st%value = orb%vec(6)
    case ('time');        st%value = orb%t
    case ('sx');          st%value = orb%spin(1)
    case ('sy');          st%value = orb%spin(2)
    case ('sz');          st%value = orb%spin(3)
    case default
    end select
  enddo

  if (col%format(1:1) == 'i') then
    write (str, fmt) nint(expression_stack_value(stack, err, err_str))
  else
    write (str, fmt) expression_stack_value(stack, err, err_str)
  endif

  line = trim(line) // str 
enddo

write (iu, '(2a)') ' ', trim(line)

close (iu)

end subroutine ltt_write_custom

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_write_params_header(lttp, ltt_com, iu, n_particle, n_bunch)

type (ltt_params_struct) lttp
type (ltt_com_struct), target :: ltt_com

integer, optional :: n_bunch
integer iu, n_particle

!

write (iu,  '(3a)')      '# lattice                        = ', quote(lttp%lat_file)
write (iu,  '(3a)')      '# simulation_mode                = ', quote(lttp%simulation_mode)
if (present(n_bunch)) then
  write (iu,  '(a, i8)')   '# n_bunch                        = ', n_bunch
endif
write (iu,  '(a, i8)')   '# n_particle                     = ', n_particle
write (iu,  '(a, i8)')   '# n_turns                        = ', lttp%n_turns
write (iu,  '(a, l1)')   '# ramping_on                     = ', lttp%ramping_on
write (iu,  '(2a)')      '# ramping_start_time             = ', real_str(lttp%ramping_start_time, 6)
write (iu,  '(a, i8)')   '# particle_output_every_n_turns  = ', lttp%particle_output_every_n_turns
write (iu,  '(a, i8)')   '# averages_output_every_n_turns  = ', lttp%averages_output_every_n_turns
write (iu,  '(a, i8)')   '# averaging_window               = ', lttp%averaging_window
write (iu,  '(a, l1)')   '# Radiation_Damping_on           = ', bmad_com%radiation_damping_on
write (iu,  '(a, l1)')   '# Radiation_Fluctuations_on      = ', bmad_com%radiation_fluctuations_on
write (iu,  '(a, l1)')   '# Spin_tracking_on               = ', bmad_com%spin_tracking_on
write (iu,  '(a, l1)')   '# sr_wakes_on                    = ', bmad_com%sr_wakes_on
if (bmad_com%sr_wakes_on) then
  write (iu, '(a, i0)')  '# Number_of_wake_elements        = ', size(ltt_com%ix_wake_ele)
endif
write (iu, '(3a)')       '# Map_file_prefix                = ', quote(lttp%map_file_prefix)
if (lttp%tracking_method == 'MAP') then
  write (iu, '(a, i0)')  '# map_order                      = ', ltt_com%sec(1)%map%map_order
  write (iu, '(a, a)')   '# exclude_from_maps              = ', lttp%exclude_from_maps
  write (iu, '(a, l1)')  '# split_bends_for_radiation      = ', lttp%split_bends_for_radiation
  write (iu, '(a, l1)')  '# symplectic_map_tracking        = ', lttp%symplectic_map_tracking
endif
write (iu, '(a)') '#'

end subroutine ltt_write_params_header

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_allocate_beam_data (lttp, beam_data, n_bunch)

type (ltt_params_struct) lttp
type (ltt_beam_data_struct), target :: beam_data
integer n_bunch, i

!

select case (lttp%averages_output_every_n_turns)
case (-1);      allocate (beam_data%turn(1))    ! Stats only for end.
case (0);       allocate (beam_data%turn(0:1))  ! Stats for beginning and end
case default;   allocate (beam_data%turn(0:1+lttp%n_turns/lttp%averages_output_every_n_turns)) ! Stats for every %averages_output_every_n_turns 
end select

do i = lbound(beam_data%turn, 1), ubound(beam_data%turn, 1)
  allocate (beam_data%turn(i)%bunch(n_bunch))
enddo

end subroutine ltt_allocate_beam_data

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_calc_beam_sums (lttp, i_turn, beam, beam_data)

type (ltt_params_struct) lttp
type (beam_struct), target :: beam
type (bunch_struct), pointer :: bunch
type (ltt_beam_data_struct), target :: beam_data
type (ltt_bunch_data_struct), pointer :: bd

integer i, j, ib, ix, i_turn, this_turn, n2w

!

if (lttp%averages_output_file == '') return
if (.not. allocated(beam_data%turn)) call ltt_allocate_beam_data (lttp, beam_data, size(beam%bunch))

select case (lttp%averages_output_every_n_turns)
! Stats only for end.
case (-1)
  ix = in_this_window(i_turn, lttp%n_turns, lttp%averaging_window, beam_data)
  if (ix < 1) return   ! Out of window. Do not do any averaging

! Stats for beginning and end
case (0)
  ix = in_this_window(i_turn, lttp%n_turns, lttp%averaging_window, beam_data)
  if (ix < 0) return   ! Out of window. Do not do any averaging

! Stats for every %averages_output_every_n_turns 
case default
  ix = in_this_window(i_turn, lttp%averages_output_every_n_turns, lttp%averaging_window, beam_data)
  if (ix < 0) return   ! Out of window. Do not do any averaging
end select

!

do ib = 1, size(beam%bunch)
  bunch => beam%bunch(ib)
  bd => beam_data%turn(ix)%bunch(ib)

  bd%n_live = count(bunch%particle%state == alive$)
  bd%n_count = bd%n_count + bd%n_live
  bd%species = bunch%particle(1)%species

  do i = 1, 6
    bd%orb_sum(i) = bd%orb_sum(i) + sum(bunch%particle%vec(i), bunch%particle%state == alive$) 
    do j = i, 6
      bd%orb2_sum(i,j) = bd%orb2_sum(i,j) + sum(bunch%particle%vec(i) * bunch%particle%vec(j), bunch%particle%state == alive$) 
    enddo
  enddo

  do i = 1, 3
    bd%spin_sum(i) = bd%spin_sum(i) + sum(bunch%particle%spin(i), bunch%particle%state == alive$)
  enddo

  bd%p0c_sum = bd%p0c_sum + sum(bunch%particle%p0c, bunch%particle%state == alive$)
  bd%time_sum = bd%time_sum + sum(bunch%particle%t, bunch%particle%state == alive$)
enddo

!-------------------------------------------------------------------------------------------
contains

function in_this_window(i_turn, n_every, window_width, beam_data) result (ix_this)

type (ltt_beam_data_struct), target :: beam_data
integer i_turn, n_every, window_width, n2w, n, ix_this

!

ix_this = nint(i_turn * 1.0_rp / n_every)
n2w = window_width/2

do n = lbound(beam_data%turn,1), floor(real(i_turn-n2w)/n_every)
  if (beam_data%turn(n)%status == new$) beam_data%turn(n)%status = valid$
enddo

if (i_turn < ix_this*n_every - n2w .or. i_turn > ix_this*n_every + n2w) then
  ix_this = -1  ! Out of window
  return
endif

if (ix_this >= lbound(beam_data%turn,1) .and. ix_this <= ubound(beam_data%turn,1)) &
                                        beam_data%turn(ix_this)%i_turn = ix_this * n_every

end function in_this_window

end subroutine ltt_calc_beam_sums

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_write_beam_averages (lttp, beam_data)

type (ltt_params_struct) lttp
type (ltt_beam_data_struct), target :: beam_data
type (ltt_bunch_data_struct), pointer :: bd
type (bunch_params_struct) bunch_params

real(rp) sig1(6), sigma(6,6)
integer i, j, n, iu, ix, ib
logical error
character(200) :: file_name

!

if (lttp%averages_output_file == '') return

iu = lunget()

n = size(beam_data%turn(1)%bunch)
do ib = 1, n
  ix = index(lttp%averages_output_file, '#')

  if (n == 1) then
    file_name = lttp%averages_output_file
  elseif (ix == 0) then
    file_name = trim(lttp%averages_output_file) // int_str(ib)
  else
    file_name = lttp%averages_output_file(1:ix-1) // int_str(ib) // lttp%averages_output_file(ix+1:)
  endif

  if (beam_data%turn(1)%status == valid$) then
    open (iu, file = file_name, recl = 400)
    write (iu, '(a1, a8, a9, 2a14, 2x, 3a14, 2x, 13a14, 2x, 3a14)') '#', 'Turn', 'N_live', 'Time', 'Polarization', &
                     '<Sx>', '<Sy>', '<Sz>', 'Sig_x', 'Sig_px', 'Sig_y', 'Sig_py', 'Sig_z', 'Sig_pz', &
                     '<x>', '<px>', '<y>', '<py>', '<z>', '<pz>', '<p0c>', 'emit_a', 'emit_b', 'emit_c'
                       
  else
    open (iu, file = lttp%averages_output_file, recl = 400, access = 'append')
  endif

  !

  do ix = 1, size(beam_data%turn)
    if (beam_data%turn(ix)%status /= valid$) cycle
    bd => beam_data%turn(ix)%bunch(ib)
    if (bd%n_count == 0) exit

    do i = 1, 6
    do j = i, 6
      sigma(i,j) = bd%orb2_sum(i,j) / bd%n_count - bd%orb_sum(i) * bd%orb_sum(j) / bd%n_count**2
      sigma(j,i) = sigma(i,j)
    enddo
    enddo

    forall (i = 1:6) sig1(i) = sqrt(max(0.0_rp, sigma(i,i)))

    call calc_emittances_and_twiss_from_sigma_matrix (sigma, 0.0_rp, bunch_params, error)

    write (iu, '(i9, i9, 2f14.9, 2x, 3f14.9, 2x, 13es14.6, 2x, 3es14.6)') &
            beam_data%turn(ix)%i_turn, bd%n_live, bd%time_sum/bd%n_count, norm2(bd%spin_sum/bd%n_count), &
            bd%spin_sum/bd%n_count, sig1, bd%orb_sum/bd%n_count, bd%p0c_sum/bd%n_count, &
            bunch_params%a%emit, bunch_params%b%emit, bunch_params%c%emit
  enddo

  close (iu)
enddo

end subroutine ltt_write_beam_averages

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_write_sigma_matrix (lttp, beam_data)

type (ltt_params_struct) lttp
type (ltt_beam_data_struct), target :: beam_data
type (ltt_bunch_data_struct), pointer :: bd

real(rp) sigma(21)
integer ix, i, j, k, n, iu, ib
character(200) :: file_name

! i_turn = -1 is used with mpi version.

if (lttp%sigma_matrix_output_file == '') return

!

iu = lunget()

n = size(beam_data%turn(1)%bunch)
do ib = 1, n
  ix = index(lttp%sigma_matrix_output_file, '#')

  if (n == 1) then
    file_name = lttp%sigma_matrix_output_file
  elseif (ix == 0) then
    file_name = trim(lttp%sigma_matrix_output_file) // int_str(ib)
  else
    file_name = lttp%sigma_matrix_output_file(1:ix-1) // int_str(ib) // lttp%sigma_matrix_output_file(ix+1:)
  endif

  if (beam_data%turn(1)%status == valid$) then
    open (iu, file = file_name, recl = 400)
    write (iu, '(a1, a8, a9, a12, 2x, 22a14)') '#', 'Turn', 'N_live', 'Time', '<p0c>', &
      '<x.x>', '<x.px>', '<x.y>', '<x.py>', '<x.z>', '<x.pz>', '<px.px>', '<px.y>', '<px.py>', '<px.z>', '<px.pz>', &
      '<y.y>', '<y.py>', '<y.z>', '<y.pz>', '<py.py>', '<py.z>', '<py.pz>', '<z.z>', '<z.pz>', '<pz.pz>'
  else
    open (iu, file = lttp%sigma_matrix_output_file, recl = 400, access = 'append')
  endif

  !

  do ix = 1, size(beam_data%turn)
    if (beam_data%turn(ix)%status /= valid$) cycle
    bd => beam_data%turn(ix)%bunch(ib)
    if (bd%n_count == 0) exit
    k = 0
    do i = 1, 6
    do j = i, 6
      k = k + 1
      sigma(k) = bd%orb2_sum(i,j) / bd%n_count - bd%orb_sum(i) * bd%orb_sum(j) / bd%n_count**2
    enddo
    enddo
    write (iu, '(i9, i9, es14.6, 2x, 22es14.6)') beam_data%turn(ix)%i_turn, bd%n_live, &
                                      bd%time_sum/bd%n_live, bd%p0c_sum/bd%n_count, (sigma(k), k = 1, 21)
  enddo

  close (iu)
enddo

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
  write (1, '(a, i0)') '# averaging_window = ', lttp%averaging_window
  write (1, '(a)') '# NO DATA TO AVERAGE OVER!'
  return
endif

!

average = average / n_sum
sigma = sigma / n_sum - outer_product(average, average)

write (1, '(a, i0)') '# averaging_window = ', lttp%averaging_window
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

subroutine ltt_read_map (lttp, ltt_com, err_flag)

type (ltt_params_struct) lttp
type (ltt_com_struct), target :: ltt_com
type (lat_struct), pointer :: lat
type (branch_struct), pointer :: branch
type (ltt_section_struct), pointer :: sec
type (ptc_map_with_rad_struct), pointer :: map

integer i, ix, n_sec, ib, ie, creation_hash
real(rp) s
logical err_flag, err, map_file_exists
character(200) map_file 

!

err_flag = .true.

call ltt_map_file_name (lttp, ltt_com, map_file)

lat => ltt_com%tracking_lat
branch => lat%branch(ltt_com%ix_branch)

inquire (file = map_file, exist = map_file_exists)
if (.not. map_file_exists) then
  if (ltt_com%mpi_rank == master_rank$) print '(a)', 'MAP FILE DOES NOT EXIST: ' // map_file
  return
endif

open (1, file = map_file, form = 'unformatted')
read (1, err = 9000, end = 9000) creation_hash

read (1, err = 9000, end = 9000) n_sec
if (allocated(ltt_com%sec)) deallocate(ltt_com%sec)
allocate (ltt_com%sec(0:n_sec))

do i = 0, n_sec
  sec => ltt_com%sec(i)
  read (1, err = 9000, end = 9000) ix, sec%type, ib, ie
  sec%ele => lat%branch(ib)%ele(ie)
  if (sec%type == map$) then
    allocate (sec%map)
    call ptc_read_map_with_radiation(sec%map, err, file_unit = 1)
    if (err) then
      print '(a)', 'ERROR READING MAP: ' // map_file
      return
    endif
    map => sec%map
  endif
enddo    

close (1)

if (ltt_com%mpi_rank == master_rank$) then
  if (ltt_com%lat%creation_hash == creation_hash) then
    print '(2a)', 'Map read in from file: ', trim(map_file)
    print '(2a)', 'Lattice file used for map: ', trim(map%lattice_file)
  else
    print '(a)', 'NOTE: LATTICE HAS BEEN MODIFIED SINCE MAP FILE ' // trim(map_file) // ' WAS CREATED.'
    print '(a)', '      WILL MAKE A NEW MAP.'
  endif
endif

if (ltt_com%lat%creation_hash == creation_hash) err_flag = .false.
return

9000 continue
close (1)
print '(a)', 'ERROR READING MAP FILE: ' // map_file

end subroutine ltt_read_map

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_write_map (lttp, ltt_com)

type (ltt_params_struct) lttp
type (ltt_com_struct), target :: ltt_com
type (ltt_section_struct), pointer :: sec

integer i, n_sec, ib, ie
character(200) map_file 

!

call ltt_map_file_name (lttp, ltt_com, map_file)

!

open (1, file = map_file, form = 'unformatted')
write (1) ltt_com%lat%creation_hash

do n_sec = ubound(ltt_com%sec, 1), 0, -1
  if (associated(ltt_com%sec(n_sec)%ele)) exit
enddo

write (1) n_sec

do i = 0, n_sec
  sec => ltt_com%sec(i)
  write (1) i, sec%type, sec%ele%ix_branch, sec%ele%ix_ele
  if (sec%type == map$) call ptc_write_map_with_radiation(sec%map, file_unit = 1)
enddo

close (1)

end subroutine ltt_write_map

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_map_file_name (lttp, ltt_com, map_file)

type (ltt_params_struct) lttp
type (ltt_com_struct), target :: ltt_com

integer hash

character(*) map_file
character(200) string
logical damping_on

!

if (lttp%map_file_prefix == '') then
  map_file = 'ltt'
else
  map_file = lttp%map_file_prefix
endif

damping_on = (bmad_com%radiation_damping_on .and. .not. lttp%split_bends_for_radiation) 
write (string, '(2a,2l1,i0,l1)') trim(lttp%exclude_from_maps), trim(lttp%ele_start), &
                                   damping_on, lttp%split_bends_for_radiation, lttp%map_order, lttp%rfcavity_on
if (lttp%simulation_mode == 'CHECK') string = trim(string) // lttp%ele_stop

hash = djb_hash(string)
write (string, '(i0)') hash
if (hash < 0) string(1:1) = 'n'    ! Looks nicer

map_file = trim(map_file) // '.' // trim(string) // '.map'

end subroutine ltt_map_file_name

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_make_map (lttp, ltt_com)

type (ltt_params_struct) lttp
type (ltt_com_struct), target :: ltt_com
type (lat_struct), pointer :: lat
type (ltt_section_struct), pointer :: sec
type (ele_struct), pointer :: ele0, ele, ele_start, ele_stop
type (branch_struct), pointer :: branch

real(rp) time0, time1, time_now
integer i, n_sec, ib, ie, ix_branch
logical err, in_map_section, damping_on

!

print '(a)', 'Creating map(s)...'
call run_timer ('ABS', time0)
time1 = time0

lat => ltt_com%tracking_lat
call ltt_pointer_to_map_ends (lttp, ltt_com%tracking_lat, ele_start, ele_stop)
branch => lat%branch(ele_start%ix_branch)
call lat_to_ptc_layout(lat)

i = 0
n_sec = 2
! Count slices and superimpose radiation points.
do i = 1, branch%n_ele_track
  ele => branch%ele(i)
  if (ele%key == sbend$) then
    n_sec = n_sec + 1
  elseif (ele%ix_pointer == 1) then
    n_sec = n_sec + 2
  endif
enddo

if (allocated(ltt_com%sec)) deallocate (ltt_com%sec)
allocate (ltt_com%sec(0:n_sec))
n_sec = 0
in_map_section = .false.
damping_on = bmad_com%radiation_damping_on
if (lttp%split_bends_for_radiation) damping_on = .false.

ltt_com%sec(0)%ele => ele_start
ltt_com%sec(0)%type = ele$

ele => ltt_com%sec(0)%ele
! Section setup
n_sec = 0
do i = 1, branch%n_ele_track+1
  ele0 => ele
  ele => pointer_to_next_ele(ele, skip_beginning = .false.)
  if (ele%key == marker$ .and. .not. in_map_section) cycle

  if (ele%ix_pointer == 1) then
    if (in_map_section) then
      n_sec = n_sec + 1
      allocate(ltt_com%sec(n_sec)%map)
      call ptc_setup_map_with_radiation (ltt_com%sec(n_sec)%map, ltt_com%sec(n_sec-1)%ele, ele0, &
                                              lttp%map_order, damping_on, lttp%symplectic_map_tracking)
      ltt_com%sec(n_sec)%type = map$
      ltt_com%sec(n_sec)%ele => ele
    endif

    n_sec = n_sec + 1
    ltt_com%sec(n_sec)%type = ele$
    ltt_com%sec(n_sec)%ele => ele
    in_map_section = .false.
  elseif (ele%key /= marker$ .or. ele%key /= beginning_ele$) then
    in_map_section = .true.
  endif

  if (ele%ix_ele == ele_stop%ix_ele) then
    if (in_map_section) then
      n_sec = n_sec + 1
      allocate(ltt_com%sec(n_sec)%map)
      call ptc_setup_map_with_radiation (ltt_com%sec(n_sec)%map, ltt_com%sec(n_sec-1)%ele, ele, &
                                             lttp%map_order, damping_on, lttp%symplectic_map_tracking)
      ltt_com%sec(n_sec)%type = map$
      ltt_com%sec(n_sec)%ele => ele
    endif
    exit
  endif

  call run_timer ('ABS', time_now)
  if (time_now - time1 > lttp%timer_print_dtime) then
    print '(a)', ' Map setup at (min) ' // real_str((time_now-time0)/60, 3) // ' at element: ' // int_str(i) // ' of: ' // int_str(branch%n_ele_track)
    time1 = time_now
  endif
enddo

call run_timer ('ABS', time_now)
print '(a, f8.2)', ' Map setup time (min)', (time_now-time0)/60

end subroutine ltt_make_map

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_make_tracking_lat (lttp, ltt_com)

type (ltt_params_struct) lttp
type (ltt_com_struct), target :: ltt_com
type (ele_struct) marker
type (ele_struct), pointer :: ele, slave
type (branch_struct), pointer :: branch
type (ele_pointer_struct), allocatable :: eles(:)

integer i, ie, is, n_loc
logical err

! Mark elements excluded from maps.
! Need to mark uncut lat in case element indexes are used.

call lat_ele_locator (lttp%exclude_from_maps, ltt_com%lat, eles, n_loc, err)
if (err) then
  print '(a)', 'ERROR LOCATING LATTICE ELEMENTS SPECIFIED IN LTT%EXCULDE_FROM_MAPS: ' // lttp%exclude_from_maps
  stop
endif

do ie = 1, n_loc
  eles(ie)%ele%ix_pointer = 1  ! Mark to exclude from any maps
enddo

!

ltt_com%tracking_lat = ltt_com%lat
branch => ltt_com%tracking_lat%branch(ltt_com%ix_branch)

if (lttp%split_bends_for_radiation) then
  call init_ele(marker, marker$)
  marker%name = 'RADIATION_PT'
  marker%ix_pointer = 1

  i = 0
  do
    i = i + 1
    if (i > branch%n_ele_track) exit
    ele => branch%ele(i)
    if (ele%key /= sbend$) cycle
    marker%s = (ele%s_start + ele%s) / 2.0_rp
    call add_superimpose(branch%lat, marker, branch%ix_branch, err)
    i = i + 2
  enddo

  call lattice_bookkeeper(ltt_com%tracking_lat)
endif

! Mark slaves of lords that are not to be part of any map

do ie = ltt_com%tracking_lat%n_ele_track+1, ltt_com%tracking_lat%n_ele_max
  ele => ltt_com%tracking_lat%ele(ie)
  if (ele%ix_pointer == 0) cycle
  do is = 1, ele%n_slave
    slave => pointer_to_slave(ele, is)
    slave%ix_pointer = 1
  enddo
enddo

end subroutine ltt_make_tracking_lat

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_pointer_to_map_ends (lttp, lat, ele_start, ele_stop)

type (ltt_params_struct) lttp
type (ltt_com_struct), target :: ltt_com
type (lat_struct) :: lat
type (ele_struct), pointer :: ele_start
type (ele_struct), pointer, optional :: ele_stop
type (ele_pointer_struct), allocatable :: eles(:)

integer n_loc
logical err


! @START_ELE and @STOP_ELE strings are set in ltt_init_params

if (lttp%ele_start == '') then
  ele_start => lat%ele(0)
else
  call lat_ele_locator ('type::@START_ELE', lat, eles, n_loc, err)
  ele_start => eles(1)%ele
endif

if (present(ele_stop)) then
  if (lttp%ele_stop == '' .or. lttp%ele_stop == lttp%ele_start) then
    ele_stop => ele_start
  else
    call lat_ele_locator ('type::@STOP_ELE', lat, eles, n_loc, err)
    ele_stop => eles(1)%ele
  endif
endif

end subroutine ltt_pointer_to_map_ends

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_track_map (lttp, ltt_com, orbit, is_lost)

type (ltt_params_struct) lttp
type (ltt_com_struct), target :: ltt_com
type (coord_struct) orbit
type (ltt_section_struct), pointer :: sec
type (branch_struct), pointer :: branch

integer i
logical is_lost, damping_on, fluct_on

!

damping_on = bmad_com%radiation_damping_on
fluct_on = bmad_com%radiation_fluctuations_on
if (lttp%split_bends_for_radiation) then
  damping_on = .false.
  fluct_on = .false.
endif

branch => ltt_com%tracking_lat%branch(0)


do i = 1, ubound(ltt_com%sec, 1)
  sec => ltt_com%sec(i)
  if (.not. associated(sec%ele)) exit

  if (allocated(sec%map)) then
    call ptc_track_map_with_radiation (orbit, sec%map, rad_damp = damping_on, rad_fluct = fluct_on)
    is_lost = (orbit_too_large(orbit) .or. abs(orbit%vec(1)) > lttp%ptc_aperture(1) .or. abs(orbit%vec(3)) > lttp%ptc_aperture(2))
    if (is_lost) exit
  elseif (sec%ele%name == 'RADIATION_PT') then
    call track1_radiation_center(orbit, pointer_to_next_ele(sec%ele, -1), pointer_to_next_ele(sec%ele), branch%param)
  else  ! EG: beambeam
    call track1 (orbit, sec%ele, branch%param, orbit)
  endif
enddo


end subroutine ltt_track_map

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_print_mpi_info (lttp, ltt_com, line, do_print)

type (ltt_params_struct) lttp
type (ltt_com_struct), target :: ltt_com

real(rp) time_now
character(*) line
character(20) time_str
logical, optional :: do_print

!

if (.not. logic_option(lttp%debug, do_print)) return

call run_timer ('ABS', time_now)
call date_and_time_stamp (time_str)
print '(a, f8.2, 2a, 2x, i0, 2a)', 'dTime:', (time_now-ltt_com%time_start)/60, &
                                        ' Now: ', time_str, ltt_com%mpi_rank, ': ', trim(line)

end subroutine ltt_print_mpi_info

end module
