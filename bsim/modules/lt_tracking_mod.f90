!+
! Module lt_tracking_mod
!
! Routines used by the long_term_tracking program.
! Also the dynamic_aperture program.
!-

module lt_tracking_mod

use beam_mod
use twiss_and_track_mod
use ptc_map_with_radiation_mod
use high_energy_space_charge_mod
use superimpose_mod, only: add_superimpose
use radiation_mod
use expression_mod
use mode3_mod, only: make_N
use random_mod
use s_fitting_new, only: probe, internal_state, track_probe, assignment(=), operator(+), default, spin0

implicit none

integer, parameter :: master_rank$  = 0
integer, parameter :: in_map$ = 0
integer, parameter :: not_in_map$ = 1
integer, parameter :: core_max$ = 40
integer, private :: i_loop

! Essentially: The ltt_params_struct holds user setable parameters while the ltt_com_struct holds
! parameters that are not setable.

type ltt_column_struct
  character(120) :: param = ''
  character(40) :: header_str = ''
  character(20) :: format = ''
end type

! User settable parameters

type ltt_params_struct
  character(20) :: simulation_mode = ''       ! CHECK, SINGLE, BEAM, STAT, INDIVIDUAL
  character(20) :: tracking_method = 'BMAD'   ! MAP, PTC, BMAD
  character(20) :: extraction_method = ''     ! or ON
  character(100) :: exclude_from_maps = 'beambeam::*'
  character(40) :: ele_start = ''
  character(40) :: ele_stop = ''
  character(40) :: ele_extract = ''
  character(100) :: ele_write_at = ''
  character(200) :: lat_file = ''
  character(200) :: beam_output_file = ''
  character(200) :: custom_output_file = ''
  character(200) :: map_file_prefix = ''
  character(200) :: map_ascii_output_file = ''
  character(200) :: averages_output_file = ''
  character(200) :: phase_space_output_file = ''
  character(200) :: action_angle_output_file = ''
  character(200) :: per_particle_output_file = ''
  character(8) :: ran_engine = ''
  type (ltt_column_struct) column(100)
  integer :: ix_particle_record = -1    ! Experimental
  integer :: ix_turn_record = -1        ! Experimental
  integer :: ix_turn_start = 0
  integer :: ix_turn_stop = -1
  integer :: n_turns = -1
  integer :: map_order = -1
  integer :: particle_output_every_n_turns = -1
  integer :: averages_output_every_n_turns = -1
  integer :: beam_output_every_n_turns = -1
  integer :: random_seed = 0
  integer :: output_only_last_turns = -1
  real(rp) :: random_sigma_cut = -1  ! If positive, cutoff for Gaussian random num generator.
  real(rp) :: core_emit_cutoff(core_max$) = [0.5_rp, (-1.0_rp, i_loop = 2, core_max$)]
  real(rp) :: ramping_start_time = 0
  real(rp) :: ptc_aperture(2) = 0.1_rp
  real(rp) :: print_on_dead_loss = -1
  real(rp) :: timer_print_dtime = 120
  real(rp) :: dead_cutoff = 1
  real(rp) :: a_emittance = 0   ! Used for space charge calculation.
  real(rp) :: b_emittance = 0   ! Used for space charge calculation.
  logical :: output_combined_bunches = .true.
  logical :: action_angle_calc_uses_1turn_matrix = .false.
  logical :: core_emit_combined_calc = .true.
  logical :: only_live_particles_out = .true.
  logical :: ramping_on = .false.
  logical :: ramp_update_each_particle = .false.
  logical :: ramp_particle_energy_without_rf = .false.
  logical :: rfcavity_on = .true.
  logical :: add_closed_orbit_to_init_position = .true.
  logical :: symplectic_map_tracking = .false.
  logical :: split_bends_for_stochastic_rad = .false.
  logical :: use_rf_clock = .false.
  logical :: debug = .false.
  logical :: regression_test = .false.          ! Only used for regression testing. Not of general interest.
  logical :: set_beambeam_crossing_time = .false.
  logical :: set_beambeam_z_crossing = .false.
  logical :: print_info_messages = .true.          ! Informational messages printed?
  !
  character(200) :: sigma_matrix_output_file = ''  ! No longer used
  character(200) :: master_output_file = ''        ! No longer used
  character(200) :: particle_output_file = ''      ! Replaced by phase_space_output_file
  character(200) :: beam_binary_output_file = ''   ! No longer used
  integer :: averaging_window = 1                  ! No longer used
  integer :: mpi_runs_per_subprocess = 4           ! No longer used
  logical :: split_bends_for_radiation = .false.   ! No longer used
end type

! A section structure is either:
!   1) A map between two points.
!   2) A pointer to an element where:
!         a) Radiation is to be put in. Or
!         b) The element is to be tracked.

integer, parameter :: map$ = 1, ele$ = 2

type ltt_section_struct
  integer :: type = 0
  type (ptc_rad_map_struct), allocatable :: map
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
  type (ele_struct), pointer :: ele_start_ptr            ! Pointer to start ele in original (not tracking) lattice.
  type (beam_init_struct) beam_init
  type (beam_init_struct) beam_init_used
  type (random_state_struct) ramper_ran_state
  integer :: n_ramper_loc = 0
  integer, allocatable :: ix_wake_ele(:)     ! List of element indexes where a wake is applied.
  integer :: ix_branch = 0                   ! Lattice branch being tracked.
  integer :: random_seed_actual = 0
  integer :: num_maps = 0
  real(rp) :: ptc_closed_orb(6) = 0
  real(rp) :: time_start = 0
  logical :: wrote_phase_space_file_header = .false.
  logical :: wrote_action_angle_file_header = .false.
  logical :: debug = .false.
  integer :: n_particle      ! Num particles per bunch. Needed with MPI.
  integer :: mpi_rank = master_rank$
  integer :: mpi_run_index = 0                    ! Run index
  integer :: mpi_ix0_particle = 0                 ! Index of first particle
  logical :: using_mpi = .false.
  logical :: track_bypass = .false.               ! Used by DA program
  logical :: ltt_tracking_happening_now = .false. ! Toggled True by program after all pre-tracking 
                                                  !   done (like the closed orbit calc). 
  character(200) :: master_input_file = ''
!!  character(200) :: last_beam_output_file = ''
  character(40) :: ps_fmt = '(2i7, 9es16.8, 3x, 3f13.9, 4x, a)'
end type

integer, parameter :: new$ = 0,  valid$ = 1, written$ = 2

type ltt_bunch_data_struct
  real(rp) :: n_live = 0        ! Number alive. Use real to avoid roundoff error.
  real(rp) :: sig1(6) = 0, sigma_vec(27) = 0
  real(rp) :: kurt(3) = 0, skew(3) = 0
  real(rp) :: core_emit(core_max$,3) = 0
  type (bunch_params_struct) :: params = bunch_params_struct()
end type

type (ltt_params_struct), pointer, save :: ltt_params_global   ! Needed for ltt_track1_preprocess and ltt_track1_bunch_hook
type (ltt_com_struct),    pointer, save :: ltt_com_global      ! Needed for ltt_track1_preprocess and ltt_track1_bunch_hook

contains

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_read_params(ltt, ltt_com)

type (ltt_params_struct), target :: ltt
type (ltt_com_struct), target :: ltt_com
type (beam_init_struct) beam_init

integer i, ix
character(200) arg
character(40) m_name
character(*), parameter :: r_name = 'ltt_read_params'

namelist / params / bmad_com, beam_init, ltt

!

ltt_params_global => ltt
ltt_com_global => ltt_com
ltt_com%master_input_file = ''
ltt_com%ltt_tracking_happening_now = .false.
ltt%ix_turn_start = int_garbage$
beam_init%ix_turn = int_garbage$

! Parse command line

i = 0
do while (i < command_argument_count())
  i = i + 1
  call get_command_argument(i, arg)
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

open(1, file = ltt_com%master_input_file, status = 'old', action = 'read')
read (1, nml = params)
close(1)

if (.not. ltt%print_info_messages) call output_direct (-1, .false., s_blank$, s_success$) ! Do not print
if (.not. ltt%print_info_messages) call output_direct (-1, .false., s_important$, s_important$) ! Do not print

if (.not. ltt_com%using_mpi .or. ltt_com%mpi_rank == master_rank$) then
  call out_io (s_blank$, r_name, 'Initialization file: ' // trim(ltt_com%master_input_file))
endif

if (ltt%ix_turn_start /= int_garbage$ .and. beam_init%ix_turn /= int_garbage$) then
  print *, 'ltt_com%ix_turn_start and beam_init%ix_turn are the same thing and cannot both be set. Set one or the other.'
  stop
elseif (beam_init%ix_turn /= int_garbage$) then
  ltt%ix_turn_start = beam_init%ix_turn
elseif (ltt%ix_turn_start /= int_garbage$) then
  beam_init%ix_turn = ltt%ix_turn_start
else
  ltt%ix_turn_start = 0
  beam_init%ix_turn = 0
endif

ltt_com%beam_init = beam_init

call bmad_parser (ltt%lat_file, ltt_com%lat)

! Read the master input file again so that bmad_com parameters set in the file
! take precedence over bmad_com parameters set in the lattice file.

open(1, file = ltt_com%master_input_file, status = 'old', action = 'read')
read (1, nml = params)  
close(1)

! Flag the setting of parameters that have been deprecated

if (ltt%beam_binary_output_file /= '') then
  call out_io (s_info$, r_name, &
        '"ltt%beam_binary_output_file" has been replaced by "ltt%beam_output_file". ', &
        'Beam file names with ".h5" or ".hdf5" suffix will be binary and otherwise will have an ASCII format.', &
        'Please change your init file! Stopping here.')
  stop
endif

if (ltt%split_bends_for_radiation) then
  call out_io (s_info$, r_name, &
        '"ltt%split_bends_for_radiation" has been renamed "ltt%split_bends_for_stochastic_rad".', &
        'And the manual has been updated to correctly describe what this parameter does.', &
        'Please change in the init file! Stopping here.')
  stop
endif

if (ltt%master_output_file /= '') then
  call out_io (s_info$, r_name, &
        'Note: The master_output_file is no longer created since data in this file is written to other files.', &
        'Please correct the init file. Will stop here.')
  stop
endif

if (ltt%sigma_matrix_output_file /= '') then
  call out_io (s_info$, r_name, &
        'Note: ltt%sigma_matrix_output_file is no longer used.', &
        '  Essentially now ltt%averages_output_file sets the name of the sigma matrix file.', &
        '  See manual for details.', &
        '  Please correct the init file. Will stop here.')
  stop
endif

if (ltt%averaging_window /= 1) then
  call out_io (s_info$, r_name, &
        'Note: ltt%averaging_window is no longer used.', &
        '  Reason: The advantage of averaging over multiple turns was not worth the complications', &
        '  to the program internal bookkeeping', &
        'Please correct the init file. Will stop here.')
  stop
endif

if (ltt%mpi_runs_per_subprocess /= 4) then
  call out_io (s_info$, r_name, 'Note: ltt%mpi_runs_per_subprocess is no longer used.', &
                                'Please correct the init file. Will stop here.')
  stop
endif

if (any(ltt%core_emit_cutoff > 1.00000001_rp .or. &
              (ltt%core_emit_cutoff < 0 .and. ltt%core_emit_cutoff /= -1.0_rp))) then
  call out_io (s_fatal$, r_name, 'ltt%core_emit_cutoff MUST BE GREATER THAN ZERO AND LESS THAN OR EQUAL TO ONE. STOPPING HERE.')
  stop
endif

if (ltt%particle_output_file /= '') then
  call out_io (s_error$, r_name, 'Note: ltt%particle_output_file replaced by ltt%phase_space_output_file.', &
                                 ' Also see ltt%action_angle_output_file.', &
                                 'Please correct the init file. Will stop here.')
  stop
endif

end subroutine ltt_read_params

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_init_params(ltt, ltt_com)

type (ltt_params_struct), target :: ltt
type (ltt_com_struct), target :: ltt_com
type (lat_struct), pointer :: lat
type (branch_struct), pointer :: branch
type (ele_pointer_struct), allocatable :: eles(:)
type (ele_struct), pointer :: ele
type (controller_struct), pointer :: cr
type (ele_pointer_struct), allocatable :: ramper(:)

real(rp) dummy
integer i, n, n_loc, ir, iv, it, is
character(*), parameter :: r_name = 'ltt_init_params'
logical err, found

! Lattice init

if (ltt%ramping_on) global_com%mp_threading_is_safe = .false.

lat => ltt_com%lat
bmad_com%auto_bookkeeper = .false.

call ran_seed_put (ltt%random_seed)
call ptc_ran_seed_put (ltt%random_seed)
call ran_seed_get (ltt_com%random_seed_actual)
call ran_gauss_converter(set_sigma_cut = ltt%random_sigma_cut)
if (ltt%ran_engine /= '') call ran_engine (set = ltt%ran_engine)

if (ltt_com%using_mpi) then
  call ran_seed_get (ir)
  call ran_seed_put (ir + ltt_com%mpi_rank)
  call ptc_ran_seed_put (ir + ltt_com%mpi_rank)
endif

if (.not. ltt%rfcavity_on) call set_on_off (rfcavity$, lat, off$)

! Sanity checks

call upcase_string(ltt%simulation_mode)

select case (ltt%simulation_mode)
case ('BEAM', 'STAT')
case ('CHECK', 'SINGLE', 'INDIVIDUAL')
  ltt%ramp_update_each_particle = .true.
case ('BUNCH')
  print '(a)', '"BUNCH" SETTING FOR LTT%SIMULATION_MODE HAS BEEN CHANGED TO "BEAM"'
case default
  print '(a)', 'UNKNOWN LTT%SIMULATION_MODE: ' // ltt%simulation_mode
  stop
end select

if (ltt%simulation_mode /= 'CHECK') then
  select case (ltt%tracking_method)
  case ('', 'MAP', 'PTC', 'BMAD', 'OLD')
  case default
    print '(a)', 'UNKNOWN LTT%TRACKING_METHOD: ' // ltt%tracking_method
    stop
  end select
endif

!

call upcase_string(ltt%extraction_method)
select case (ltt%extraction_method)
case ('', 'ON')
case default
  print '(a)', 'UNKNOWN LTT%EXTRACTION_METHOD: ' // ltt%extraction_method
  stop
end select

! Start/stop turn indexes

if (ltt%n_turns > 0 .and. ltt%ix_turn_stop > 0) then
  print *, 'Setting both ltt%n_turns and ltt%ix_turn_stop not allowed.'
  stop
elseif (ltt%ix_turn_stop > 0) then
  ltt%n_turns = ltt%ix_turn_stop - ltt%ix_turn_start
else
  ltt%ix_turn_stop = ltt%ix_turn_start + ltt%n_turns
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
  ltt_com%ele_start_ptr => eles(1)%ele
else
  ltt_com%ix_branch = 0
  ltt_com%ele_start_ptr => lat%branch(0)%ele(0)
endif

ltt_com%ele_start_ptr%type = '@START_ELE'
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

! Ramper setup

if (ltt%ramping_on) then
  if (ltt%tracking_method == 'PTC' .or. ltt%tracking_method == 'MAP') THEN
    print *, 'NOTE: ltt%tracking_method MUST BE SET TO "BMAD" IF LTT%RAMPING_ON IS SET TO TRUE.'
  endif

  call lat_ele_locator ('RAMPER::*', lat, ramper, ltt_com%n_ramper_loc, err)
  if (ltt_com%n_ramper_loc == 0) then
    print '(2a)', 'Warning! NO RAMPER ELEMENTS FOUND IN LATTICE.'
    stop
  endif

  do i = 1, ltt_com%n_ramper_loc
    ele => ramper(i)%ele
    cr => ele%control
    found = .false.

    do iv = 1, size(cr%var)
      if (cr%var(iv)%name == 'TIME') found = .true.
    enddo

    do is = 1, size(cr%ramp)
      if (.not. allocated(cr%ramp(is)%stack)) cycle
      do it = 1, size(cr%ramp(is)%stack)
        if (cr%ramp(is)%stack(it)%type == ran$ .or. cr%ramp(is)%stack(it)%type == ran_gauss$) found = .true.
      enddo
    enddo

    if (.not. found) then
      print *, 'Note! This ramper does not use "time" as a control variable and does not use any random functions: ' // trim(ele%name)
      print *, '      This ramper will not be directly varied in the simulation (but may be indirectly used).'
    endif
  enddo
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

call fullfilename (ltt%per_particle_output_file, ltt%per_particle_output_file)
call fullfilename (ltt%phase_space_output_file, ltt%phase_space_output_file)
call fullfilename (ltt%action_angle_output_file, ltt%action_angle_output_file)
call fullfilename (ltt%custom_output_file, ltt%custom_output_file)

end subroutine ltt_init_params

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_init_tracking (lttp, ltt_com, beam)

use f95_lapack, only: dpotrf_f95

type (ltt_params_struct) lttp
type (ltt_com_struct), target :: ltt_com
type (beam_struct), optional :: beam
type (lat_struct), pointer :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele
type (rad_map_ele_struct), pointer :: ri

real(rp) closed_orb(6), f_tol, time, ff

integer i, iv, n, ix_branch, ib, n_slice, ie, ir, info, n_rf_included, n_rf_excluded

logical err, map_file_exists, ramping_on

character(40) start_name, stop_name
character(*), parameter :: r_name = 'ltt_init_tracking'

!

if (lttp%split_bends_for_stochastic_rad .and. lttp%ramping_on) then
  call out_io (s_fatal$, r_name, 'LTT%SPLIT_BENDS_FOR_STOCHASTIC_RAD NOT COMPATIBLE WITH RAMPING.')
  stop
endif

! Apply rampers and then turn off ramper use for twiss_and_track since rampers simulating 
! noise (using random numbers) will drive the closed orbit calc crazy.

call ltt_make_tracking_lat(lttp, ltt_com)
lat => ltt_com%tracking_lat
branch => lat%branch(ltt_com%ix_branch)

ramping_on = lttp%ramping_on

if (lttp%ramping_on) then
  ! Using the simulation starting time for setting the ramper time to calculate the closed orbit is somewhat arbitrary.
  ele => ltt_com%ele_start_ptr
  time = 0.5_rp * (ele%ref_time + ele%value(ref_time_start$)) + &
           lttp%ix_turn_start * (branch%ele(branch%n_ele_track)%ref_time - branch%ele(0)%ref_time) + &
           ltt_params_global%ramping_start_time

  do ie = 0, branch%n_ele_max
    call ltt_apply_rampers_to_slave (ltt_com, branch%ele(ie), time, err)
  enddo

  lttp%ramping_on = .false.
endif  

call twiss_and_track (lat, ltt_com%bmad_closed_orb, ix_branch = ltt_com%ix_branch)
call radiation_integrals (lat, ltt_com%bmad_closed_orb, ltt_com%modes, ix_branch = ltt_com%ix_branch)

! PTC has an internal aperture of 1.0 meter. To be safe, set default aperture at 0.9 meter

lttp%ptc_aperture = min([0.9_rp, 0.9_rp], lttp%ptc_aperture)

if ((lttp%tracking_method == 'PTC' .or. lttp%tracking_method == 'MAP') .and. &
                                             abs(ltt_com%bmad_closed_orb(0)%vec(5)) > 0.8) then
  call out_io(s_fatal$, r_name, 'Z PHASE SPACE COORDINATE HAS MAGNITUDE TOO LARGE FOR PTC TRACKING (PTC HAS AN', &
                                'APERTURE LIMIT OF 1 METER). SOLUTION: SHIFT ALL RF PHI0 TO MAKE Z SMALLER.')
  stop
endif

! Setup RF clock if wanted.

if (lttp%use_rf_clock) then
  call out_io (s_info$, r_name, 'LTT%USE_RF_CLOCK NOT LONGER USED AND WILL BE IGNORED.')
endif

if (lttp%simulation_mode == 'BMAD' .and. .not. bmad_com%absolute_time_tracking) then
  call out_io (s_warn$, r_name, 'NOTE: It is recommended that bmad_com%absolute_time_tracking be set to True to', &
                                '      avoid problems with frequencies incommensurate with the revolution frequency.')
endif

!

if (lttp%split_bends_for_stochastic_rad) then
  do ie = 1, branch%n_ele_track
    ele => branch%ele(ie)
    if (ele%name /= 'RADIATION_POINT') cycle
    if (.not. associated(ele%rad_map)) allocate (ele%rad_map)
    ri => ele%rad_map
    call tracking_rad_map_setup(branch%ele(ie-1), 1e-4_rp, downstream_end$, ri%rm0, err)
    call tracking_rad_map_setup(branch%ele(ie+1), 1e-4_rp, upstream_end$,   ri%rm1, err)

    ! Combine matrices for quicker evaluation
    ri%rm1%stoc_mat = matmul(ri%rm0%stoc_mat, transpose(ri%rm0%stoc_mat)) + matmul(ri%rm1%stoc_mat, transpose(ri%rm1%stoc_mat))
    call dpotrf_f95 (ri%rm1%stoc_mat, 'L', info = info)
    if (info /= 0) then
      ri%rm1%stoc_mat = 0  ! Cholesky failed
    endif

    do i = 2, 6
      ri%rm1%stoc_mat(1:i-1, i) = 0
    enddo

    ri%rm1%xfer_damp_mat = ri%rm0%xfer_damp_mat + ri%rm1%xfer_damp_mat
    ri%rm1%xfer_damp_vec = ri%rm0%xfer_damp_vec + ri%rm1%xfer_damp_vec
  enddo
endif

!

if (lttp%tracking_method == 'PTC' .or. lttp%simulation_mode == 'CHECK') then
  if (.not. associated(lat%branch(0)%ptc%m_t_layout)) call lat_to_ptc_layout(lat)
  call ptc_setup_tracking_with_damping_and_excitation(lat%branch(0), bmad_com%radiation_damping_on, &
                                  bmad_com%radiation_fluctuations_on, ltt_com%ptc_state, ltt_com%ptc_closed_orb)
  if (bmad_com%spin_tracking_on) ltt_com%ptc_state = ltt_com%ptc_state + SPIN0
endif


if (lttp%set_beambeam_z_crossing) then
  call out_io(s_error$, r_name, '"lttp%set_beambeam_z_crossing" is now "lttp%set_beambeam_crossing_time".', &
                              'Please change your init file.')
  stop
endif

if (lttp%set_beambeam_crossing_time) then
  do ie = 1, branch%n_ele_track
    ele => branch%ele(ie)
    if (ele%key /= beambeam$) cycle
    ff = -1.0_rp / (c_light * ltt_com%bmad_closed_orb(ie)%beta)
    if (lttp%tracking_method == 'PTC') then
      ele%value(crossing_time$) = ff * (ltt_com%bmad_closed_orb(ie)%vec(5) + (ltt_com%ptc_closed_orb(5) - ltt_com%bmad_closed_orb(0)%vec(5)))
    else
      ele%value(crossing_time$) = ff * ltt_com%bmad_closed_orb(ie)%vec(5)
    endif
  enddo
endif

!

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

  do i = 0, ubound(ltt_com%sec, 1)
    if (allocated(ltt_com%sec(i)%map)) ltt_com%num_maps = ltt_com%num_maps + 1
  enddo
endif

lttp%ramping_on = ramping_on 

! Beam setup

if (present(beam) .and. (lttp%simulation_mode == 'INDIVIDUAL' .or. lttp%simulation_mode == 'BEAM')) then
  call ltt_init_beam_distribution(lttp, ltt_com, beam)
endif

end subroutine ltt_init_tracking

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_write_per_particle_file_header(lttp, ltt_com, beam)

type (ltt_params_struct) lttp
type (ltt_com_struct), target :: ltt_com
type (beam_struct), target :: beam

integer ib, ip, iu, np
character(200) file
character(*), parameter :: r_name = 'ltt_write_per_particle_file_header'
!

if (lttp%per_particle_output_file == '') return

iu = lunget()

do ib = 1, size(beam%bunch)
  np = size(beam%bunch(ib)%particle)
  if (np > 1000) then
    call out_io (s_fatal$, r_name, 'IT DOES NOT MAKE SENSE TO HAVE PER_PARTICLE OUTPUT WITH THIS NUMBER OF PARTICLES: ' // int_str(np), &
                                   'WILL STOP HERE')
    stop
  endif

  do ip = 1, np
    file = ltt_per_particle_file_name(lttp, ib, ip, size(beam%bunch), np)
    open(iu, file = file, recl = 300)
    call ltt_write_preheader(iu, lttp, ltt_com, .false.)
    write (iu, '(a, i0)') '# ix_bunch        = ', ib
    write (iu, '(a, i0)') '# ix_particle     = ', ip
    write (iu, '(a)')  '##   Turn |           x              px               y              py               z              pz              pc             p0c            time   |    spin_x      spin_y      spin_z     State'
    close(iu)
  enddo
enddo

end subroutine ltt_write_per_particle_file_header

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

function ltt_per_particle_file_name (lttp, ix_bunch, ix_particle, n_bunch, n_particle) result (file_name)

type (ltt_params_struct) lttp

integer ix_bunch, ix_particle, n_bunch, n_particle, np, ih, ix
character(200) file_name
character(16) fmt, str

!

np = int(log10(1.001*n_particle)) + 1

if (n_bunch == 1) then
  fmt = '(i' // int_str(np) // '.' // int_str(np) // ')'
  write (str, fmt) ix_particle
else
  fmt = '(2a, i' // int_str(np) // '.' // int_str(np) // ')'
  write (str, fmt) int_str(ix_bunch), '-', ix_particle
endif

ih = index(lttp%per_particle_output_file, '#')
if (ih == 0) then
  file_name = trim(lttp%per_particle_output_file) // '.' // str
else
  file_name = lttp%per_particle_output_file(:ih-1) // trim(str) // lttp%per_particle_output_file(ih+1:)
endif

end function ltt_per_particle_file_name

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_write_averages_header (lttp, ltt_com, extracting)

type (ltt_params_struct) lttp
type (ltt_com_struct), target :: ltt_com
type (branch_struct), pointer :: branch

integer i, ix, iu, nb, ib
logical, optional :: extracting

character(8) t_str
character(200) :: file_in, file_name
character(2000) :: line1, line2

!

nb = max(1, ltt_com%beam_init%n_bunch)
iu = lunget()
if (lttp%averages_output_file == '') return

if (logic_option(.false., extracting)) then
  t_str = 'S_pos'
  file_in = 'extraction-' // lttp%averages_output_file
else
  t_str = 'Turn'
  file_in = lttp%averages_output_file
endif

!

do ib = 0, nb
  if (ib == 0 .and. nb == 1) cycle  ! ib = zero is for all bunch averages if there is more than one bunch
  call ltt_averages_file_name(file_in, 'ave', ib, nb, file_name)
  open(iu, file = file_name, recl = 400)

  call ltt_write_preheader(iu, lttp, ltt_com, .false.)

  write (iu, '(a2, a8, a9, 2a14, 3a14, 2x, 8a14, 2x, 6a14)') '##', '1', '2', '3', '4', '5', '6', '7', &
                  '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21'
  write (iu, '(a2, a8, a9, 2a14, 3a14, 2x, 8a14, 2x, 6a14)') '##', trim(t_str), 'N_live', 'Time', 'Polarization', &
                   '<Sx>', '<Sy>', '<Sz>', '<x>', '<px>', '<y>', '<py>', '<z>', '<pz>', '<pc>', '<p0c>', &
                   'Sig_x', 'Sig_px', 'Sig_y', 'Sig_py', 'Sig_z', 'Sig_pz'

  close(iu)

  !

  call ltt_averages_file_name(file_in, 'sigma', ib, nb, file_name)
  open(iu, file = file_name, recl = 400)

  call ltt_write_preheader(iu, lttp, ltt_com, .false.)

  write (iu, '(a2, a8, a9, a12, 2x, 21a14)') '##', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', &
                  '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24'

  write (iu, '(a2, a8, a9, a12, 2x, 21a14)') '##', trim(t_str), 'N_live', 'Time', &
    '<x.x>', '<x.px>', '<x.y>', '<x.py>', '<x.z>', '<x.pz>', '<px.px>', '<px.y>', '<px.py>', '<px.z>', '<px.pz>', &
    '<y.y>', '<y.py>', '<y.z>', '<y.pz>', '<py.py>', '<py.z>', '<py.pz>', '<z.z>', '<z.pz>', '<pz.pz>'

  close(iu)

  !

  call ltt_averages_file_name(file_in, 'emit', ib, nb, file_name)
  open(iu, file = file_name, recl = 2000)

  call ltt_write_preheader(iu, lttp, ltt_com, .false.)

  write (line1, '(a2, a8, a10, a12, 2x, 3a14, 2x, 3a14, 2x, 3a14, a)') '##', '1', '2', '3', &
                '4', '5', '6', '7', '8', '9', '10', '11', '12', '  |'

  write (line2, '(a2, a8, a10, a12, 2x, 3a14, 2x, 3a14, 2x, 3a14, a)') '##', trim(t_str), 'N_live', 'Time', &
                'Emit_a', 'Emit_b', 'Emit_c', 'Kurtosis_x', 'Kurtosis_y', 'Kurtosis_z', 'Skew_x', 'Skew_y', 'Skew_z', '  |'


  do i = 1, core_max$
    if (lttp%core_emit_cutoff(i) <= 0) exit
    line1 = trim(line1) // '                   ZZ% Core                |'
    line2 = trim(line2) // '      Emit_a        Emit_b        Emit_c   |'
    ix = index(line1, 'ZZ')
    line1 = line1(1:ix-1) // int_str(nint(lttp%core_emit_cutoff(i)*100)) // line1(ix+2:)
  enddo

  write (iu, '(a)') trim(line1)
  write (iu, '(a)') trim(line2)

  close(iu)
enddo

end subroutine ltt_write_averages_header

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_write_preheader(iu, lttp, ltt_com, print_this)

type (ltt_params_struct) lttp
type (ltt_com_struct), target :: ltt_com
type (branch_struct), pointer :: branch
type(beam_init_struct), pointer :: bi

integer iu, species
real(rp) e_tot, a_gam, t0
logical print_this

!

branch => ltt_com%tracking_lat%branch(ltt_com%ix_branch)
species = branch%ele(0)%ref_species
e_tot = branch%ele(0)%value(e_tot$)
a_gam = anomalous_moment_of(species) * e_tot / mass_of(species)
t0 = branch%ele(branch%n_ele_track)%ref_time
call ltt_write_line('# e_tot                                   = ' // real_str(e_tot, 6), lttp, iu, print_this)
call ltt_write_line('# t_1turn                                 = ' // real_str(t0, 6), lttp, iu, print_this)
call ltt_write_line('# anom_moment_times_gamma                 = ' // real_str(a_gam, 6), lttp, iu, print_this)
call ltt_write_line('# master_input_file                       = ' // quote(ltt_com%master_input_file), lttp, iu, print_this)
call ltt_write_line('# ltt%lat_file                            = ' // quote(lttp%lat_file), lttp, iu, print_this)
call ltt_write_line('# ltt%averages_output_file                = ' // quote(lttp%averages_output_file), lttp, iu, print_this)
call ltt_write_line('# ltt%beam_output_file                    = ' // quote(lttp%beam_output_file), lttp, iu, print_this)
call ltt_write_line('# ltt%custom_output_file                  = ' // quote(lttp%custom_output_file), lttp, iu, print_this)
call ltt_write_line('# ltt%per_particle_output_file            = ' // quote(lttp%per_particle_output_file), lttp, iu, print_this)
call ltt_write_line('# ltt%phase_space_output_file             = ' // quote(lttp%phase_space_output_file), lttp, iu, print_this)
call ltt_write_line('# ltt%action_angle_output_file            = ' // quote(lttp%action_angle_output_file), lttp, iu, print_this)
call ltt_write_line('# ltt%map_file_prefix                     = ' // quote(lttp%map_file_prefix), lttp, iu, print_this)
call ltt_write_line('# ltt%simulation_mode                     = ' // quote(lttp%simulation_mode), lttp, iu, print_this)
call ltt_write_line('# ltt%tracking_method                     = ' // quote(lttp%tracking_method), lttp, iu, print_this)
call ltt_write_line('# ltt%ele_start                           = ' // quote(lttp%ele_start), lttp, iu, print_this)

call ltt_write_line('# ltt%extraction_method                   = ' // quote(lttp%extraction_method), lttp, iu, print_this)
if (lttp%extraction_method /= '') then
  call ltt_write_line('# ltt%ele_extract                         = ' // quote(lttp%ele_extract), lttp, iu, print_this)
  call ltt_write_line('# ltt%ele_write_at                        = ' // quote(lttp%ele_write_at), lttp, iu, print_this)
endif

if (lttp%simulation_mode == 'CHECK') then
  call ltt_write_line('# ltt%ele_stop                            = ' // quote(lttp%ele_stop), lttp, iu, print_this)
endif

if (lttp%tracking_method == 'MAP' .or. lttp%simulation_mode == 'CHECK') then
  call ltt_write_line('# ltt%map_order                           = ' // int_str(lttp%map_order), lttp, iu, print_this)
  call ltt_write_line('# ltt%exclude_from_maps                   = ' // quote(lttp%exclude_from_maps), lttp, iu, print_this)
  call ltt_write_line('# ltt%symplectic_map_tracking             = ' // logic_str(lttp%symplectic_map_tracking), lttp, iu, print_this)
  call ltt_write_line('# Number_of_maps                          = ' // int_str(ltt_com%num_maps), lttp, iu, print_this)
endif

call ltt_write_line('# ltt%split_bends_for_stochastic_rad      = ' // logic_str(lttp%split_bends_for_stochastic_rad), lttp, iu, print_this)
call ltt_write_line('# ltt%core_emit_combined_calc             = ' // logic_str(lttp%core_emit_combined_calc), lttp, iu, print_this)
call ltt_write_line('# ltt%action_angle_calc_uses_1turn_matrix = ' // logic_str(lttp%action_angle_calc_uses_1turn_matrix), lttp, iu, print_this)

if (bmad_com%sr_wakes_on) then
  call ltt_write_line('# Number_of_wake_elements                 = ' // int_str(size(ltt_com%ix_wake_ele)), lttp, iu, print_this)
endif

call ltt_write_line('# ltt%random_sigma_cut                    = ' // real_str(lttp%random_sigma_cut, 6), lttp, iu, print_this)
call ltt_write_line('# ltt%n_turns                             = ' // int_str(lttp%n_turns), lttp, iu, print_this)
call ltt_write_line('# ltt%ix_turn_start                       = ' // int_str(lttp%ix_turn_start), lttp, iu, print_this)
call ltt_write_line('# ltt%ix_turn_stop                        = ' // int_str(lttp%ix_turn_stop), lttp, iu, print_this)
call ltt_write_line('# ltt%beam_output_every_n_turns           = ' // int_str(lttp%beam_output_every_n_turns), lttp, iu, print_this)
call ltt_write_line('# ltt%particle_output_every_n_turns       = ' // int_str(lttp%particle_output_every_n_turns), lttp, iu, print_this)
call ltt_write_line('# ltt%averages_output_every_n_turns       = ' // int_str(lttp%averages_output_every_n_turns), lttp, iu, print_this)
call ltt_write_line('# ltt%output_only_last_turns              = ' // int_str(lttp%output_only_last_turns), lttp, iu, print_this)
call ltt_write_line('# ltt%ramping_on                          = ' // logic_str(lttp%ramping_on), lttp, iu, print_this)
call ltt_write_line('# ltt%output_combined_bunches             = ' // logic_str(lttp%output_combined_bunches), lttp, iu, print_this)
call ltt_write_line('# ltt%ramp_update_each_particle           = ' // logic_str(lttp%ramp_update_each_particle), lttp, iu, print_this)
call ltt_write_line('# ltt%ramp_particle_energy_without_rf     = ' // logic_str(lttp%ramp_particle_energy_without_rf), lttp, iu, print_this)
call ltt_write_line('# ltt%ramping_start_time                  = ' // real_str(lttp%ramping_start_time, 6), lttp, iu, print_this)
call ltt_write_line('# ltt%set_beambeam_crossing_time          = ' // logic_str(lttp%set_beambeam_crossing_time), lttp, iu, print_this)
call ltt_write_line('# ltt%random_seed                         = ' // int_str(lttp%random_seed), lttp, iu, print_this)

if (lttp%random_seed == 0) then
  call ltt_write_line('# random_seed_actual                      = ' // int_str(ltt_com%random_seed_actual), lttp, iu, print_this)
endif

call ltt_write_line('# ltt%rfcavity_on                         = ' // logic_str(lttp%rfcavity_on), lttp, iu, print_this)
call ltt_write_line('# is_RF_on                                = ' // logic_str(rf_is_on(branch)) // '  #  M65 /= 0 ?', lttp, iu, print_this)
call ltt_write_line('# bmad_com%radiation_damping_on           = ' // logic_str(bmad_com%radiation_damping_on), lttp, iu, print_this)
call ltt_write_line('# bmad_com%radiation_fluctuations_on      = ' // logic_str(bmad_com%radiation_fluctuations_on), lttp, iu, print_this)
call ltt_write_line('# bmad_com%spin_tracking_on               = ' // logic_str(bmad_com%spin_tracking_on), lttp, iu, print_this)
call ltt_write_line('# bmad_com%sr_wakes_on                    = ' // logic_str(bmad_com%sr_wakes_on), lttp, iu, print_this)

bi => ltt_com%beam_init_used
if (bi%a_emit /= real_garbage$) then
  call ltt_write_line('# a_emit                                  = ' // real_str(bi%a_emit, 6), lttp, iu, print_this)
  call ltt_write_line('# b_emit                                  = ' // real_str(bi%b_emit, 6), lttp, iu, print_this)
  call ltt_write_line('# a_norm_emit                             = ' // real_str(bi%a_norm_emit, 6), lttp, iu, print_this)
  call ltt_write_line('# b_norm_emit                             = ' // real_str(bi%a_norm_emit, 6), lttp, iu, print_this)
  call ltt_write_line('# sig_z                                   = ' // real_str(bi%sig_z, 6), lttp, iu, print_this)
  call ltt_write_line('# sig_pz                                  = ' // real_str(bi%sig_pz, 6), lttp, iu, print_this)
endif
call ltt_write_line('#--------------------------------------', lttp, iu, print_this)

end subroutine ltt_write_preheader

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_write_line (line, lttp, iu, print_this)

type (ltt_params_struct), optional :: lttp
integer :: iu
integer iv
logical, optional :: print_this
character(*) line
character(*), parameter :: r_name = 'ltt_write_line'

!

if (logic_option(.true., print_this)) call out_io(s_blank$, r_name, trim(line))
if (iu /= 0) write (iu, '(a)') trim(line)

end subroutine ltt_write_line

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_run_check_mode (lttp, ltt_com)

type (ltt_params_struct) lttp
type (ltt_com_struct), target :: ltt_com
type (lat_struct), pointer :: lat
type (coord_struct) orb_map, orb_map1, orb_map2, orb_bmad, orb_bmad1, orb_bmad2, orb_init
type (ele_struct), pointer :: ele_start, ele_stop
type (probe) prb_ptc, prb_ptc1, prb_ptc2

real(rp) mat_map(6,6), mat_bmad(6,6), mat_ptc(6,6), spin_ptc(3)
integer i
integer ix_branch, ix_start, ix_stop

! Run serial in check mode.

bmad_com%radiation_fluctuations_on = .false.
ltt_com%ltt_tracking_happening_now = .true.    ! Start of main tracking
lat => ltt_com%tracking_lat

call ltt_pointer_to_map_ends(lttp, lat, ele_start, ele_stop)
ix_branch = ltt_com%ix_branch
ix_start = ele_start%ix_ele
ix_stop  = ele_stop%ix_ele

call ltt_init_coord(lttp, ltt_com, orb_init, ele_start)

if (lttp%add_closed_orbit_to_init_position) orb_init%vec = orb_init%vec + ltt_com%bmad_closed_orb(ix_start)%vec

call track_check_all (orb_init, 0, 0, orb_map, orb_bmad, prb_ptc, spin_ptc, ele_start, ele_stop)

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
  call track_check_all (orb_init, i, -1, orb_map1, orb_bmad1, prb_ptc1, spin_ptc, ele_start, ele_stop)
  call track_check_all (orb_init, i, +1, orb_map2, orb_bmad2, prb_ptc2, spin_ptc, ele_start, ele_stop)
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

subroutine track_check_all (orb_init, i_ps, sgn, orb_map, orb_bmad, prb_ptc, spin_ptc, ele_start, ele_stop)

type (coord_struct) orb_init, orb_map, orb_bmad, orb_start
type (probe) prb_ptc
type (lat_struct), pointer :: lat
type (ele_struct), pointer :: ele_start, ele_stop

real(rp) spin_ptc(3)
integer i_ps, sgn

!

lat => ltt_com%tracking_lat
orb_start = orb_init
if (i_ps > 0) orb_start%vec(i_ps) = orb_start%vec(i_ps) + sgn * bmad_com%d_orb(i)

! Map

orb_map = orb_start
call ltt_track_map (lttp, ltt_com, orb_map)

! Bmad

orb_bmad = orb_start
call ltt_track_bmad_single (lttp, ltt_com, ele_start, ele_stop, orb_bmad)

! PTC

prb_ptc = orb_start%vec
prb_ptc%q%x = [1, 0, 0, 0]   ! Unit quaternion
if (ix_stop == ix_start) then
  call track_probe (prb_ptc, ltt_com%ptc_state, fibre1 = pointer_to_fibre(ele_start))
else
  call track_probe (prb_ptc, ltt_com%ptc_state, fibre1 = pointer_to_fibre(ele_start), &
                                                fibre2 = pointer_to_fibre(ele_stop))
endif
spin_ptc = quat_rotate(prb_ptc%q%x, orb_start%spin)

end subroutine track_check_all

end subroutine ltt_run_check_mode

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
! Note: When running with MPI, the slaves will call ltt_run_single_mode directly.
! That is, this routine is not called when using MPI.

subroutine ltt_run_individual_mode (lttp, ltt_com, beam)

type (ltt_params_struct) lttp
type (ltt_com_struct), target :: ltt_com
type (beam_struct), target :: beam
type (coord_struct), pointer :: particle
type (coord_struct) orbit

integer ib, ip

!

ltt_com%ltt_tracking_happening_now = .true.    ! Start of main tracking
call ltt_write_per_particle_file_header(lttp, ltt_com, beam)

do ib = 1, size(beam%bunch)
  do ip = 1, size(beam%bunch(ib)%particle)
    particle => beam%bunch(ib)%particle(ip)
    call ltt_run_single_mode(lttp, ltt_com, particle, beam, ib, ip)
  enddo
enddo

!

call ltt_write_particle_data (lttp, ltt_com, lttp%n_turns, beam)
if (lttp%beam_output_file /= '') call write_beam_file (lttp%beam_output_file, beam, .true.) ! Not ltt_write_beam_file

end subroutine ltt_run_individual_mode

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_run_extraction (lttp, ltt_com, beam)

type (ltt_params_struct) lttp
type (ltt_com_struct), target :: ltt_com
type (beam_struct), target :: beam
type (bunch_struct), pointer :: bunch
type (branch_struct), pointer :: branch
type (coord_struct), pointer :: p
type (ele_pointer_struct), allocatable :: eles(:)
type (ele_struct), pointer :: septum, fork, ele
type (lat_struct), pointer :: lat

real(rp) s_pos
integer ib, ip, ii, nt, n_loc, physical_end
logical err
character(200) file

!

if (lttp%extraction_method == '') return   ! No extraction wanted

lat => ltt_com%tracking_lat
branch => lat%branch(ltt_com%ix_branch)

if (lttp%beam_output_file /= '') then
  file = 'extraction-start-' // trim(lttp%beam_output_file)
  call write_beam_file (file, beam, .true.)  ! Note: Not ltt_write_beam_file
endif

! With INDIVIDUAL mode, only particles that have hit the septum aperture are to be
! tracked through the extraction line.

if (lttp%simulation_mode == 'INDIVIDUAL') then
  call lat_ele_locator(lttp%ele_extract, lat, eles, n_loc, err, .true.)
  if (n_loc == 0) then
    print '(a)', 'CANNOT FIND LTT%ELE_EXTRACT ELEMENT: ' // quote(lttp%ele_extract)
    stop
  endif

  if (n_loc > 1) then
    print '(a)', 'MULTIPLE ELEMENTS MATCH LTT%ELE_EXTRACT NAME'
    stop
  endif

  septum => eles(1)%ele

  if (septum%ix_branch /= ltt_com%ix_branch) then
    print '(a)', 'SEPTUM ELEMENT NOT IN RING TRACKING BRANCH!'
    stop
  endif

  ! Particles that have hit the aperture at the septum are resurrected and vice versa.

  do ib = 1, size(beam%bunch)
    bunch => beam%bunch(ib)
    do ip = 1, size(bunch%particle)
      p => bunch%particle(ip)
      if (p%state /= alive$ .and. p%ix_ele == septum%ix_ele) then
        p%state = alive$

        if (p%location == upstream_end$) then
          bmad_com%aperture_limit_on = .false.   ! Otherwise particle rehits upstream aperture.
          call track1(p, septum, branch%param, p)
          bmad_com%aperture_limit_on = .true.

          physical_end = physical_ele_end (second_track_edge$, p, septum%orientation)
          if (at_this_ele_end (physical_end, septum%aperture_at)) then
            call check_aperture_limit(p, septum, second_track_edge$, branch%param)
            if (p%state == alive$) then  ! Passed through septum wall from outside to inside!
              p%state = lost$           
            else
              p%state = alive$
            endif
          endif

        else
          if (septum%aperture_at == both_ends$) p%state = lost$  ! particle went through septum wall
        endif

      elseif (p%state == alive$) then   ! Did not make it to the septum.
        p%state = lost$
      endif
    enddo
  enddo

  ele => septum  ! Starting point for extraction tracking

! Else with bunch mode just start tracking from where we left off
else
  call ltt_pointer_to_map_ends(lttp, lat, ele)  ! ele => starting point for extraction tracking
  if (ele%key == fork$) ele => pointer_to_next_ele(ele, -1)
endif

! What places to print at?

call lat_ele_locator(lttp%ele_write_at, lat, eles, n_loc, err, .true.)

! Find where fork is so relative s-positions can be calculated

fork => pointer_to_next_ele(ele, skip_beginning = .true.)
do ii = 1, branch%n_ele_track
  if (fork%key == fork$) exit
  fork => pointer_to_next_ele(fork, skip_beginning = .true.)
enddo

!

do ib = 1, size(beam%bunch)
  call remove_dead_from_bunch(beam%bunch(ib), beam%bunch(ib))
enddo

! Track to the fork

call ltt_write_averages_header (lttp, ltt_com, .true.)

do ii = 1, branch%n_ele_track
  call ltt_extraction_write(lttp, ltt_com, beam, ele, eles(1:n_loc), ele%s - fork%s)
  ele => pointer_to_next_ele(ele, skip_beginning = .true.)
  if (ele%key == fork$) exit
  call track1_beam(beam, ele, err)
enddo

if (ii > branch%n_ele_track) then
  print '(a)', 'NO FORK ELEMENT IN TRACKING BRANCH!'
  stop
endif

! Extraction branch.

s_pos = 0

fork_loop: do   ! Loop in case extraction branch forks at the end to yet another branch, etc.
  branch => lat%branch(nint(ele%value(ix_to_branch$)))
  ele => branch%ele(nint(ele%value(ix_to_element$)))   ! Starting point in extraction branch.
  nt = branch%n_ele_track

  do
    call ltt_extraction_write(lttp, ltt_com, beam, ele, eles(1:n_loc), s_pos)
    if (ele%ix_ele == nt) exit fork_loop

    ele => pointer_to_next_ele(ele)
    if (ele%key == fork$ .and. (ele%ix_ele == nt .or. &
                (ele%ix_ele == nt-1 .and. branch%ele(nt)%key == marker$))) cycle fork_loop  ! Fork to new branch
    call track1_beam(beam, ele, err)
    s_pos = s_pos + ele%value(l$)
  enddo
enddo fork_loop

end subroutine ltt_run_extraction

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_extraction_write (lttp, ltt_com, beam, ele, eles, s_pos)

type (ltt_params_struct) lttp
type (ltt_com_struct), target :: ltt_com
type (beam_struct), target :: beam
type (bunch_struct), pointer :: bunch
type (branch_struct), pointer :: branch
type (coord_struct), pointer :: p
type (ele_pointer_struct) :: eles(:)
type (ele_struct), pointer :: septum, ele

real(rp) s_pos
integer i, ix
logical found
character(200) file

! Always add to averages file

call ltt_write_averages_data (lttp, -1, beam, s_pos, ele)

! Beam files

found = .false.
do i = 1, size(eles)
  if (associated(eles(i)%ele, ele)) found = .true.
enddo
if (.not. found) return

!

if (lttp%beam_output_file == '') return
file = int_str(ele%ix_branch) //  '.' // int_str(ele%ix_ele) // '_' // trim(ele%name) // '-' // trim(lttp%beam_output_file)
call fullfilename(file, file)
call write_beam_file (file, beam, .true.)

end subroutine ltt_extraction_write

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
! orb_in, etc. args are used with individual mode

subroutine ltt_run_single_mode (lttp, ltt_com, orb_in, beam, ix_bunch, ix_particle)

type (ltt_params_struct) lttp
type (ltt_com_struct), target :: ltt_com
type (beam_struct), optional :: beam
type (branch_struct), pointer :: branch
type (lat_struct), pointer :: lat
type (coord_struct), allocatable :: orb(:)
type (coord_struct) :: orbit
type (coord_struct), optional :: orb_in
type (ele_struct), pointer :: ele_start, ele0, ele1
type (probe) prb

real(rp) average(6), sigma(6,6)
integer i, n_sum, iu_part, i_turn, ix_branch
integer, optional :: ix_bunch, ix_particle
logical do_write
character(*), parameter :: r_name = 'ltt_run_single_mode'

! Run a single particle in single mode.

ltt_com%ltt_tracking_happening_now = .true.    ! Start of main tracking
lat => ltt_com%tracking_lat

n_sum = 0
sigma = 0
average = 0
ix_branch = ltt_com%ix_branch
branch => lat%branch(ix_branch)
call ltt_pointer_to_map_ends(lttp, lat, ele_start)

call ltt_setup_high_energy_space_charge(lttp, ltt_com, branch)

if (lttp%tracking_method == 'BMAD') call reallocate_coord (orb, lat)

if (present(orb_in)) then
  orbit = orb_in
else
  call ltt_init_coord(lttp, ltt_com, orbit, ele_start)
endif

!

if (lttp%simulation_mode == 'INDIVIDUAL') then
  iu_part = -1
  if (lttp%per_particle_output_file /= '') call ltt_per_particle_file_write1(orbit, lttp, lttp%ix_turn_start, &
                                             ix_bunch, ix_particle, size(beam%bunch), size(beam%bunch(ix_bunch)%particle))
else
  iu_part = lunget()
  if (lttp%phase_space_output_file == '') lttp%phase_space_output_file = 'single.dat'
  open(iu_part, file = lttp%phase_space_output_file, recl = 300)
  call ltt_write_params_header(lttp, ltt_com, iu_part)
  write (iu_part, '(2a)') '## Turn ix_ele |            x              px               y              py               z              pz', &
                                       '              pc             p0c            time   |       spin_x       spin_y       spin_z  | Element'
  write (iu_part, ltt_com%ps_fmt) lttp%ix_turn_start, ele_start%ix_ele, orbit%vec, (1.0_rp+orbit%vec(6))*orbit%p0c, &
                                                                      orbit%p0c, orbit%t, orbit%spin, trim(ele_start%name)
  if (lttp%custom_output_file /= '') call ltt_write_custom (lttp, ltt_com, lttp%ix_turn_start, orbit = orbit)
endif

!

do i_turn = lttp%ix_turn_start, lttp%ix_turn_stop-1
  orbit%ix_turn = i_turn

  select case (lttp%tracking_method)
  case ('BMAD')
    call ltt_track_bmad_single (lttp, ltt_com, ele_start, ele_start, orbit, i_turn, iu_part)

  case ('PTC')
    ele0 => ele_start
    do
      prb = orbit%vec
      prb%q%x = [1, 0, 0, 0]  ! Unit quaternion
      ele1 => pointer_to_next_ele(ele0)
      call track_probe (prb, ltt_com%ptc_state, fibre1 = pointer_to_fibre(ele0), fibre2 = pointer_to_fibre(ele1))
      orbit%vec = prb%x
      orbit%spin = quat_rotate(prb%q%x, orbit%spin)
      if (abs(orbit%vec(1)) > lttp%ptc_aperture(1) .or. abs(orbit%vec(3)) > lttp%ptc_aperture(2) .or. &
                                                     orbit_too_large(orbit) .or. prb%u) orbit%state = lost$

      if (orbit%state /= alive$) exit
      if (ele1%ix_ele == ele_start%ix_ele) exit
      ele0 => ele1
    enddo

  case ('MAP')
    call ltt_track_map(lttp, ltt_com, orbit, i_turn, iu_part)

  case default
    print '(a)', 'Unknown tracking_method: ' // lttp%tracking_method
    stop
  end select

  !

  if (lttp%simulation_mode == 'INDIVIDUAL' .and. lttp%per_particle_output_file /= '') then
    call ltt_per_particle_file_write1(orbit, lttp, i_turn+1, ix_bunch, ix_particle, size(beam%bunch), size(beam%bunch(ix_bunch)%particle))
  endif

  if (orbit%state /= alive$) then
    if (lttp%simulation_mode == 'INDIVIDUAL') then
      orb_in = orbit
      return
    endif
    print '(a, i0, 8a)', 'Particle lost at turn: ', i_turn
    exit
  endif

  if (lttp%simulation_mode == 'INDIVIDUAL') cycle

  do_write = (lttp%output_only_last_turns < 1 .or. lttp%n_turns - i_turn > lttp%output_only_last_turns)
  if (lttp%particle_output_every_n_turns > 0 .and. do_write) then
    if (modulo(i_turn+1, lttp%particle_output_every_n_turns) == 0) then
      write (iu_part, ltt_com%ps_fmt) i_turn+1, ele_start%ix_ele, orbit%vec, (1.0_rp+orbit%vec(6))*orbit%p0c, orbit%p0c, orbit%t, orbit%spin, trim(ele_start%name)
    endif
  endif

  if (lttp%averages_output_file /= '') then
    average = average + orbit%vec
    sigma = sigma + outer_product(orbit%vec, orbit%vec)
    n_sum = n_sum + 1
  endif

  if (lttp%custom_output_file /= '') call ltt_write_custom (lttp, ltt_com, i_turn+1, orbit = orbit)
enddo

if (lttp%simulation_mode == 'INDIVIDUAL') then
  orb_in = orbit
  return
endif

print '(2a)', 'Phase_space output file: ', trim(lttp%phase_space_output_file)
close(iu_part)

if (lttp%averages_output_file /= '') then
  call ltt_write_single_mode_sigma_file (lttp, n_sum, average, sigma)
endif

end subroutine ltt_run_single_mode

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_init_beam_distribution (lttp, ltt_com, beam)

type (ltt_params_struct) lttp
type (beam_struct), target :: beam
type (bunch_struct), pointer :: bunch
type (ltt_com_struct), target :: ltt_com
type (lat_struct), pointer :: lat
type (ele_struct), pointer :: ele_start
type (coord_struct), pointer :: p
type (random_state_struct) ran_state
type (branch_struct), pointer :: branch

real(rp) dt_branch, time, n_alive
integer ie, ib, n
logical err_flag
character(*), parameter :: r_name = 'ltt_init_beam_distribution'

! Important to apply the rampers to the starting element to get the correct ref momentum.
! But there is a chicken and egg problem here since the ramper time is determined from the beam.
! So do the init twice when rampers are on

lat => ltt_com%tracking_lat
branch => lat%branch(ltt_com%ix_branch)
call ltt_pointer_to_map_ends(lttp, lat, ele_start)

call run_timer('ABS', ltt_com%time_start)

!

ltt_com%beam_init_used%a_emit = real_garbage$
ltt_com%beam_init%ix_turn = lttp%ix_turn_start

if (lttp%ramping_on) then
  time = 0.5_rp * (ele_start%ref_time + ele_start%value(ref_time_start$)) + &
           lttp%ix_turn_start * (branch%ele(branch%n_ele_track)%ref_time - branch%ele(0)%ref_time) + &
           ltt_params_global%ramping_start_time
  call ltt_apply_rampers_to_slave (ltt_com, ele_start, time, err_flag)
  call ran_default_state(get_state = ran_state)
endif


call init_beam_distribution (ele_start, lat%param, ltt_com%beam_init, beam, err_flag, ltt_com%modes, &
                                                   ltt_com%beam_init_used, print_p0c_shift_warning = .false.)
bunch => beam%bunch(1)
if (err_flag) stop
call out_io (s_blank$, r_name, 'n_particle: ' // int_str(size(bunch%particle)))
ltt_com%n_particle = size(bunch%particle)

!

if (lttp%ramping_on) then
  n_alive = count(bunch%particle%state == alive$)
  time = sum(bunch%particle%t, bunch%particle%state == alive$) / n_alive + &
              0.5_rp * ele_start%value(delta_ref_time$) + ltt_params_global%ramping_start_time
  call ltt_apply_rampers_to_slave (ltt_com, ele_start, time, err_flag)
  call ran_default_state(set_state = ran_state)
  call init_beam_distribution (ele_start, lat%param, ltt_com%beam_init, beam, err_flag, ltt_com%modes, &
                                                     ltt_com%beam_init_used, print_p0c_shift_warning = .false.)
  bunch => beam%bunch(1)
endif

!

if (bmad_com%spin_tracking_on .and. all(ltt_com%beam_init%spin == 0) .and. all(bunch%particle%spin(1) == 0) .and. &
                            all(bunch%particle%spin(2) == 0) .and. all(bunch%particle%spin(3) == 0)) then
  ie = ele_start%ix_ele
  do ib = 1, size(beam%bunch)
    forall (n = 1:3) beam%bunch(ib)%particle%spin(n) = ltt_com%bmad_closed_orb(ie)%spin(n)
  enddo
endif

call ltt_setup_high_energy_space_charge(lttp, ltt_com, branch)

do ib = 1, size(beam%bunch)
  bunch => beam%bunch(ib)
  if (all(bunch%particle%charge == 0)) bunch%particle%charge = 1
  do n = 1, size(beam%bunch(ib)%particle)
    p => bunch%particle(n)
    if (lttp%add_closed_orbit_to_init_position) then
      select case (lttp%tracking_method)
      case ('MAP');    p%vec = p%vec + ltt_com%bmad_closed_orb(ele_start%ix_ele)%vec
      case ('BMAD');   p%vec = p%vec + ltt_com%bmad_closed_orb(ele_start%ix_ele)%vec
      case ('PTC');    p%vec = p%vec + ltt_com%ptc_closed_orb
      end select
    endif
  enddo
enddo

end subroutine ltt_init_beam_distribution

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_run_beam_mode (lttp, ltt_com, ix_start_turn, ix_end_turn, beam)

type (ltt_params_struct) lttp
type (beam_struct), target :: beam
type (bunch_struct), pointer :: bunch
type (branch_struct), pointer :: branch
type (ltt_com_struct), target :: ltt_com
type (lat_struct), pointer :: lat
type (coord_struct), allocatable :: orb(:)
type (coord_struct), pointer :: p
type (ele_struct), pointer :: ele_start
type (probe) prb

real(rp) time0, time1, time_now, dref_time, beta_new

integer ix_start_turn, ix_end_turn
integer n, n_print_dead_loss, i_turn, ie, n_part_tot
integer ip, ib, n_live_old, n_live, ix_branch

logical err_flag

character(16) prefix_str
character(200) file_name
character(*), parameter :: r_name = 'ltt_run_beam_mode'

! Init

ltt_com%ltt_tracking_happening_now = .true.    ! Start of main tracking
lat => ltt_com%tracking_lat
time0 = ltt_com%time_start
time1 = time0

ix_branch = ltt_com%ix_branch
branch => lat%branch(ix_branch)
call ltt_pointer_to_map_ends(lttp, lat, ele_start)
call ltt_write_per_particle_file_header(lttp, ltt_com, beam)

if (ltt_com%using_mpi) then
  prefix_str = 'Thread ' // int_str(ltt_com%mpi_rank) // ':'
else
  prefix_str = ''
endif

n_part_tot = 0
do ib = 1, size(beam%bunch)
  n_part_tot = n_part_tot + size(beam%bunch(ib)%particle)
enddo
n_print_dead_loss = max(1, nint(lttp%print_on_dead_loss * n_part_tot))
n_live_old = n_part_tot

! If using mpi then beam data will be written out by the master thread.

if (ix_start_turn == 0 .and. .not. ltt_com%using_mpi) then
  call ltt_write_particle_data (lttp, ltt_com, ix_start_turn, beam)
  call ltt_write_averages_data(lttp, ix_start_turn, beam)
  call ltt_write_custom (lttp, ltt_com, ix_start_turn, beam = beam)
endif

!

do i_turn = ix_start_turn, ix_end_turn-1
  do ib = 1, size(beam%bunch)
    bunch => beam%bunch(ib)
    bunch%ix_turn = i_turn   ! Turn index before tracking this turn

    select case (lttp%tracking_method)
    case ('MAP')
      do ip = 1, size(bunch%particle)
        p => bunch%particle(ip)
        if (p%state /= alive$) cycle
        call ltt_track_map (lttp, ltt_com, p)
      enddo

    case ('BMAD')
      call ltt_track_bmad_bunch (lttp, ltt_com, ele_start, bunch)

    case ('PTC')
      dref_time = branch%ele(branch%n_ele_track)%ref_time - branch%ele(0)%ref_time
      do ip = 1, size(bunch%particle)
        p => bunch%particle(ip)
        if (p%state /= alive$) cycle

        prb = p%vec
        prb%q%x = [1, 0, 0, 0]  ! Unit quaternion
        call track_probe (prb, ltt_com%ptc_state, fibre1 = ele_start%ptc_fibre)

        call convert_pc_to ((1 + prb%x(6)) * p%p0c, p%species, beta = beta_new)
        p%t = p%t + (p%vec(5) / p%beta - prb%x(5) / beta_new) / c_light + dref_time
        p%beta = beta_new

        p%vec = prb%x
        p%spin = quat_rotate(prb%q%x, p%spin)
        if (abs(p%vec(1)) > lttp%ptc_aperture(1) .or. abs(p%vec(3)) > lttp%ptc_aperture(2) .or. prb%u) p%state = lost$
      enddo

    case default
      print '(a)', 'Unknown tracking_method: ' // lttp%tracking_method
      stop
    end select

    bunch%ix_turn = i_turn + 1  ! Turn index after tracking this turn
    bunch%n_live = count(bunch%particle%state == alive$)

    if (i_turn+1 == lttp%ix_turn_record .and. lttp%ix_particle_record > 0) then
      bunch%particle(lttp%ix_particle_record)%ix_user = 1
    else
      bunch%particle%ix_user = -1
    endif
  enddo

  n_live = sum(beam%bunch%n_live)
  if (n_live_old - n_live >= n_print_dead_loss) then
    print '(2a, i0, a, i0)', trim(prefix_str), ' Cumulative number dead at turn ', i_turn+1, ': ', n_part_tot - n_live
    n_live_old = n_live
  endif

  if (n_live  == 0 .and. .not. ltt_com%using_mpi) then
    print '(2a)', trim(prefix_str), ' NO LIVE PARTICLES! STOPPING NOW.'
    exit
  endif

  if (n_part_tot - n_live >= nint(lttp%dead_cutoff * n_part_tot) .and. .not. ltt_com%using_mpi) then
    print '(2a)', trim(prefix_str), ' PARTICLE LOSS GREATER THAN SET BY DEAD_CUTOFF. STOPPING NOW.'
    exit
  endif

  call run_timer('ABS', time_now)
  if (time_now-time1 > lttp%timer_print_dtime .and. .not. ltt_com%using_mpi) then
    call out_io (s_blank$, r_name, trim(prefix_str) // ' Ellapsed time (min): ' // &
                                real_str((time_now-time0)/60, 10, 2) // ', At turn: ' // int_str(i_turn+1))
    time1 = time_now
  endif

  ! If using mpi then long_term_tracking_mpi will assemble the beam and call the write routines.

  if (.not. ltt_com%using_mpi) then
    call ltt_write_particle_data (lttp, ltt_com, i_turn+1, beam)
    call ltt_write_averages_data(lttp, i_turn+1, beam)
    call ltt_write_custom (lttp, ltt_com, i_turn+1, beam = beam)
    call ltt_write_beam_file(lttp, ltt_com, i_turn+1, beam)
  endif
end do

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
open(unit=22, file = "closed_orbit.dat")

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
call calc_z_tune (branch)

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

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine ltt_track1_radiation_center (orbit, ele, rad_damp, rad_fluct)
!
! Used for elements that have been split in half: This routine applies a kick to a particle 
! to account for radiation dampling and/or fluctuations.
!
! Input:
!   orbit     -- coord_struct: Particle at center of element before radiation applied.
!   ele       -- ele_struct: Radiation point
!   rad_damp  -- logical, optional: If present, override setting of bmad_com%radiation_damping_on.
!   rad_fluct -- logical, optional: If present, override setting of bmad_com%radiation_fluctuations_on.
!
! Output:
!   orbit     -- coord_struct: Particle position after radiation has been applied.
!-

subroutine ltt_track1_radiation_center (orbit, ele, rad_damp, rad_fluct)

use random_mod

implicit none

type (coord_struct) :: orbit
type (ele_struct), target :: ele
type (ele_struct), pointer :: ele1, ele2
type (lat_param_struct) :: param
type (rad_map_struct), pointer :: rad_mat

real(rp) int_gx, int_gy, this_ran, mc2, int_g2, int_g3, gxi, gyi, g2i, g3i, ran6(6)
real(rp) gamma_0, dE_p, fact_d, fact_f, q_charge2, p_spin, spin_norm(3), norm, rel_p
real(rp), parameter :: rad_fluct_const = 55.0_rp * classical_radius_factor * h_bar_planck * c_light / (24.0_rp * sqrt_3)
real(rp), parameter :: spin_const = 5.0_rp * sqrt_3 * classical_radius_factor * h_bar_planck * c_light / 16
real(rp), parameter :: damp_const = 2 * classical_radius_factor / 3
real(rp), parameter :: c1_spin = 2.0_rp / 9.0_rp, c2_spin = 8.0_rp / (5.0_rp * sqrt_3)

logical, optional :: rad_damp, rad_fluct
logical r_damp, r_fluct
character(*), parameter :: r_name = 'ltt_track1_radiation_center'

!

r_damp  = logic_option(bmad_com%radiation_damping_on, rad_damp)
r_fluct = logic_option(bmad_com%radiation_fluctuations_on, rad_fluct)
if (.not. r_damp .and. .not. r_fluct) return

rad_mat => ele%rad_map%rm1
if (r_damp) then
  orbit%vec = orbit%vec + bmad_com%synch_rad_scale * matmul(rad_mat%xfer_damp_mat, orbit%vec - rad_mat%ref_orb)
  if (.not. bmad_com%radiation_zero_average) orbit%vec = orbit%vec + bmad_com%synch_rad_scale * rad_mat%xfer_damp_vec
endif

if (r_fluct) then
  call ran_gauss (ran6)
  orbit%vec = orbit%vec + bmad_com%synch_rad_scale * matmul(rad_mat%stoc_mat, ran6)
endif

end subroutine ltt_track1_radiation_center 

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_setup_high_energy_space_charge(lttp, ltt_com, branch)

type (ltt_params_struct) lttp
type (branch_struct) branch
type (ltt_com_struct), target :: ltt_com
type (normal_modes_struct) modes
type (ele_struct), pointer :: ele_start, ele_stop

real(rp) n_particle, gamma
logical err

!

if (.not. bmad_com%high_energy_space_charge_on) return

if (lttp%tracking_method == 'MAP') then
  print '(a)', 'WARNING! Space effects are not present when using a map tracking!'
  return
endif

if (.not. lttp%rfcavity_on) then
  print '(a)', 'WARNING! RF is not on. Cannot calculate a longitudinal bunch length.'
  print '(a)', '      Therefore no space charge kick will be applied.'
  return
endif

gamma = branch%ele(0)%value(e_tot$) / mass_of(branch%ele(0)%ref_species)
modes = ltt_com%modes

call ltt_pointer_to_map_ends(lttp, ltt_com%tracking_lat, ele_start, ele_stop)
ltt_com%beam_init_used = set_emit_from_beam_init(ltt_com%beam_init, ele_start, ele_start%ref_species, ltt_com%modes, err)
if (err) stop

modes%a%emittance = ltt_com%beam_init_used%a_emit
modes%b%emittance = ltt_com%beam_init_used%b_emit

if (modes%a%emittance == 0 .or. modes%b%emittance == 0) then
  print *, 'WARNING! No a-mode or b-mode emittance set in beam_init structrue. Cannot compute high energy space charge kick.'
  print '(a)', '      Therefore no space charge kick will be applied.'
  return
endif

n_particle = abs(ltt_com%beam_init%bunch_charge / (e_charge * charge_of(ltt_com%bmad_closed_orb(0)%species)))

if (n_particle == 0) then
  print *, 'WARNING! beam_init%bunch_charge not set. Cannot compute high energy space charge kick.'
  print '(a)', '      Therefore no space charge kick will be applied.'
  return
endif
  

call setup_high_energy_space_charge_calc (.true., branch, n_particle, modes)

end subroutine ltt_setup_high_energy_space_charge

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_write_particle_data (lttp, ltt_com, i_turn, beam)

type (ltt_params_struct) lttp
type (ltt_com_struct), target :: ltt_com
type (beam_struct), target :: beam
type (bunch_struct), pointer :: bunch
type (coord_struct), pointer :: p

integer i_turn, ib, ip, nb
character(200) file

!

if (lttp%particle_output_every_n_turns > 0 .and. lttp%output_only_last_turns > 0 .and. &
                                        lttp%n_turns - i_turn >= lttp%output_only_last_turns) return

select case (lttp%particle_output_every_n_turns)
case (-1);    if (i_turn /= lttp%n_turns) return
case (0);     if (i_turn /= 0 .and. i_turn /= lttp%n_turns) return
case default; if (modulo(i_turn, lttp%particle_output_every_n_turns) /= 0) return
end select

!

nb = size(beam%bunch)
do ib = 1, nb
  bunch => beam%bunch(ib)

  if (lttp%phase_space_output_file /= '') call write_this_data (lttp, ltt_com, 'phase_space', &
            lttp%phase_space_output_file, bunch, i_turn, ib, nb, ltt_com%wrote_phase_space_file_header)

  if (lttp%action_angle_output_file /= '') call write_this_data (lttp, ltt_com, 'action_angle', &
            lttp%action_angle_output_file, bunch, i_turn, ib, nb, ltt_com%wrote_action_angle_file_header)

  if (lttp%per_particle_output_file /= '') then
    do ip = 1, size(bunch%particle)
      p => bunch%particle(ip)
      call ltt_per_particle_file_write1(p, lttp, i_turn, ib, ip, nb, size(bunch%particle))
    enddo
  endif
enddo

!--------------------------------------------------
contains

subroutine write_this_data (lttp, ltt_com, who, base_name, bunch, i_turn, ix_bunch, n_bunch, wrote_header)

type (ltt_params_struct) lttp
type (ltt_com_struct), target :: ltt_com
type (bunch_struct), target :: bunch
type (coord_struct), pointer :: p
type (bunch_params_struct) b_params
type (lat_struct), pointer :: lat
type (ele_struct), pointer :: ele_start

real(rp) n_inv_mat(6,6), jvec(6), jamp(3), jphase(3), m(6,6), vec0(6)
integer i_turn, ix_bunch, n_bunch, j, ix, iu, ip
logical wrote_header, error
character(*) who, base_name
character(200) file_name
character(40) str, fmt

!

if (base_name == '') return

iu = lunget()

if (lttp%n_turns > 0) then
  j = int(log10(real(lttp%n_turns, rp)) + 1 + 1d-10)
  str = int_str(i_turn, j)
else
  str = '000'
endif

if (n_bunch > 1) str = 'bunch' // int_str(ix_bunch) // '-' // str

ix = index(base_name, '#')
if (ix == 0) then
  file_name = trim(base_name) // str
else
  file_name = base_name(1:ix-1) // trim(str) // trim(base_name(ix+1:))
endif

! Currently str is always non-blank and a new file is created (no appending is done).
! Modifying the code to put everything in one file is a consideration but it is not clear if this is useful.

if (wrote_header .and. str == '') then
  open(iu, file = file_name, recl = 300, access = 'append')
else
  open(iu, file = file_name, recl = 300)
  call ltt_write_params_header(lttp, ltt_com, iu, ltt_com%n_particle, count(bunch%particle%state == alive$), bunch%ix_bunch)
  if (who == 'phase_space') then
    write (iu, '(a)')  '##     Ix     Turn |           x              px               y              py               z              pz              pc             p0c            time   |        spin_x       spin_y       spin_z    State'
  else
    write (iu, '(a)')  '##     Ix     Turn |          Ja         Angle_a              Jb         Angle_b              Jc         Angle_c   |     spin_x    spin_y    spin_z    State'
  endif
  wrote_header = .true.
endif

!

if (who == 'action_angle') then
  if (lttp%action_angle_calc_uses_1turn_matrix) then
    lat => ltt_com%tracking_lat
    call ltt_pointer_to_map_ends(lttp, lat, ele_start)

    call transfer_matrix_calc (lat, m, vec0, ele_start%ix_ele, ele_start%ix_ele, ele_start%ix_branch, .true.)
    if (all(m(1:5,6) == 0)) then
      print *, 'ltt%action_angle_calc_uses_1turn_matrix is set True but RF is off. No computations can be done. Stopping here.'
      stop
    endif

    call make_N(m, n_inv_mat, error)
    if (error) then
      print *, 'ltt%action_angle_calc_uses_1turn_matrix is set True but cannot get 1-turn matrix eigen vectors. Stopping here.'
      stop
    endif
    call mat_inverse (n_inv_mat, n_inv_mat)    

  else
    call calc_bunch_params (bunch, b_params, error, n_mat = n_inv_mat)
    if (error) then
      n_inv_mat = 0
    else
      call mat_inverse (n_inv_mat, n_inv_mat)
    endif
  endif
endif

do ip = 1, size(bunch%particle)
  p => bunch%particle(ip)
  ix = ip + ltt_com%mpi_ix0_particle
  if (lttp%only_live_particles_out .and. p%state /= alive$) cycle
  if (who == 'phase_space') then
    write (iu, '(i9, i9, 9es16.8, 3x, 3f13.9, 4x, a)')  ix, i_turn, p%vec, (1.0_rp+p%vec(6))*p%p0c, p%p0c, p%t, p%spin, trim(coord_state_name(p%state))
  else
    jvec = matmul(n_inv_mat, p%vec-b_params%centroid%vec)
    jamp = 0.5_rp * [jvec(1)**2 + jvec(2)**2, jvec(3)**2 + jvec(4)**2, jvec(5)**2 + jvec(6)**2]
    jphase = [atan2(jvec(2), jvec(1)), atan2(jvec(4), jvec(3)), atan2(jvec(6), jvec(5))]
    write (iu, '(i9, i9, 6es16.8, 3x, 3f13.9, 4x, a)')  ix, i_turn, jamp(1), jphase(1), jamp(2), jphase(2), jamp(3), jphase(3), p%spin, trim(coord_state_name(p%state))
  endif
enddo

close(iu)

end subroutine write_this_data

end subroutine ltt_write_particle_data

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_per_particle_file_write1 (p, lttp, i_turn, ix_bunch, ix_part, n_bunch, n_particle)

type (coord_struct) p
type (ltt_params_struct) lttp

integer i_turn, ix_bunch, ix_part, n_bunch, n_particle, iu
character(200) file

!

iu = lunget()
file = ltt_per_particle_file_name(lttp, ix_bunch, ix_part, n_bunch, n_particle)
open(iu, file = file, recl = 300, access = 'append')
write (iu, '(i9, 8es16.8, es16.8, 1x, 3f12.8, 4x, a)')  i_turn, p%vec, (1.0_rp+p%vec(6))*p%p0c, p%p0c, p%t, p%spin, trim(coord_state_name(p%state))
close(iu)

end subroutine ltt_per_particle_file_write1

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
! The optional argument "beam" is passed when doing beam tracking. Otherwise orbit is passed. Never both.

subroutine ltt_write_custom (lttp, ltt_com, i_turn, orbit, beam)

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
type (bunch_params_struct) bunch_params

integer i_turn
integer i, n, ib, iu, ix, ix1, ix2, is, multi, power, width, digits, n_loc, n_stack, n_bunch
integer n_particle, n_alive, ix_bunch

logical err

character(1000) line
character(100) err_str
character(40) fmt, str, ele_name, param
character(4) code

!

if (lttp%custom_output_file == '') return

if (lttp%particle_output_every_n_turns > 0 .and. lttp%output_only_last_turns > 0 .and. &
                                       lttp%n_turns - i_turn >= lttp%output_only_last_turns) return

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
n_bunch = 0
n_particle = 1

if (present(beam)) then
  n_bunch = size(beam%bunch)
  n_particle = size(beam%bunch(1)%particle)
  n_alive = count(beam%bunch(1)%particle%state == alive$)
else
  n_alive = count([orbit%state == alive$])
endif

if (i_turn == 0 .or. lttp%averages_output_every_n_turns == -1) then
  open(iu, file = lttp%custom_output_file, recl = 2000)
  call ltt_write_params_header(lttp, ltt_com, iu)
  line = '#'
  do i = 1, size(lttp%column)
    col => lttp%column(i)
    if (col%param == '') exit
    call parse_fortran_format (col%format, multi, power, code, width, digits)
    line = trim(line) // adjustr(col%header_str(1:width))
  enddo
  write (iu, '(a)') trim(line)

else
  open(iu, file = lttp%custom_output_file, access = 'append', recl = 2000)
endif

!

if (present(orbit)) then
  allocate (orb_b(1))
  orb_b(1) = orbit
else  ! beam is passed
  ix_bunch = -1
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

    ib = 1
    ix = index(st%name, '@')
    if (ix == 0) then
      str = st%name
    elseif (is_integer(st%name(1:ix-1), ib)) then
      str = st%name(ix+1:)
    endif

    if (present(orbit)) then
      orb = orbit
    else
      if (ib /= ix_bunch) call calc_bunch_params(beam%bunch(ib), bunch_params, err)
      ix_bunch = ib
      orb = bunch_params%centroid
    endif

    ! expression_string_to_stack will (unfortunately) upper case names. So we just live with it.
    select case (str)
    case ('N_TURN');      st%value = i_turn
    case ('N_PARTICLE');  st%value = n_particle
    case ('N_ALIVE');     st%value = n_alive
    case ('X');           st%value = orb%vec(1)
    case ('PX');          st%value = orb%vec(2)
    case ('Y');           st%value = orb%vec(3)
    case ('PY');          st%value = orb%vec(4)
    case ('Z');           st%value = orb%vec(5)
    case ('PZ');          st%value = orb%vec(6)
    case ('TIME');        st%value = orb%t
    case ('P0C');         st%value = (1.0_rp + orb%vec(6)) * orb%p0c
    case ('E_TOT');       st%value = (1.0_rp + orb%vec(6)) * orb%p0c / mass_of(orb%species)
    case ('SX');          st%value = orb%spin(1)
    case ('SY');          st%value = orb%spin(2)
    case ('SZ');          st%value = orb%spin(3)
    case default
      print '(2a)', 'Unknown parameter: ', trim(str)
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

close(iu)

end subroutine ltt_write_custom


!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_write_beam_file (lttp, ltt_com, i_turn, beam)

type (ltt_params_struct) lttp
type (ltt_com_struct), target :: ltt_com
type (beam_struct), target :: beam

integer i_turn, iu, ios, nt, ix, n
character(200) file

!

if (lttp%beam_output_file == '') return

nt = lttp%beam_output_every_n_turns
if (nt > 0 .and. lttp%output_only_last_turns > 0 .and. &
                                        lttp%n_turns - i_turn >= lttp%output_only_last_turns) return
if ((nt < 0 .or. mod(i_turn, nt) /= 0) .and. i_turn /= lttp%n_turns) return

file = lttp%beam_output_file
ix = str_last_in_set(file, '.')
n = int(log10(lttp%n_turns + 0.1_rp)) + 1
if (ix == 0) then
  file = trim(file) // '-' // int_str(i_turn, n)
else
  file = file(:ix-1) // '-' // int_str(i_turn, n) // file(ix:)
endif

call fullfilename(file, file)

call write_beam_file (file, beam)

!if (ltt_com%last_beam_output_file /= '') then
!  iu = lunget()
!  open(iu, file = ltt_com%last_beam_output_file, status = 'old')
!  close(iu, status = 'delete')
!endif

!ltt_com%last_beam_output_file = file

end subroutine ltt_write_beam_file

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_write_params_header(lttp, ltt_com, iu, n_particle, n_alive, ix_bunch)

type (ltt_params_struct) lttp
type (ltt_com_struct), target :: ltt_com

integer, optional :: n_particle, n_alive, ix_bunch
integer iu

!

write (iu,  '(3a)')      '# lattice                             = ', quote(lttp%lat_file)
write (iu,  '(3a)')      '# simulation_mode                     = ', quote(lttp%simulation_mode)
if (present(n_particle)) then
  write (iu,  '(a, i8)')   '# ix_bunch                            = ', ix_bunch
  write (iu,  '(a, i8)')   '# n_turns                             = ', lttp%n_turns
  write (iu,  '(a, i8)')   '# n_particle                          = ', n_particle
  write (iu,  '(a, i8)')   '# n_alive                             = ', n_alive
endif
write (iu,  '(a, l1)')   '# ramping_on                          = ', lttp%ramping_on
write (iu,  '(a, l1)')   '# ramp_update_each_particle           = ', lttp%ramp_update_each_particle
write (iu,  '(a, l1)')   '# ramp_particle_energy_without_rf     = ', lttp%ramp_particle_energy_without_rf
write (iu,  '(2a)')      '# ramping_start_time                  = ', real_str(lttp%ramping_start_time, 6)
write (iu,  '(a, i8)')   '# particle_output_every_n_turns       = ', lttp%particle_output_every_n_turns
write (iu,  '(a, i8)')   '# averages_output_every_n_turns       = ', lttp%averages_output_every_n_turns
write (iu,  '(a, i0)')   '# random_seed                         = ', lttp%random_seed
write (iu,  '(a, i0)')   '# random_seed_actual                  = ', ltt_com%random_seed_actual
write (iu,  '(a, l1)')   '# output_only_last_turns              = ', lttp%output_only_last_turns
write (iu,  '(a, l1)')   '# Radiation_Damping_on                = ', bmad_com%radiation_damping_on
write (iu,  '(a, l1)')   '# Radiation_Fluctuations_on           = ', bmad_com%radiation_fluctuations_on
write (iu,  '(a, l1)')   '# Spin_tracking_on                    = ', bmad_com%spin_tracking_on
write (iu,  '(a, l1)')   '# sr_wakes_on                         = ', bmad_com%sr_wakes_on
write (iu,  '(a, l1)')   '# RF_is_on                            = ', rf_is_on(ltt_com%tracking_lat%branch(ltt_com%ix_branch))
write (iu,  '(a, l1)')   '# action_angle_calc_uses_1turn_matrix = ', lttp%action_angle_calc_uses_1turn_matrix
if (bmad_com%sr_wakes_on) then
  write (iu, '(a, i0)')  '# Number_of_wake_elements             = ', size(ltt_com%ix_wake_ele)
endif
write (iu, '(3a)')       '# Map_file_prefix                     = ', quote(lttp%map_file_prefix)
if (lttp%tracking_method == 'MAP') then
  write (iu, '(a, i0)')  '# map_order                           = ', ltt_com%sec(1)%map%map_order
  write (iu, '(a, a)')   '# exclude_from_maps                   = ', lttp%exclude_from_maps
  write (iu, '(a, l1)')  '# split_bends_for_stochastic_rad      = ', lttp%split_bends_for_stochastic_rad
  write (iu, '(a, l1)')  '# symplectic_map_tracking             = ', lttp%symplectic_map_tracking
endif
write (iu, '(a)') '#'

end subroutine ltt_write_params_header

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_calc_bunch_data (lttp, ix_turn, bunch, bd)

type (ltt_params_struct) lttp
type (bunch_struct), target :: bunch
type (ltt_bunch_data_struct) :: bd

real(rp) ave, z, n_inv_mat(6,6)
real(rp) :: orb2_sum(6,6), orb3_sum(6), orb4_sum(6), orb_sum(6), core_emit_cutoff

integer ix_turn, i, j, k, n, ib, ic, ix, i_turn, this_turn, n2w
logical error

!

bd = ltt_bunch_data_struct()
call calc_bunch_params(bunch, bd%params, error, .true., n_inv_mat)
bd%n_live = bd%params%n_particle_live
if (bd%n_live < 10) then
  !! print *, 'Warning: Some parameters not computed since less than 10 particles are alive. From turn: ', ix_turn
  return
endif

orb_sum = 0;  orb2_sum = 0;  orb3_sum = 0;  orb4_sum = 0

do i = 1, 6
  orb_sum(i) = orb_sum(i)  + sum(bunch%particle%vec(i), bunch%particle%state == alive$)
  do j = i, 6
    orb2_sum(i,j) = orb2_sum(i,j) + sum(bunch%particle%vec(i) * bunch%particle%vec(j), bunch%particle%state == alive$) 
  enddo

  ave = sum(bunch%particle%vec(i), bunch%particle%state == alive$) / bd%n_live
  orb3_sum(i) = orb3_sum(i) + sum((bunch%particle%vec(i)-ave)**3, bunch%particle%state == alive$) 
  orb4_sum(i) = orb4_sum(i) + sum((bunch%particle%vec(i)-ave)**4, bunch%particle%state == alive$) 
enddo

bd%params%t = sum(bunch%particle%t, bunch%particle%state == alive$) / bd%n_live

!

k = 0
do i = 1, 6
do j = i, 6
  k = k + 1
  ! Test if the value of sigma(i,j) is significant. If not set to zero. 
  ! This is to avoid problems due to round-off errors. 
  ! One notible case is at the beginning of tracking if a mode has zero emittance.
  bd%sigma_vec(k) = orb2_sum(i,j) / bd%n_live - orb_sum(i) * orb_sum(j) / bd%n_live**2
  if (abs(bd%sigma_vec(k)) < 1e-15_rp * abs(orb_sum(i)*orb_sum(j)) / bd%n_live) bd%sigma_vec(k) = 0
  bd%params%sigma(i,j) = bd%sigma_vec(k)
  bd%params%sigma(j,i) = bd%sigma_vec(k)
  if (i == j) bd%sig1(i) = sqrt(max(0.0_rp, bd%sigma_vec(k)))
enddo
enddo

if (bd%sig1(1) /= 0) bd%skew(1) = orb3_sum(1) / (bd%n_live * bd%sig1(1)**3)
if (bd%sig1(1) /= 0) bd%kurt(1) = orb4_sum(1) / (bd%n_live * bd%sig1(1)**4) - 3.0_rp
if (bd%sig1(3) /= 0) bd%skew(2) = orb3_sum(3) / (bd%n_live * bd%sig1(3)**3)
if (bd%sig1(3) /= 0) bd%kurt(2) = orb4_sum(3) / (bd%n_live * bd%sig1(3)**4) - 3.0_rp
if (bd%sig1(5) /= 0) bd%skew(3) = orb3_sum(5) / (bd%n_live * bd%sig1(5)**3)
if (bd%sig1(5) /= 0) bd%kurt(3) = orb4_sum(5) / (bd%n_live * bd%sig1(5)**4) - 3.0_rp

! Calc core emit

bd%core_emit = 0

if (bd%params%a%emit == 0 .or. bd%params%b%emit == 0 .or. bd%params%c%emit == 0) then
  print '(a)', 'WARNING IN LTT_CALC_BUNCH_DATA: ZERO BEAM EMITTANCE MEANS I CANNOT DO CORE EMIT CALC.'
  return
endif

call mat_inverse (n_inv_mat, n_inv_mat)

do ic = 1, core_max$
  core_emit_cutoff = lttp%core_emit_cutoff(ic)
  if (core_emit_cutoff <= 0) exit
  if (bd%params%a%emit == 0 .or. bd%params%b%emit == 0 .or. bd%params%c%emit == 0) exit
  call this_core_calc(bunch, bd, core_emit_cutoff, n_inv_mat, bd%core_emit(ic,:))
enddo

!------------------------------------------------
contains

subroutine this_core_calc(bunch, bd, core_emit_cutoff, n_inv_mat0, core_emit)

type (bunch_struct), target :: bunch, core_bunch
type (ltt_bunch_data_struct) :: bd
type (bunch_params_struct) b_params

real(rp) core_emit_cutoff, core_emit(3)
real(rp) f, sig_cut, n_inv_mat0(6,6), n_inv_mat(6,6), cutoff

integer i, n, n_cut

logical error

!

n = bd%n_live
n_cut = int(core_emit_cutoff * n)
if (n_cut == 0) return

allocate (core_bunch%particle(n_cut))
cutoff = min(1.0_rp-1e-8_rp, core_emit_cutoff)
 
if (lttp%core_emit_combined_calc) then
  sig_cut = -log(inverse(beam_fraction, cutoff, 1e-12_rp, 1.0_rp, 1e-8_rp))
  f = cutoff / (1 - (1 + sig_cut + sig_cut**2/2 + sig_cut**3/6) * exp(-sig_cut))

  call core_bunch_construct(0, bunch, bd%params%centroid%vec, n_inv_mat0, n_cut, core_bunch, bd%params)

  call calc_bunch_params(core_bunch, b_params, error, n_mat = n_inv_mat);  if (error) return
  call mat_inverse (n_inv_mat, n_inv_mat)
  call core_bunch_construct(0, bunch, b_params%centroid%vec, n_inv_mat, n_cut, core_bunch, b_params)

  call calc_bunch_params(core_bunch, b_params, error);  if (error) return
  core_emit(1) = f * b_params%a%emit
  core_emit(2) = f * b_params%b%emit
  core_emit(3) = f * b_params%c%emit

else
  sig_cut = -log(1 - cutoff)
  f =  cutoff / (1 - (1+sig_cut)*exp(-sig_cut))

  do i = 1, 3
    call core_bunch_construct(i, bunch, bd%params%centroid%vec, n_inv_mat0, n_cut, core_bunch)

    call calc_bunch_params(core_bunch, b_params, error, n_mat = n_inv_mat); if (error) return

    call mat_inverse (n_inv_mat, n_inv_mat)
    call core_bunch_construct(i, bunch, b_params%centroid%vec, n_inv_mat, n_cut, core_bunch)

    call calc_bunch_params(core_bunch, b_params, error);  if (error) return
    select case (i)
    case (1); core_emit(1) = f * b_params%a%emit
    case (2); core_emit(2) = f * b_params%b%emit
    case (3); core_emit(3) = f * b_params%c%emit
    end select
  enddo
endif

end subroutine this_core_calc

!------------------------------------------------
! contains
! Fraction of the beam within region jx+jy+jz < jc where jx,jy,jz are normalized action coordinates.

function beam_fraction(exp_njc) result (fract)
real(rp) exp_njc, jc, fract
jc = -log(exp_njc)
fract = 1.0_rp - exp_njc * (1 + jc + 0.5_rp*jc**2)
end function beam_fraction

!------------------------------------------------
! contains

subroutine core_bunch_construct(ix_mode, bunch, center, n_inv_mat, n_cut, core_bunch, b_params)

type (bunch_struct), target :: bunch, core_bunch
type (bunch_params_struct), optional :: b_params
type (coord_struct), pointer :: p

real(rp) n_inv_mat(6,6), center(6), jvec(6)
real(rp), allocatable :: jamp(:)

integer ix_mode, n_cut
integer n, k, j, ip, ix
integer, allocatable :: indx(:)

!

n = size(bunch%particle)
allocate(jamp(n), indx(n))
jamp = 1e100_rp  ! Something large for dead particles

do ip = 1, size(bunch%particle)
  p => bunch%particle(ip)
  if (p%state /= alive$) cycle
  jvec = matmul(n_inv_mat, p%vec-center)

  if (ix_mode == 0) then
    jamp(ip) = (jvec(1)**2 + jvec(2)**2) / b_params%a%emit + (jvec(3)**2 + jvec(4)**2) / b_params%b%emit + &
                                                             (jvec(5)**2 + jvec(6)**2) / b_params%c%emit
  else
    ix = ix_mode * 2 - 1
    jamp(ip) = 0.5_rp * (jvec(ix)**2 + jvec(ix+1)**2)
  endif
enddo

call indexer(jamp, indx)

do ip = 1, n_cut
  j = indx(ip)
  if (bunch%particle(j)%state /= alive$) then
    print '(a, 2i8, f10.3)', 'ERROR IN CORE EMIT CALC. PLEASE REPORT!', &
                    count(bunch%particle%state == alive$), n_cut, core_emit_cutoff 
    stop
  endif
  core_bunch%particle(ip) = bunch%particle(j)
enddo

end subroutine core_bunch_construct

end subroutine ltt_calc_bunch_data

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_averages_file_name (template_name, suffix, ix_bunch, n_bunch, file_name)

character(*) template_name, suffix, file_name
integer ix_bunch, n_bunch   ! Index of bunch and number of bunches
integer ix

!

file_name = template_name
ix = index(file_name, '#')

if (n_bunch > 1) then
  if (ix == 0) then
    file_name = trim(file_name) // '.bunch' // int_str(ix_bunch) // '.#'
  else
    file_name = file_name(1:ix-1) // 'bunch' // int_str(ix_bunch) // '.#' // file_name(ix+1:) 
  endif
  ix = index(file_name, '#')
endif

if (suffix == '') return

if (ix == 0) then
  file_name = trim(file_name) // '.' // suffix
else
    file_name = file_name(1:ix-1) // trim(suffix) // file_name(ix+1:) 
endif

end subroutine ltt_averages_file_name

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_write_averages_data (lttp, ix_turn, beam, s_pos, ele)

type (ltt_params_struct) lttp
type (ltt_bunch_data_struct) bd
type (beam_struct), target :: beam
type (bunch_struct), pointer :: bunch
type (bunch_struct), target :: bunch0
type (ele_struct), optional :: ele

real(rp), optional :: s_pos
integer i, ix_turn, k, nb, iu1, iu2, iu3, ix, ib, n, n1
logical error

character(10) t_str
character(50) e_str
character(200) :: file, file_name
character(2000) :: line

! s_pos and ele are present in extraction tracking.

if (lttp%averages_output_file == '') return

if (lttp%averages_output_every_n_turns > 0 .and. lttp%output_only_last_turns > 0 .and. &
                                        lttp%n_turns - ix_turn >= lttp%output_only_last_turns) return

select case (lttp%averages_output_every_n_turns)
case (-1);    if (ix_turn /= lttp%n_turns) return
case (0);     if (modulo(ix_turn, lttp%n_turns) /= 0) return
case default; if (modulo(ix_turn, lttp%averages_output_every_n_turns) /= 0) return
end select

if (present(s_pos)) then
  write (t_str, '(f10.6)') s_pos
  e_str = '  ' // ele_full_name(ele, '!# @N')
  file = 'extraction-' // lttp%averages_output_file
else
  write (t_str, '(i10)') ix_turn
  e_str = ''
  file = lttp%averages_output_file
endif

nb = size(beam%bunch)
do ib = 0, nb

  if (ib == 0) then
    if (.not. lttp%output_combined_bunches .or. nb == 1) cycle
    n = 0
    do ix = 1, nb
      n = n + size(beam%bunch(ix)%particle)
    enddo
    allocate(bunch0%particle(n))
    bunch => bunch0
    n = 0
    do ix = 1, nb
      n1 = size(beam%bunch(ix)%particle)
      bunch%particle(n+1:n+n1) = beam%bunch(ix)%particle
      n = n + n1
    enddo
  else
    bunch => beam%bunch(ib)
    if (bunch%n_live == 0) exit
  endif

  call ltt_calc_bunch_data(lttp, ix_turn, bunch, bd)

  call ltt_averages_file_name(file, 'ave', ib, nb, file_name)
  iu1 = lunget(); open(iu1, file = file_name, recl = 400, access = 'append')
  call ltt_averages_file_name(file, 'sigma', ib, nb, file_name)
  iu2 = lunget(); open(iu2, file = file_name, recl = 400, access = 'append')
  call ltt_averages_file_name(file, 'emit', ib, nb, file_name)
  iu3 = lunget(); open(iu3, file = file_name, recl = 2000, access = 'append')

  !

  write (iu1, '(a, i9, es14.6, f14.9, 2x, 3f14.9, 2x, 8es14.6, 2x, 6es14.6, a)') t_str, nint(bd%n_live), &
                   bd%params%t, norm2(bd%params%centroid%spin), bd%params%centroid%spin, bd%params%centroid%vec, &
                   (1.0_rp+bd%params%centroid%vec(6))*bd%params%centroid%p0c, bd%params%centroid%p0c, bd%sig1, trim(e_str)

  write (iu2, '(a, i9, es14.6, 2x, 21es14.6, a)') t_str, nint(bd%n_live), bd%params%t, (bd%sigma_vec(k), k = 1, 21), trim(e_str)

  write (line, '(a, i9, es14.6, 2x, 3es14.6, 2x, 3es14.6, 2x, 3es14.6, a)') t_str, nint(bd%n_live), &
                             bd%params%t, bd%params%a%emit, bd%params%b%emit, bd%params%c%emit, bd%kurt, bd%skew, trim(e_str)

  do i = 1, core_max$
    if (lttp%core_emit_cutoff(i) <= 0) exit
    write (line, '(a, 2x, 3es14.6)') trim(line), bd%core_emit(i,:)
  enddo

  write (iu3, '(a)') trim(line)

  close(iu1)
  close(iu2)
  close(iu3)
enddo

end subroutine ltt_write_averages_data

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_write_single_mode_sigma_file (lttp, n_sum, average, sigma)

type (ltt_params_struct) lttp

real(rp) average(6), sigma(6,6)
integer i, n_sum

!

if (lttp%averages_output_file == '') return

open(1, file = lttp%averages_output_file)

if (n_sum == 0) then
  write (1, '(a)') '# NO DATA TO AVERAGE OVER!'
  return
endif

!

average = average / n_sum
sigma = sigma / n_sum - outer_product(average, average)

write (1, '(a)') '# Average:'
write (1, '(5x, 6es16.8)') average
write (1, *)
write (1, '(a)') '# Sigma:'
do i = 1, 6
  write (1, '(5x, 6es16.8)') sigma(i,:)
enddo

close(1)

print '(2a)', 'Sigma matrix data file: ', trim(lttp%averages_output_file)

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
type (ptc_rad_map_struct), pointer :: map

integer i, j, ix, n_sec, ib, ie, creation_hash
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
  if (ltt_com%mpi_rank == master_rank$) print '(a)', 'MAP FILE DOES NOT EXIST: ' // trim(map_file)
  return
endif

open(1, file = map_file, form = 'unformatted')
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
      print '(a)', 'ERROR READING MAP: ' // trim(map_file)
      return
    endif
    map => sec%map
  endif
enddo

close(1)

if (ltt_com%mpi_rank == master_rank$) then
  if (ltt_com%lat%creation_hash == creation_hash) then
    print '(2a)', 'Map read in from file: ', trim(map_file)
    print '(2a)', 'Lattice file used for map: ', trim(map%lattice_file)
  else
    print '(a)', 'NOTE: LATTICE HAS BEEN MODIFIED SINCE MAP FILE ' // trim(map_file) // ' WAS CREATED.'
    print '(a)', '      WILL MAKE A NEW MAP.'
  endif
endif

if (lttp%map_ascii_output_file /= '') then
  open(2, file = lttp%map_ascii_output_file)
  do i = 1, n_sec
    if (.not. allocated(ltt_com%sec(i)%map)) cycle
    map => ltt_com%sec(i)%map
    if (i > 1) write (2, *)
    write (2, '(a, i0, a)') 'Map: ', i, ' !-------------------------'
    write (2, *)
    ix = map%ix_ele_start
    write (2, '(a, i5, 2x, a30, f14.6)') 'Start Ele:', ix, branch%ele(ix)%name, branch%ele(ix)%s
    write (2, '(a, 6f12.6)') 'Orb start:', map%ref0
    ix = map%ix_ele_end
    write (2, '(a, i5, 2x, a30, f14.6)') 'End Ele:  ', ix, branch%ele(ix)%name, branch%ele(ix)%s
    write (2, '(a, 6f12.6)') 'Orb stop:', map%ref1
    write (2, *)
    call mat_type (map%nodamp_mat, 2, 'T-Matrix without Damping. Symplectic Error: ' // real_str(mat_symp_error(map%nodamp_mat), 6), '(4x, 6es16.8)')
    write (2, *)
    call mat_type (map%damp_mat, 2, 'D-Damping Matrix. Damp Factor: ' // real_str((1-determinant(map%damp_mat))/10, 6), '(4x, 6es16.8)')
    write (2, *)
    call mat_type (map%stoc_mat, 2, 'S-Radiation Matrix:', '(4x, 6es16.8)')
  enddo
  close(2)
  print *, 'Written: ', trim(lttp%map_ascii_output_file)
endif

if (ltt_com%lat%creation_hash == creation_hash) err_flag = .false.
return

9000 continue
close(1)
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

open(1, file = map_file, form = 'unformatted')
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

close(1)

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

!

if (lttp%map_file_prefix == '') then
  map_file = 'ltt'
else
  map_file = lttp%map_file_prefix
endif

! The "2" signifies version 2 of the map with radiation file storage format.

write (string, '(2a,2l1,i0,l1,i0,l1,es12.3,2l1)') trim(lttp%exclude_from_maps), trim(lttp%ele_start), &
              bmad_com%radiation_damping_on , lttp%split_bends_for_stochastic_rad, lttp%map_order, lttp%rfcavity_on, &
              ptc_com%max_fringe_order, ptc_com%old_integrator, ptc_com%vertical_kick, ptc_com%exact_misalign, ptc_com%exact_model
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
logical err, in_map_section

!

print '(a)', 'Creating map(s) for: ' // trim(ltt_com%lat%input_file_name)
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
  if (ele%ix_pointer == not_in_map$) then  ! Covers bends not in maps
    n_sec = n_sec + 2
  elseif (ele%key == sbend$) then
    n_sec = n_sec + 1
  endif
enddo

if (allocated(ltt_com%sec)) deallocate (ltt_com%sec)
allocate (ltt_com%sec(0:n_sec))
n_sec = 0
in_map_section = .false.

ltt_com%sec(0)%ele => ele_start
ltt_com%sec(0)%type = ele$

ele => ltt_com%sec(0)%ele
! Section setup
n_sec = 0
do i = 1, branch%n_ele_track+1
  ele0 => ele
  ele => pointer_to_next_ele(ele, skip_beginning = .false.)
  if (ele%key == marker$ .and. .not. in_map_section) cycle

  if (ele%ix_pointer == not_in_map$) then
    if (in_map_section) then
      n_sec = n_sec + 1
      allocate(ltt_com%sec(n_sec)%map)
      call ptc_setup_map_with_radiation (ltt_com%sec(n_sec)%map, ltt_com%sec(n_sec-1)%ele, ele0, &
                   lttp%map_order, bmad_com%radiation_damping_on, lttp%symplectic_map_tracking, err_flag = err)
      if (err) stop
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
                                             lttp%map_order, bmad_com%radiation_damping_on, lttp%symplectic_map_tracking)
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
call ltt_write_line(' Map setup time (min)' // real_str((time_now-time0)/60, 4), lttp, 0)

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
  eles(ie)%ele%ix_pointer = not_in_map$  ! Mark to exclude from any maps
enddo

!

ltt_com%tracking_lat = ltt_com%lat
branch => ltt_com%tracking_lat%branch(ltt_com%ix_branch)

if (lttp%split_bends_for_stochastic_rad) then
  call init_ele(marker, marker$)
  marker%name = 'RADIATION_POINT'
  marker%ix_pointer = not_in_map$

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

! Warning if beambeam element is included in a map

do ie = 1, branch%n_ele_track
  ele => branch%ele(ie)
  if (ele%ix_pointer == in_map$ .and. ele%key == beambeam$) then
    print '(a)', 'WARNING! Beambeam element is included in a map.'
    print '(a)', '          This is incaccurate at amplitudes larger than 1 sigma!'
  endif
enddo

! Mark slaves of lords that are not to be part of any map

do ie = ltt_com%tracking_lat%n_ele_track+1, ltt_com%tracking_lat%n_ele_max
  ele => ltt_com%tracking_lat%ele(ie)
  if (ele%ix_pointer == in_map$) cycle
  do is = 1, ele%n_slave
    slave => pointer_to_slave(ele, is)
    slave%ix_pointer = not_in_map$
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

subroutine ltt_track_bmad_single (lttp, ltt_com, ele_start, ele_stop, orbit, i_turn, iu_part)

type (ltt_params_struct) lttp
type (ltt_com_struct), target :: ltt_com
type (coord_struct) orbit
type (ele_struct), target :: ele_start, ele_stop
type (ele_struct), pointer :: ele

integer, optional :: i_turn, iu_part
integer i, it
logical err_flag, rad_fluct

!

if (lttp%split_bends_for_stochastic_rad) rad_fluct = set_parameter (bmad_com%radiation_fluctuations_on, .false.)
ele => ele_start
it = integer_option(0, i_turn)

do
  ele => pointer_to_next_ele(ele)
  if (ele%ix_ele == 0) it = it + 1
  if (ele%name == 'RADIATION_POINT') then
    call ltt_track1_radiation_center(orbit, ele, .false., rad_fluct)
  else
    call track1 (orbit, ele, ele%branch%param, orbit, err_flag = err_flag)
  endif

  if (orbit%state /= alive$) return

  if (integer_option(-1, iu_part) > 0 .and. lttp%particle_output_every_n_turns < 1) then
    write (iu_part, ltt_com%ps_fmt) it, ele%ix_ele, orbit%vec, (1.0_rp+orbit%vec(6))*orbit%p0c, orbit%p0c, orbit%t, orbit%spin, trim(ele%name)
  endif

  if (ele%ix_ele == ele_stop%ix_ele) exit
enddo

if (lttp%split_bends_for_stochastic_rad) bmad_com%radiation_fluctuations_on = rad_fluct

end subroutine ltt_track_bmad_single

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_track_bmad_bunch (lttp, ltt_com, ele_start, bunch)

type (ltt_params_struct) lttp
type (ltt_com_struct), target :: ltt_com
type (bunch_struct) bunch
type (ele_struct), target :: ele_start
type (ele_struct), pointer :: ele

integer i
logical err_flag, rad_fluct

!

if (lttp%split_bends_for_stochastic_rad) then
  rad_fluct = set_parameter (bmad_com%radiation_fluctuations_on, .false.)
  ele => ele_start

  do
    ele => pointer_to_next_ele(ele)
    if (ele%name == 'RADIATION_POINT') then
      do i = 1, size(bunch%particle)
        call ltt_track1_radiation_center(bunch%particle(i), ele, .false., rad_fluct)
      enddo
    else
      call track1_bunch (bunch, ele, err_flag)
    endif
    if (ele%ix_ele == ele_start%ix_ele) exit
  enddo

  bmad_com%radiation_fluctuations_on = rad_fluct

else
  call track_bunch (ele_start%branch%lat, bunch, ele_start, ele_start, err_flag)
endif

end subroutine ltt_track_bmad_bunch

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_track_map (lttp, ltt_com, orbit, i_turn, iu_part)

type (ltt_params_struct) lttp
type (ltt_com_struct), target :: ltt_com
type (coord_struct) orbit
type (ltt_section_struct), pointer :: sec
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele

integer, optional :: i_turn, iu_part
integer i
logical map_fluct_on, pt_fluct_on

!

map_fluct_on = (bmad_com%radiation_fluctuations_on .and. .not. lttp%split_bends_for_stochastic_rad) 
pt_fluct_on = bmad_com%radiation_fluctuations_on
branch => ltt_com%tracking_lat%branch(ltt_com%ix_branch)

do i = 1, ubound(ltt_com%sec, 1)
  sec => ltt_com%sec(i)
  if (.not. associated(sec%ele)) exit

  if (allocated(sec%map)) then
    call ptc_track_map_with_radiation (orbit, sec%map, rad_damp = bmad_com%radiation_damping_on, rad_fluct = map_fluct_on)
    if (abs(orbit%vec(1)) > lttp%ptc_aperture(1) .or. abs(orbit%vec(3)) > lttp%ptc_aperture(2)) orbit%state = lost$
    ele => branch%ele(sec%map%ix_ele_end)
  elseif (sec%ele%name == 'RADIATION_POINT') then
    call ltt_track1_radiation_center(orbit, sec%ele, .false., pt_fluct_on)
    ele => sec%ele
  else  ! EG: beambeam
    call track1 (orbit, sec%ele, branch%param, orbit)
    ele => sec%ele
  endif

  if (present(iu_part) .and. lttp%particle_output_every_n_turns < 1) then
    write (iu_part, ltt_com%ps_fmt) i_turn, ele%ix_ele, orbit%vec, (1.0_rp+orbit%vec(6))*orbit%p0c, orbit%p0c, orbit%t, orbit%spin, trim(ele%name)
  endif

  if (orbit%state /= alive$) return
enddo

end subroutine ltt_track_map

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ltt_print_mpi_info (lttp, ltt_com, line, do_print)

type (ltt_params_struct) lttp
type (ltt_com_struct), target :: ltt_com

real(rp) time_now
logical, optional :: do_print

character(*) line
character(20) time_str
character(200) str
character(*), parameter :: r_name = 'ltt_print_mpi_info'

!

if (.not. logic_option(lttp%debug, do_print)) return

call run_timer ('ABS', time_now)
call date_and_time_stamp (time_str)
write(str, '(a, f8.2, 2a, 2x, i0, 2a)') 'dTime:', (time_now-ltt_com%time_start)/60, &
                                        ' Now: ', time_str, ltt_com%mpi_rank, ': ', trim(line)

call out_io(s_blank$, r_name, str)
end subroutine ltt_print_mpi_info

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

! See documentation in long_term_tracking_mpi.f90 for why this routine is needed.

subroutine ltt_apply_rampers_to_slave (ltt_com, slave, time, err_flag)

type (ltt_com_struct), target :: ltt_com
type (ele_struct), target :: slave
type (ele_struct), pointer :: ramper
type (random_state_struct) ran_state_save
type (random_state_struct), pointer :: ran_state_ptr
type (lat_struct), pointer :: lat

real(rp) time
integer ir, iv
logical err_flag

! An element may be affected by ramping via a lord of the element.
! So it is only safe to not do any bookkeeping if the element has no rampers or lords.

if (slave%n_lord_ramper == 0 .and. slave%n_lord == 0) return

! Set rampers.

lat => ltt_com%tracking_lat
do ir = lat%n_ele_track+1, lat%n_ele_max
  ramper => lat%ele(ir)
  if (ramper%key /= ramper$) cycle
  do iv = 1, size(ramper%control%var)
    if (ramper%control%var(iv)%name /= 'TIME') cycle
    ramper%control%var(iv)%value = time
  enddo
enddo

! Swap out current "radiation" ran state for ramper ran state.
! This is for mpi running.

if (ltt_com_global%using_mpi) then
  ran_state_ptr => pointer_to_ran_state()
  ran_state_save = ran_state_ptr
  ran_state_ptr = ltt_com_global%ramper_ran_state
endif

! Apply rampers

call apply_rampers_to_slave (slave, err_flag)

! Swap out ramper ran state for radiation ran state

if (ltt_com_global%using_mpi) then
  ltt_com_global%ramper_ran_state = ran_state_ptr ! Update
  ran_state_ptr = ran_state_save
endif

end subroutine ltt_apply_rampers_to_slave

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine ltt_track1_bunch_hook (bunch, ele, err, centroid, direction, finished, bunch_track)
!
! This routine handles ramper element bookkeeping.
!
! This routine is called by the customized version of track1_bunch_hook that is linked with
! the long_term_tracking program.
!
! Input:
!   bunch_start   -- Bunch_struct: Starting bunch position.
!   ele          -- Ele_struct: Element to track through.
!   centroid(0:) -- coord_struct, optional: Approximate centroid orbit. Only needed if CSR is on.
!                     Hint: Calculate this before bunch tracking by tracking a single particle.
!   direction    -- integer, optional: +1 (default) -> Track forward, -1 -> Track backwards.
!
! Output:
!   bunch_end    -- bunch_struct: Ending bunch position.
!   err         -- Logical: Set true if there is an error. 
!                    EG: Too many particles lost for a CSR calc.
!   finished    -- logical: When set True, the standard track1_bunch code will not be called.
!-

subroutine ltt_track1_bunch_hook (bunch, ele, err, centroid, direction, finished, bunch_track)

type (bunch_struct), target :: bunch
type (ele_struct), target :: ele
type (ele_struct), pointer :: ele0
type (bunch_track_struct), optional :: bunch_track
type (coord_struct), optional :: centroid(0:)
type (coord_struct), pointer :: orb

real(rp) t, r

integer, optional :: direction
integer ip, ir, ie, n, iv, n_alive

logical err, finished

! Rampers are only applied to the element once per bunch. That is, it is assumed 
! that the ramper control function variation is negligible over the time scale of a bunch passage. 
! To evaluate multiple times in a bunch passage would, in general, be wrong if using ran() or ran_gauss().

err = .false.
finished = .false.
if (.not. ltt_params_global%ramping_on) return 
if (ltt_params_global%ramp_update_each_particle) return 

n = ltt_com_global%n_ramper_loc
if (n == 0) return

n_alive = count(bunch%particle%state == alive$)
if (n_alive == 0) return

t = sum(bunch%particle%t, bunch%particle%state == alive$) / n_alive + &
            0.5_rp * ele%value(delta_ref_time$) + ltt_params_global%ramping_start_time

call ltt_apply_rampers_to_slave (ltt_com_global, ele, t, err)

! The beginning element is never tracked through. If there is energy ramping and the user is writing out 
! p0c or E_tot from the beginning element, the user may be confused since these values will not change. 
! So adjust the beginning element's p0c and E_tot to keep users happy.

if (ele%ix_ele == 1) then
  ele0 => pointer_to_next_ele(ele, -1)
  ele0%value(p0c$) = ele%value(p0c_start$)
  ele0%value(E_tot$) = ele%value(E_tot_start$)
endif

! Adjust particle reference energy if needed.

if (bunch%particle(1)%p0c == ele%value(p0c_start$)) return

do ip = 1, size(bunch%particle)
  orb => bunch%particle(ip)
  if (orb%state /= alive$) cycle
  r = orb%p0c / ele%value(p0c_start$)
  orb%vec(2) = r * orb%vec(2)
  orb%vec(4) = r * orb%vec(4)
  ! Normally need to adjust pz (vec(6)) so that the particle energy (1 + orb%vec(6)) * orb%p0c is not changed.
  if (.not. ltt_params_global%ramp_particle_energy_without_rf) then
    orb%vec(6) = r * orb%vec(6) + (orb%p0c - ele%value(p0c_start$)) / ele%value(p0c_start$)
  endif
  orb%p0c = ele%value(p0c_start$)
enddo

end subroutine ltt_track1_bunch_hook

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine ltt_track1_preprocess (start_orb, ele, param, err_flag, finished, radiation_included, track)
!
! Routine for pre-processing at the start of the track1 routine.
!
! Also see:
!   track1_postprocess
!   track1_custom
!
! The radiation_included argument should be set to True if this routine (or a modified version of track1_custom)
! takes into account radiation damping and/or excitation. This will prevent track1 from calling track1_radiation.
! Note: If symp_lie_bmad is being called by this routine, symp_lie_bmad does take into account radiation effects.
! 
! General rule: Your code may NOT modify any argument that is not listed as an output agument below.
!
! Input:
!   start_orb  -- coord_struct: Starting position at the beginning of ele.
!   ele        -- ele_struct: Element.
!   param      -- lat_param_struct: Lattice parameters.
!
! Output:
!   start_orb   -- coord_struct: Modified starting position for track1 to use.
!   err_flag    -- logical: Set true if there is an error. False otherwise.
!   finished    -- logical: When set True, track1 will halt further processing and use the modified 
!                     start_orb as the final end position of the particle.
!   radiation_included
!               -- logical: Should be set True if radiation damping/excitation is included in the tracking.
!   track       -- track_struct, optional: Structure holding the track information if the 
!                    tracking method does tracking step-by-step.
!-

subroutine ltt_track1_preprocess (start_orb, ele, param, err_flag, finished, radiation_included, track)

type (coord_struct) :: start_orb
type (ele_struct), target :: ele
type (ele_struct), pointer :: ele0
type (lat_param_struct) :: param
type (track_struct), optional :: track

real(rp) r, t
integer n, iu
logical err_flag, finished, radiation_included, is_there

character(*), parameter :: r_name = 'ltt_track1_preprocess'

! This routine may be called by bmad_parser (via ele_compute_ref_energy_and_time) which is 
! before LTT ramping has been setup. To avoid problems, return if setup has not been done.

err_flag = .false.
if (.not. associated(ltt_com_global%tracking_lat%ele)) return

! Recording a particle track?

if (start_orb%ix_user > 0 .and. start_orb%state == alive$) then
  iu = lunget()
  if (ele%ix_ele <= 1) then
    open (iu, file = 'particle_track.dat')
  else
    open (iu, file = 'particle_track.dat', access = 'append')
  endif
  write (iu, '(i8, 2x, a20, 6es16.8)') ele%ix_ele, ele%name, start_orb%vec
  close (iu)
endif

if (.not. ltt_params_global%ramping_on) return
if (.not. ltt_params_global%ramp_update_each_particle .and. ltt_params_global%simulation_mode == 'BEAM') return 

! If bunch tracking, ramper bookkeeping is handled by track1_bunch_hook.

t = start_orb%t + 0.5_rp * ele%value(delta_ref_time$) + ltt_params_global%ramping_start_time
call ltt_apply_rampers_to_slave (ltt_com_global, ele, t, err_flag)

! The beginning element (with index 0) is never tracked through. If there is energy ramping and the user is 
! writing out p0c or E_tot from the beginning element, the user may be confused since these values will not change. 
! So adjust the beginning element's p0c and E_tot to keep users happy.

if (ele%ix_ele == 1) then
  ele0 => pointer_to_next_ele(ele, -1)
  ele0%value(p0c$) = ele%value(p0c_start$)
  ele0%value(E_tot$) = ele%value(E_tot_start$)
endif

! Adjust particle reference energy if needed.

if (start_orb%p0c == ele%value(p0c_start$)) return

r = start_orb%p0c / ele%value(p0c_start$)
start_orb%vec(2) = r * start_orb%vec(2)
start_orb%vec(4) = r * start_orb%vec(4)
! Normally need to adjust pz (vec(6)) so that the particle energy (1 + orb%vec(6)) * orb%p0c is not changed.
if (.not. ltt_params_global%ramp_particle_energy_without_rf) then
  start_orb%vec(6) = r * start_orb%vec(6) + (start_orb%p0c - ele%value(p0c_start$)) / ele%value(p0c_start$)
endif
start_orb%p0c = ele%value(p0c_start$)

end subroutine ltt_track1_preprocess

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
!-

subroutine ltt_init_coord (lttp, ltt_com, orbit, ele_start)

type (ltt_params_struct), target :: lttp
type (ltt_com_struct), target :: ltt_com
type (coord_struct) orbit
type (ele_struct) ele_start
type (branch_struct), pointer :: branch
type (random_state_struct) ran_state

real(rp) dt, time
integer n
logical err_flag

!

branch => ltt_com%tracking_lat%branch(ele_start%ix_branch)
dt = lttp%ix_turn_start * (branch%ele(branch%n_ele_track)%ref_time - branch%ele(0)%ref_time)

if (ltt_com%beam_init%use_particle_start) then
  call init_coord (orbit, ltt_com%lat%particle_start, ele_start, downstream_end$, branch%param%particle, t_offset = dt)
else
  call init_coord (orbit, ltt_com%beam_init%center, ele_start, downstream_end$, &
                                branch%param%particle, t_offset = dt, spin = ltt_com%beam_init%spin)
endif

if (lttp%ramping_on) then
  time = orbit%t + 0.5_rp * ele_start%value(delta_ref_time$) + ltt_params_global%ramping_start_time
  call ltt_apply_rampers_to_slave (ltt_com, ele_start, time, err_flag)
  call ran_default_state(set_state = ran_state)

  if (ltt_com%beam_init%use_particle_start) then
    call init_coord (orbit, ltt_com%lat%particle_start, ele_start, downstream_end$, branch%param%particle, t_offset = dt)
  else
    call init_coord (orbit, ltt_com%beam_init%center, ele_start, downstream_end$, &
                                  branch%param%particle, t_offset = dt, spin = ltt_com%beam_init%spin)
  endif
endif



if (lttp%add_closed_orbit_to_init_position) then
  if (lttp%simulation_mode == 'CHECK') then
    orbit%vec = orbit%vec + ltt_com%bmad_closed_orb(ele_start%ix_ele)%vec
  else
    select case (lttp%tracking_method)
    case ('MAP');     orbit%vec = orbit%vec + ltt_com%bmad_closed_orb(ele_start%ix_ele)%vec
    case ('BMAD');    orbit%vec = orbit%vec + ltt_com%bmad_closed_orb(ele_start%ix_ele)%vec
    case ('PTC');     orbit%vec = orbit%vec + ltt_com%ptc_closed_orb
    end select
  endif
endif

end subroutine ltt_init_coord

end module
