module tao_init_mod

use tao_interface
 
implicit none

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_init_global (init_file)
!
! Subroutine to initialize the tao global structures.
!
! Input:
!   init_file  -- Character(*): Tao initialization file.
!                  If blank, there is no file so just use the defaults.
!-

subroutine tao_init_global (init_file)

use random_mod
use opti_de_mod, only: opti_de_param
use input_mod

type (tao_global_struct) :: global

integer ios, iu, i, j, k, ix, num
integer n_data_max, n_var_max
integer n_d2_data_max, n_v1_var_max ! Deprecated variables
integer n, iostat

character(*) init_file
character(*), parameter :: r_name = 'tao_init_global'
character(200) file_name
character(40) name, universe

character(100) line

logical err, xxx

namelist / tao_params / global, bmad_com, csr_param, opti_de_param, &
          n_data_max, n_var_max, n_d2_data_max, n_v1_var_max
  
!-----------------------------------------------------------------------
! First time through capture the default global (could have been set via command line arg.)

global = s%global    ! establish defaults

call tao_hook_init_global (init_file, global)

! read global structure from tao_params namelist
! init_file == '' means there is no lattice file so just use the defaults.

if (init_file == '') then
  call end_bookkeeping()
  return
endif

call out_io (s_blank$, r_name, '*Init: Opening Init File: ' // init_file)
call tao_open_file (init_file, iu, file_name, s_blank$)
if (iu == 0) then
  call out_io (s_blank$, r_name, "Note: Cannot open initialization file for reading")
  call end_bookkeeping()
  return
endif

! Read tao_params

call out_io (s_blank$, r_name, 'Init: Reading tao_params namelist')
read (iu, nml = tao_params, iostat = ios)
if (ios > 0) then
  call out_io (s_error$, r_name, 'ERROR READING TAO_PARAMS NAMELIST.')
  rewind (iu)
  read (iu, nml = tao_params)  ! To give error message
endif
if (ios < 0) call out_io (s_blank$, r_name, 'Note: No tao_params namelist found')

! transfer global to s%global
s%global = global

close (iu)

call end_bookkeeping()

!-----------------------------------------------------------------------
contains

subroutine end_bookkeeping ()

! Tao does its own bookkeeping

bmad_com%auto_bookkeeper = .false.
s%com%valid_plot_who(1:5) = ['model ', 'base  ', 'ref   ', 'design', 'meas  ']

! Seed random number generator

call ran_seed_put (s%global%random_seed)
call ran_engine (s%global%random_engine)
call ran_gauss_converter (s%global%random_gauss_converter, s%global%random_sigma_cutoff)

call set_this_logical_command_arg (s%init%rf_on_arg, .false., s%global%rf_on)
call set_this_logical_command_arg (s%init%no_stopping_arg, .true., s%global%stop_on_error)
call set_this_logical_command_arg (s%init%noplot_arg, .true., s%global%plot_on)
call set_this_logical_command_arg (s%init%no_rad_int_arg, .true., s%global%rad_int_calc_on)
call set_this_logical_command_arg (s%init%disable_smooth_line_calc_arg, .false., s%global%disable_smooth_line_calc)
if (s%init%prompt_color_arg /= '')  s%global%prompt_color = s%init%prompt_color_arg
if (s%init%quiet_arg /= '')         s%global%quiet = s%init%quiet_arg

s%com%n_alias = 0
s%com%ix_key_bank = 0             ! For single mode.
s%com%use_saved_beam_in_tracking = .false.
if (.not. allocated(s%com%cmd_file)) allocate (s%com%cmd_file(0:0))

s%global%optimizer_allow_user_abort = (isatty(0) == 1)  ! Allow abort if input from tty (instead of a file).

call getenv ('ACC_PLOT_DISPLAY_TYPE', name)
if (name /= '') s%plot_page%plot_display_type = name

call readline_read_history(s%global%history_file, ios)

end subroutine end_bookkeeping

!-----------------------------------------------------------------------
! contains

subroutine set_this_logical_command_arg(cmd_arg, default, global_val)
character(*) cmd_arg
logical default, global_val
!
select case(cmd_arg)
case ('<present>');   global_val = .not. default
case ('<negated>');   global_val = default
end select
end subroutine set_this_logical_command_arg

end subroutine tao_init_global

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_init_beams (init_file)
!
! Subroutine to initialize beam stuff.
!
! Input:
!   init_file  -- Character(*): Tao initialization file.
!                  If blank, there is no file so just use the defaults.
!-

subroutine tao_init_beams (init_file)

use tao_input_struct

type (tao_universe_struct), pointer :: u
type (beam_init_struct) beam_init
type (tao_beam_branch_struct), pointer :: bb

integer i, k, iu, ios, ib, n_uni
integer n, iostat, ix_universe, ix_branch

character(*) init_file
character(40) :: r_name = 'tao_init_beams'
character(40) track_start, track_end, beam_track_start, beam_track_end
character(200) beam_saved_at, beam_dump_at, beam_dump_file
character(200) file_name, beam0_file, beam_track_data_file
character(200) beam_init_file_name           ! old style syntax
character(200) beam_position0_file           ! old style syntax
character(60), target :: save_beam_at(100)   ! old style syntax

logical err

namelist / tao_beam_init / ix_universe, ix_branch, beam0_file, beam_init, beam_init_file_name, &
            beam_saved_at, track_start, track_end, beam_track_start, beam_track_end, beam_position0_file, &
            beam_track_data_file, beam_dump_at, beam_dump_file
         
!-----------------------------------------------------------------------
! Init Beams

call tao_hook_init_beam ()

do i = lbound(s%u, 1), ubound(s%u, 1)
  bb => s%u(i)%model_branch(0)%beam
  bb%beam_init = beam_init_struct()
  bb%beam_init%file_name = ''
  bb%beam_init%position_file = s%init%beam_init_position_file_arg
  bb%track_start    = ''
  bb%track_end      = ''
  bb%ix_track_start = not_set$
  bb%ix_track_end   = not_set$
  s%u(i)%beam%track_data_file = s%init%beam_track_data_file_arg
  s%u(i)%beam%track_beam_in_universe = .false.
enddo

if (.not. s%com%init_beam .or. init_file == '') return

!

call out_io (s_blank$, r_name, '*Init: Opening Beam File: ' // init_file)
call tao_open_file (init_file, iu, file_name, s_fatal$)
if (iu == 0) then
  call out_io (s_fatal$, r_name, 'CANNOT OPEN BEAM INIT FILE: ' // init_file)
  return
endif

do
  ! defaults
  ix_universe = -1
  ix_branch = 0
  beam_init = beam_init_struct()
  beam_saved_at = ''
  beam_dump_at = ''
  beam_dump_file = ''
  save_beam_at  = ''
  track_start = ''        ! Old style
  track_end = ''          ! Old style
  beam_track_start = ''
  beam_track_end = ''
  beam_init_file_name = ''
  beam0_file = ''
  beam_position0_file = ''
  beam_track_data_file = ''

  ! Read beam parameters

  read (iu, nml = tao_beam_init, iostat = ios)
  if (ios > 0) then
    call out_io (s_abort$, r_name, 'INIT: TAO_BEAM_INIT NAMELIST READ ERROR!')
    rewind (iu)
    do
      read (iu, nml = tao_beam_init)  ! generate an error message
    enddo
  endif

  if (beam0_file /= '') then
    beam_init%position_file = beam0_file
    call out_io (s_important$, r_name, &
          'Note: Parameter beam0_file in the tao_beam_init structure has been replaced by beam_init%position_file.', &
          'PLEASE MODIFY YOUR INPUT FILE. This is just a warning. Tao will run normally...')
  endif

  if (beam_position0_file /= '') then
    beam_init%position_file = beam_position0_file
    call out_io (s_important$, r_name, &
          'Note: Parameter beam_position0_file in the tao_beam_init structure has been replaced by beam_init%position_file.', &
          'PLEASE MODIFY YOUR INPUT FILE. This is just a warning. Tao will run normally...')
  endif

  if (beam_init_file_name /= '') then
    beam_init%position_file = beam_init_file_name
    call out_io (s_important$, r_name, &
          'Note: Parameter beam_init_file_name in the tao_beam_init structure has been replaced by beam_init%position_file.', &
          'PLEASE MODIFY YOUR INPUT FILE. This is just a warning. Tao will run normally...')
  endif

  if (beam_init%file_name /= '') then
    beam_init%position_file = beam_init%position_file
    call out_io (s_important$, r_name, &
          'Note: Parameter beam_init%file_name in the tao_beam_init structure has been replaced by beam_init%position_file.', &
          'PLEASE MODIFY YOUR INPUT FILE. This is just a warning. Tao will run normally...')
  endif

  if (s%init%beam_init_position_file_arg /= '') beam_init%position_file = s%init%beam_init_position_file_arg

  if (beam_init%sig_e /= 0 .and. beam_init%sig_pz /= 0) then   ! sig_e is superceeded by sig_pz
    beam_init%sig_pz = beam_init%sig_e
    beam_init%sig_e = 0
  endif

  if (s%init%beam_track_data_file_arg /= '') beam_track_data_file = s%init%beam_track_data_file_arg  ! From the command line
  if (track_start /= '') beam_track_start = track_start   ! For backwards compatibility
  if (track_end /= '')   beam_track_end   = track_end     ! For backwards compatibility

  ! transfer info from old style save_beam_at(:) to beam_saved_at

  do i = 1, size(save_beam_at)
    if (save_beam_at(i) == '') cycle
    call out_io (s_error$, r_name, 'OLD STYLE "save_beam_at" REPLACED BY "beam_saved_at".', &
                                   'PLEASE MODIFY ACCORDINGLY. TAO WILL RUN NORMALLY FOR NOW...')
    beam_saved_at = trim(beam_saved_at) // ', ' // trim(save_beam_at(i))
  enddo

  if (ios < 0 .and. ix_universe == -1) exit  ! Exit on end-of-file and no namelist read

  ! init

  call out_io (s_blank$, r_name, 'Init: Read tao_beam_init namelist for universe \i3\ ', ix_universe)
  if (ix_universe == -1) then
    do i = lbound(s%u, 1), ubound(s%u, 1)
      s%u(i)%beam%dump_file = beam_dump_file
      call init_beam(s%u(i), beam_init)
    enddo
  else
    if (ix_universe < lbound(s%u, 1) .or. ix_universe > ubound(s%u, 1)) then
      call out_io (s_error$, r_name, 'BAD IX_UNIVERSE IN TAO_BEAM_INIT NAMELIST: \i0\ ', ix_universe)
      return
    endif
    s%u(ix_universe)%beam%dump_file = beam_dump_file
    call init_beam(s%u(ix_universe), beam_init)
  endif

enddo

close (iu)

!----------------------------------------------------------------
!----------------------------------------------------------------
contains

! Initialize the beams. Determine which element to track beam to

subroutine init_beam (u, beam_init)

use beam_file_io

type (tao_universe_struct), target :: u
type (beam_init_struct) beam_init
type (ele_pointer_struct), allocatable, save, target :: eles(:)
type (ele_struct), pointer :: ele
type (branch_struct), pointer :: branch

real(rp) v(6), bunch_charge, gamma
integer i, j, ix, iu, n_part, ix_class, n_bunch, n_particle, n_loc
character(60) at, class, ele_name, line

! Set tracking start/stop

u%beam%track_beam_in_universe = .true.
bb => u%model_branch(0)%beam
bb%track_start = beam_track_start
bb%track_end   = beam_track_end

if (beam_track_start /= '') then
  call lat_ele_locator (beam_track_start, u%design%lat, eles, n_loc, err)
  if (err .or. n_loc == 0) then
    call out_io (s_error$, r_name, 'BEAM_TRACK_START ELEMENT NOT FOUND: ' // beam_track_start, &
                                   'WILL NOT TRACK A BEAM.')
    s%global%track_type = 'single'
    return
  endif
  if (n_loc > 1) then
    call out_io (s_error$, r_name, 'MULTIPLE BEAM_TRACK_START ELEMENTS FOUND: ' // beam_track_start, &
                                   'WILL NOT TRACK A BEAM.')
    s%global%track_type = 'single'
    return
  endif
  bb%ix_track_start = eles(1)%ele%ix_ele
endif

if (beam_track_end /= '') then
  call lat_ele_locator (beam_track_end, u%design%lat, eles, n_loc, err)
  if (err .or. n_loc == 0) then
    call out_io (s_error$, r_name, 'BEAM_TRACK_END ELEMENT NOT FOUND: ' // beam_track_end, &
                                   'WILL NOT TRACK A BEAM.')
    s%global%track_type = 'single'
    return
  endif
  if (n_loc > 1) then
    call out_io (s_error$, r_name, 'MULTIPLE BEAM_TRACK_END ELEMENTS FOUND: ' // beam_track_end, &
                                   'WILL NOT TRACK A BEAM.')
    s%global%track_type = 'single'
    return
  endif
  bb%ix_track_end = eles(1)%ele%ix_ele
endif

bb%beam_init = beam_init
u%beam%track_data_file = beam_track_data_file

! Find where to save the beam at.
! Note: Beam will automatically be saved at fork elements and at the ends of the beam tracking.

do i = 0, ubound(u%model%lat%branch, 1)
  branch => u%design%lat%branch(i)
  u%model_branch(i)%ele%save_beam_internally = .false.
  u%model_branch(i)%ele%save_beam_to_file = .false.
enddo

u%beam%saved_at = beam_saved_at
if (beam_saved_at /= '') then
  call tao_locate_elements (beam_saved_at, u%ix_uni, eles, err, ignore_blank = .false.)
  if (err) then
    call out_io (s_error$, r_name, 'BAD "beam_saved_at" ELEMENT: ' // beam_saved_at)
  else
    do k = 1, size(eles)
      ele => eles(k)%ele
      u%model_branch(ele%ix_branch)%ele(ele%ix_ele)%save_beam_internally = .true.
    enddo
  endif
endif

u%beam%dump_at = beam_dump_at
if (beam_dump_at /= '') then
  call tao_locate_elements (beam_dump_at, u%ix_uni, eles, err, ignore_blank = .false.)
  if (err) then
    call out_io (s_error$, r_name, 'BAD "beam_dump_at" ELEMENT: ' // beam_dump_at)
  else
    do k = 1, size(eles)
      ele => eles(k)%ele
      u%model_branch(ele%ix_branch)%ele(ele%ix_ele)%save_beam_to_file = .true.
    enddo
  endif
endif

! If beam_track_data_file is set, read in the beam tracking data.

if (u%beam%track_data_file /= '') then
  call out_io (s_fatal$, r_name, 'beam_track_data_file not yet implemented. Please contact David Sagan...')
  stop
endif

if (allocated(eles)) deallocate (eles)

end subroutine init_beam

end subroutine tao_init_beams


!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_init_dynamic_aperture (init_file)
!
! Routine to initalize dynamic aperture simulations.
!
! Input:
!   init_file   -- character(*): File setting dynamic_aperture parameters.
!-

subroutine tao_init_dynamic_aperture(init_file)

use tao_input_struct

type (tao_dynamic_aperture_input) :: da_init(100)
type (tao_universe_struct), pointer :: u

integer :: ios, iu, i, j, n_pz

character(*) init_file
character(200) file_name
character(40) :: r_name = 'tao_init_dynamic_aperture'

namelist / tao_dynamic_aperture / da_init

!

if (init_file == '') return

call tao_open_file (init_file, iu, file_name, s_blank$)
if (iu == 0) then
  call out_io (s_blank$, r_name, "Note: Cannot open init file for tao_dynamic_aperture namelist read")
  return
endif

! Read tao_dynamic_aperture
call out_io (s_blank$, r_name, 'Init: Reading tao_dynamic_aperture namelist')
read (iu, nml = tao_dynamic_aperture, iostat = ios)
if (ios > 0) then
  call out_io (s_error$, r_name, 'ERROR READING TAO_DYNAMIC_APERTURE NAMELIST.')
  rewind (iu)
  read (iu, nml = tao_dynamic_aperture)  ! To give error message
endif

close(iu)

if (ios < 0) then
  call out_io (s_blank$, r_name, 'Note: No tao_dynamic_aperture namelist found')
  return
endif

do i = lbound(s%u, 1), ubound(s%u, 1)
  ! Count the list of pz
  do n_pz = 0, size(da_init(i)%pz)-1
    if (da_init(i)%pz(n_pz+1) == real_garbage$) exit
  enddo

  ! Default if no pz set
  if (n_pz == 0) then
    n_pz = 1
    da_init(i)%pz(1) = 0
  endif

  ! Set
  u => s%u(i)
  allocate(u%dynamic_aperture%scan(n_pz))
  allocate(u%dynamic_aperture%pz(n_pz))
  u%dynamic_aperture%param%n_turn    = da_init(i)%n_turn
  u%dynamic_aperture%param%x_init    = da_init(i)%x_init
  u%dynamic_aperture%param%y_init    = da_init(i)%y_init
  u%dynamic_aperture%param%accuracy  = da_init(i)%accuracy
  u%dynamic_aperture%param%min_angle = da_init(i)%min_angle
  u%dynamic_aperture%param%max_angle = da_init(i)%max_angle
  u%dynamic_aperture%param%n_angle   = da_init(i)%n_angle
  u%dynamic_aperture%pz(1:n_pz)      = da_init(i)%pz(1:n_pz)  
enddo

end subroutine tao_init_dynamic_aperture

end module
