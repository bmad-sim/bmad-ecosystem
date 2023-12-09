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

namelist / tao_params / global, bmad_com, space_charge_com, opti_de_param, &
          n_data_max, n_var_max, n_d2_data_max, n_v1_var_max
  
!-----------------------------------------------------------------------
! First time through capture the default global (could have been set via command line arg.)

global = s%global    ! establish defaults

if (associated(tao_hook_init_global_ptr)) call tao_hook_init_global_ptr(init_file, global)

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

! Not using auto bookkeeping will speed up calculations.

bmad_com%auto_bookkeeper = .false.
s%com%valid_plot_who(1:5) = ['model ', 'base  ', 'ref   ', 'design', 'meas  ']

call set_this_logical_command_arg (s%init%disable_smooth_line_calc_arg, .false., s%global%disable_smooth_line_calc)
call set_this_logical_command_arg (s%init%no_stopping_arg, .true., s%global%stop_on_error)
call set_this_logical_command_arg (s%init%noplot_arg, .true., s%global%plot_on)
call set_this_logical_command_arg (s%init%no_rad_int_arg, .true., s%global%rad_int_calc_on)
call set_this_logical_command_arg (s%init%rf_on_arg, .false., s%global%rf_on)
call set_this_logical_command_arg (s%init%symbol_import_arg, .false., s%global%symbol_import)

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
type (branch_struct), pointer :: branch

real(rp) comb_ds_save, comb_max_ds_save

integer i, k, iu, ios, ib, n_uni, ib0, ie0
integer n, iostat, ix_universe

character(*) init_file
character(40) track_start, track_end, beam_track_start, beam_track_end
character(200) file_name
character(200) beam0_file, beam_init_file_name, beam_position0_file    ! Very old style syntax
character(200) beam_saved_at, beam_dump_at, beam_dump_file
character(200) saved_at, dump_at, dump_file
character(*), parameter :: r_name = 'tao_init_beams'

logical err, always_reinit

namelist / tao_beam_init / ix_universe, beam_init, always_reinit, &
            beam0_file, beam_init_file_name, beam_position0_file, &
            beam_track_start, beam_track_end, beam_saved_at, beam_dump_at, beam_dump_file, &
            track_start, track_end, saved_at, dump_at, dump_file, comb_ds_save, comb_max_ds_save

!-----------------------------------------------------------------------
! Init Beams
! Some init that will be needed with lat sigma tracking

if (associated(tao_hook_init_beam_ptr)) call tao_hook_init_beam_ptr()

do i = lbound(s%u, 1), ubound(s%u, 1)
  u => s%u(i)
  do ib = 0, ubound(u%model%lat%branch, 1)
    branch => u%model%lat%branch(ib)
    bb => u%model_branch(ib)%beam

    u%model_branch(ib)%ele%save_beam_internally = .false.
    u%model_branch(ib)%ele%save_beam_to_file = .false.

    if (branch%ix_from_branch < 0) then  ! Root branch
      bb%beam_init      = beam_init_struct()
      bb%beam_init_used = beam_init_struct()
      bb%beam_init_used%a_emit = real_garbage$  ! Tag to see when structure is set.
      bb%beam_init%position_file = s%init%beam_init_position_file_arg
      bb%track_start    = ''
      bb%track_end      = ''
      bb%ix_track_start = 0
      bb%ix_track_end   = branch%n_ele_track
    else                                ! Branch injected into
      bb%ix_track_start = branch%ix_to_ele 

      if (branch%param%geometry == open$) then
        bb%ix_track_end = branch%n_ele_track
      else
        bb%ix_track_end = bb%ix_track_start - 1
        if (bb%ix_track_end == -1) bb%ix_track_end = branch%n_ele_track
      endif

      bb%track_start = branch%ele(bb%ix_track_start)%name
      bb%track_end   = branch%ele(bb%ix_track_end)%name
    endif
  enddo

  u%beam%track_beam_in_universe = .false.
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
  always_reinit = .false.
  beam_init = beam_init_struct()

  beam0_file = ''             ! Very old style
  beam_position0_file = ''    ! Very old style
  beam_init_file_name = ''    ! Very old style

  beam_saved_at = ''          ! Old style
  beam_dump_at = ''           ! Old style
  beam_dump_file = ''         ! Old style
  beam_track_start = ''       ! Old style
  beam_track_end = ''         ! Old style

  saved_at = ''
  dump_at = ''
  dump_file = ''
  track_start = ''
  track_end = ''
  comb_ds_save = -1
  comb_max_ds_save = -1

  ! Read beam parameters

  read (iu, nml = tao_beam_init, iostat = ios)
  if (ios > 0) then
    call out_io (s_abort$, r_name, 'INIT: TAO_BEAM_INIT NAMELIST READ ERROR!')
    rewind (iu)
    do
      read (iu, nml = tao_beam_init)  ! generate an error message
    enddo
  endif

  if (ios < 0 .and. ix_universe == -1) exit  ! Exit on end-of-file and no namelist read

  ! Convert old style to current style

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

  if (beam_saved_at /= '')        saved_at        = beam_saved_at
  if (beam_dump_at /= '')         dump_at         = beam_dump_at
  if (beam_dump_file /= '')       dump_file       = beam_dump_file
  if (beam_track_start /= '')     track_start     = beam_track_start
  if (beam_track_end /= '')       track_end       = beam_track_end

  !

  if (s%init%beam_init_position_file_arg /= '') beam_init%position_file = s%init%beam_init_position_file_arg

  if (beam_init%sig_e /= 0 .and. beam_init%sig_pz /= 0) then   ! sig_e is superceeded by sig_pz
    beam_init%sig_pz = beam_init%sig_e
    beam_init%sig_e = 0
  endif

  ! init

  call out_io (s_blank$, r_name, 'Init: Read tao_beam_init namelist for universe \i3\ ', ix_universe)
  if (ix_universe == -1) then
    do i = lbound(s%u, 1), ubound(s%u, 1)
      s%u(i)%beam = tao_beam_uni_struct(saved_at, dump_file, dump_at, .true., always_reinit)
      call tao_init_beam_in_universe(s%u(i), beam_init, track_start, track_end, comb_ds_save, comb_max_ds_save)
    enddo
  else
    if (ix_universe < lbound(s%u, 1) .or. ix_universe > ubound(s%u, 1)) then
      call out_io (s_error$, r_name, 'BAD IX_UNIVERSE IN TAO_BEAM_INIT NAMELIST: \i0\ ', ix_universe)
      return
    endif
    s%u(ix_universe)%beam = tao_beam_uni_struct(saved_at, dump_file, dump_at, .true., always_reinit)
    call tao_init_beam_in_universe(s%u(ix_universe), beam_init, track_start, track_end, comb_ds_save, comb_max_ds_save)
  endif

enddo

close (iu)

end subroutine tao_init_beams

!----------------------------------------------------------------
!----------------------------------------------------------------
!----------------------------------------------------------------

! Initialize the beams. Determine which element to track beam to

subroutine tao_init_beam_in_universe (u, beam_init, track_start, track_end, comb_ds_save, comb_max_ds_save)

type (tao_universe_struct), target :: u
type (beam_init_struct) beam_init
type (ele_struct), pointer :: ele
type (ele_pointer_struct), allocatable, target :: eles(:)
type (branch_struct), pointer :: branch
type (tao_beam_branch_struct), pointer :: bb
type (tao_lattice_branch_struct), pointer :: tao_branch

real(rp) comb_ds_save, comb_max_ds_save
integer k, n_loc

logical always_reinit, err

character(*) track_start, track_end
character(*), parameter :: r_name = 'tao_init_beam_in_universe'

! Set tracking start

u%beam%track_beam_in_universe = .true.

if (track_start == '') then
  ele => u%model%lat%branch(0)%ele(0)
else
  ele => tao_beam_track_endpoint (track_start, u%model%lat, '', 'TRACK_START')
  if (.not. associated(ele)) return
endif

bb => u%model_branch(ele%ix_branch)%beam
bb%ix_track_start = ele%ix_ele
bb%beam_init = beam_init
bb%track_start = track_start

tao_branch => u%model%tao_branch(ele%ix_branch) 
tao_branch%comb_ds_save = comb_ds_save
!! tao_branch%max_ds_save = comb_max_ds_save  ! Note: Not yet implemented

! Tracking stop

bb%track_end = track_end
branch => u%model%lat%branch(ele%ix_branch)

if (track_end == '') then
  if (branch%param%geometry == open$) then
    bb%ix_track_end = branch%n_ele_track
  else
    bb%ix_track_end = bb%ix_track_start - 1
    if (bb%ix_track_end == -1) bb%ix_track_end = branch%n_ele_track
  endif

else
  bb%ix_track_end = not_set$
  ele => tao_beam_track_endpoint (track_end, u%model%lat, int_str(ele%ix_branch), 'TRACK_END')
  if (.not. associated(ele)) return
  bb%ix_track_end = ele%ix_ele
endif

! Find where to save the beam at.
! Note: Beam will automatically be saved at fork elements and at the ends of the beam tracking.

if (u%beam%saved_at /= '') then
  call tao_locate_elements (u%beam%saved_at, u%ix_uni, eles, err, ignore_blank = .false.)
  if (err) then
    call out_io (s_warn$, r_name, 'BAD "saved_at" ELEMENT: ' // u%beam%saved_at)
  else
    do k = 1, size(eles)
      ele => eles(k)%ele
      if (ele%lord_status == super_lord$) ele => pointer_to_slave(ele, ele%n_slave)
      u%model_branch(ele%ix_branch)%ele(ele%ix_ele)%save_beam_internally = .true.
    enddo
  endif
endif

if (u%beam%dump_at /= '') then
  call tao_locate_elements (u%beam%dump_at, u%ix_uni, eles, err, ignore_blank = .false.)
  if (err) then
    call out_io (s_warn$, r_name, 'BAD "dump_at" ELEMENT: ' // u%beam%dump_at)
  else
    do k = 1, size(eles)
      ele => eles(k)%ele
      if (ele%lord_status == super_lord$) ele => pointer_to_slave(ele, ele%n_lord)
      u%model_branch(ele%ix_branch)%ele(ele%ix_ele)%save_beam_to_file = .true.
    enddo
  endif
endif

end subroutine tao_init_beam_in_universe

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

type (aperture_param_struct) :: da_param
type (tao_universe_struct), pointer :: u
type (tao_dynamic_aperture_struct), pointer :: da

real(rp) pz(100), a_emit, b_emit, ellipse_scale
integer :: ix_universe, ios, iu, i, j, n_pz, n_da

character(*) init_file
character(200) file_name
character(*), parameter :: r_name = 'tao_init_dynamic_aperture'

namelist / tao_dynamic_aperture / ix_universe, da_param, pz, a_emit, b_emit, ellipse_scale

!

if (init_file == '') return

call tao_open_file (init_file, iu, file_name, s_blank$)
if (iu == 0) then
  call out_io (s_blank$, r_name, "Note: Cannot open init file for tao_dynamic_aperture namelist read")
  return
endif

do n_da = 1, 1000
  ! Read tao_dynamic_aperture
  call out_io (s_blank$, r_name, 'Init: Begin reading tao_dynamic_aperture namelist')
  pz = real_garbage$
  a_emit = -1;  b_emit = -1
  ellipse_scale = 1
  ix_universe = -1
  da_param = aperture_param_struct()
  read (iu, nml = tao_dynamic_aperture, iostat = ios)
  if (ios > 0) then
    call out_io (s_error$, r_name, 'ERROR READING TAO_DYNAMIC_APERTURE NAMELIST.', &
                                   '[NOTE: THE SYNTAX WAS CHANGED 2021/6. PLEASE SEE THE TAO MANUAL.')
    rewind (iu)
    read (iu, nml = tao_dynamic_aperture)  ! To give error message
  endif

  if (ios < 0) then
    if (n_da == 1) call out_io (s_blank$, r_name, 'Note: No tao_dynamic_aperture namelist found')
    close(iu)
    return
  endif

  ! Count the list of pz
  ! Default if no pz set is 1 scan with pz = 0

  if (ix_universe > 0) then
    da => s%u(ix_universe)%dynamic_aperture
  else
    da => s%u(1)%dynamic_aperture
  endif

  do n_pz = 0, size(pz)-1
    if (pz(n_pz+1) == real_garbage$) exit
  enddo

  if (allocated(da%pz)) deallocate(da%pz)

  if (n_pz == 0) then
    allocate (da%pz(1))
    pz(1) = 0
  else
    allocate (da%pz(n_pz))
    da%pz = pz(1:n_pz)
  endif

  da%a_emit         = a_emit
  da%b_emit         = b_emit
  da%ellipse_scale  = ellipse_scale
  da%param          = da_param

  !

  if (ix_universe < 1) then
    do i = 2, size(s%u)
      s%u(i)%dynamic_aperture = da
    enddo
  endif
enddo

end subroutine tao_init_dynamic_aperture

end module
