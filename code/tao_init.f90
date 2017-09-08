!+
! Subroutine tao_init (err_flag)
!
! Subroutine to initialize the tao structures.
!
! Input:
!
! Output:
!-

subroutine tao_init (err_flag)

use tao_mod, dummy => tao_init
use tao_init_mod, dummy2 => tao_init
use tao_init_data_mod, dummy3 => tao_init
use tao_init_variables_mod, dummy4 => tao_init
use tao_lattice_calc_mod, dummy5 => tao_init
use tao_plot_mod, dummy6 => tao_init
use tao_data_and_eval_mod, dummy7 => tao_init
use tao_command_mod, dummy8 => tao_init

implicit none

type (tao_universe_struct), pointer :: u
type (tao_var_struct), pointer :: var
type (tao_var_slave_struct), pointer :: var_slave
type (tao_plot_struct), pointer :: p
type (tao_data_struct), pointer :: data
type (beam_struct), pointer :: beam
type (tao_lattice_struct), pointer :: tao_lat
type (branch_struct), pointer :: branch
type (tao_lattice_branch_struct), pointer :: tao_branch

real(rp) value
real(rp), pointer :: ptr_attrib

character(80) arg, arg2, startup_file, hook_init_file
character(100) plot_file, data_file, var_file, file_name
character(100) building_wall_file, beam_file, why_invalid, init_tao_file
character(40) name1, name2
character(16) :: r_name = 'tao_init'
character(16) init_name

integer i, j, i2, j2, n_universes, iu, ix, n_arg, ib, ip, ios
integer iu_log

logical err, calc_ok, valid_value, this_calc_ok
logical :: err_flag

namelist / tao_start / startup_file, building_wall_file, hook_init_file, &
               data_file, var_file, plot_file, n_universes, init_name, beam_file

! global inits

s%com%n_alias = 0
s%com%ix_key_bank = 0             ! For single mode.
s%com%use_saved_beam_in_tracking = .false.
if (.not. allocated(s%com%cmd_file)) allocate (s%com%cmd_file(0:0))

s%global%optimizer_allow_user_abort = (isatty(0) == 1)  ! Allow abort if input from tty (instead of a file).

call getenv ('ACC_PLOT_DISPLAY_TYPE', name1)
if (name1 /= '') s%plot_page%plot_display_type = name1

! Put all informational messages in the tao_init.log file.
! Only print error messages. Not standard ones.

iu_log = -1
if (s%com%log_startup) then
  iu_log = lunget()
  open (iu_log, file = 'tao_init.log', action = 'write', iostat = ios)
  if (ios == 0) then
    call out_io (s_dinfo$, r_name, 'Opening initialization logging file: tao_init.log')
    call output_direct (iu_log, .true., 0, s_blank$, s_abort$)
    call output_direct (iu_log, .false., 0, s_blank$, s_success$) ! Do not print 
  else
    iu_log = -1
    call out_io (s_error$, r_name, 'NOTE: Cannot open a file for logging initialization messages')
  endif
else
  call output_direct (-1, .false., 0, s_blank$, s_success$) ! Do not print 
endif

! Open the init file.
! If the init file name is *not* the default (that is, it has been set by
! the user) then an open failure is considered fatal.
! Additionally, if there is an open failure and no lattice file has been specified
! by the user, then there is nothing to do and is considered fatal.

err_flag = .true.

iu = 0
if (s%com%init_tao_file /= '') then
  call tao_open_file (s%com%init_tao_file, iu, file_name, s_blank$)
  if (iu == 0) then ! If open failure
    call out_io (s_info$, r_name, 'Tao initialization file not found.')
    if (s%com%lat_file == '' .or. s%com%init_tao_file_arg_set) then
      call output_direct (-1, do_print=s%com%print_to_terminal)
      call out_io (s_blank$, r_name, &
              'Note: To run Tao, you either need a Tao initialization file or', &
              '  use a lattice file using the syntax "tao -lat <lat_file_name>".', &
              '  See the Tao manual for more details...')
      call tao_print_command_line_info
      stop
    endif
    s%com%init_tao_file = ''
  endif
endif

! Set defaults.
! n_universes is present to accomodate files with the old syntax.

init_tao_file = s%com%init_tao_file

plot_file          = 'NOT SET!'         ! set default
data_file          = 'NOT SET!'         ! set default
var_file           = 'NOT SET!'         ! set default
beam_file          = 'NOT SET!'         ! set default
building_wall_file = 'NOT SET!'       
n_universes        = 1                  ! set default
init_name          = "Tao"              ! set default
startup_file       = 'NOT SET!'       
hook_init_file     = 'NOT SET!'

! Read the info

if (iu /= 0) then
  read (iu, nml = tao_start, iostat = ios)

  if (ios < 0) then
    call out_io (s_info$, r_name, 'Cannot read "tao_start" namelist in file: ' // file_name)
  endif

  if (ios > 0) then
    call out_io (s_abort$, r_name, 'ERROR IN READING "TAO_START" NAMELIST IN FILE: ' // file_name)
    rewind (iu)
    read (iu, nml = tao_start)  ! And generate error message.    
  endif

  close (iu)
  s%com%init_name = init_name
  s%com%n_universes = n_universes
endif

! Set

call set_this_file_name (plot_file, init_tao_file, s%com%plot_file)
call set_this_file_name (data_file, init_tao_file, s%com%data_file)
call set_this_file_name (var_file,  init_tao_file, s%com%var_file)
call set_this_file_name (beam_file, init_tao_file, s%com%beam_file)
call set_this_file_name (building_wall_file, '',   s%com%building_wall_file)
call set_this_file_name (startup_file, 'tao.startup', s%com%startup_file)
call set_this_file_name (hook_init_file, 'tao_hook.init', s%com%hook_init_file)
s%com%hook_init_file = hook_init_file

! Tao inits.
! Data can have variable info so init vars first.

if (allocated(s%u)) call deallocate_everything ()

global_com%exit_on_error = .false.

call tao_init_global(init_tao_file)
call tao_init_lattice (init_tao_file)
call tao_init_dynamic_aperture (init_tao_file)
call tao_init_beams (beam_file)
call tao_init_variables (var_file)
call tao_init_data (data_file)
call tao_init_building_wall (building_wall_file)

call tao_hook_init1 (init_tao_file)

! check variables
! check if vars are good

call tao_set_var_useit_opt

do i = 1, s%n_var_used
  var => s%var(i)
  if (.not. var%exists) cycle
  do j = 1, size(var%slave)
    var_slave => var%slave(j)
    u => s%u(var_slave%ix_uni)
    if (var_slave%ix_ele < 0) cycle  ! Do not check EG "beam_start".
    if (.not. attribute_free (var_slave%ix_ele, var_slave%ix_branch, var%attrib_name, u%model%lat)) then
      call out_io (s_abort$, r_name, &
                'ERROR: VARIABLE TRYING TO CONTROL AN ATTRIBUTE THAT IS NOT FREE TO VARY.', &
                '       VARIABLE:  ' // tao_var1_name(var), &
                '       ELEMENT:   ' // var%ele_name, &
                '       ATTRIBUTE: ' // var%attrib_name)
      call err_exit
    endif
  enddo
enddo

! make sure two variables do not vary the same attribute

do i = 1, s%n_var_used
  if (.not. allocated(s%var(i)%slave)) cycle
  do j = 1, size(s%var(i)%slave)
    do i2 = i, s%n_var_used
      if (.not. allocated(s%var(i2)%slave)) cycle
      do j2 = 1, size(s%var(i2)%slave)
        if (i == i2 .and. j == j2) cycle
        if (s%com%common_lattice .and. &
                          s%var(i)%slave(j)%ix_uni /= s%var(i2)%slave(j2)%ix_uni) cycle
        if (associated (s%var(i)%slave(j)%model_value, &
                          s%var(i2)%slave(j2)%model_value)) then
          write (name1, '(2a, i0, a)') trim(s%var(i)%v1%name), '[', s%var(i)%ix_v1, ']'  
          write (name2, '(2a, i0, a)') trim(s%var(i2)%v1%name), '[', s%var(i2)%ix_v1, ']'  
          call out_io (s_error$, r_name, &
               'ERROR: VARIABLE:     ' // name1, &
               '       AND VARIABLE: ' // name2, &
               '       CONTROL THE SAME EXACT THING!', &
               '       YOU HAVE BEEN WARNED!!!')            
        endif
      enddo
    enddo
  enddo
enddo

! plotting

call tao_init_plotting (plot_file)
  
! Close the log file and route all messages back to the terminal.
! Need to do this before calling tao_lattice_calc since we don't want to supress these messages.

if (iu_log > -1) close (iu_log)
call output_direct (-1, do_print=s%com%print_to_terminal)

! Set up model and base lattices.
! Must first transfer to model lattice for tao_lattice_calc to run.

if (bmad_com%radiation_fluctuations_on .and. s%global%track_type == 'single') then
  call out_io (s_info$, r_name, 'Note: Radiation fluctuations are always turned off for single particle tracking...')
endif

! Calculate radiation integrals.

do i = lbound(s%u, 1), ubound(s%u, 1)
  u => s%u(i)
  tao_lat => u%model  ! In the past tao_lat could point to design or base but no more.

  do ib = 0, ubound(tao_lat%lat%branch, 1)
    branch => tao_lat%lat%branch(ib)
    tao_branch => tao_lat%tao_branch(ib)

    call tao_data_coupling_init(branch)

    if (.not. branch%param%live_branch) cycle
    
    do j = 1, 6
      tao_lat%tao_branch(ib)%orbit%vec(j) = 0.0
    enddo

    call tao_inject_particle (u, tao_lat, ib)
    call tao_single_track (u, tao_lat, this_calc_ok, ib)

    call radiation_integrals (tao_lat%lat, tao_branch%orbit, &
                          tao_branch%modes, tao_branch%ix_rad_int_cache, ib, tao_branch%rad_int)

    tao_branch%modes_rf_on = tao_branch%modes
    tao_branch%rad_int_rf_on = tao_branch%rad_int_rf_on

    call chrom_calc (tao_lat%lat, s%global%delta_e_chrom, tao_branch%a%chrom, &
                         tao_branch%b%chrom, err, low_E_lat=tao_branch%low_E_lat, high_E_lat=tao_branch%high_E_lat)
  enddo
enddo

! Turn off RF

do i = lbound(s%u, 1), ubound(s%u, 1)
  u => s%u(i)
  do ib = 0, ubound(u%design%lat%branch, 1)
    if (u%design%lat%branch(ib)%param%geometry == closed$) then
      call calc_z_tune(u%design%lat, ib)
      if (.not. s%global%rf_on) then
        call out_io (s_info$, r_name, "Note: global%rf_on = False  -->  RFCavities will be turned off in lattices")
        call set_on_off (rfcavity$, u%model%lat, off$, ix_branch = ib)
        u%model%tao_branch(0)%orb0 = u%model%lat%beam_start
      endif
    endif
  enddo
enddo

!

s%u%calc%lattice = .true.
call tao_lattice_calc (calc_ok)

do i = lbound(s%u, 1), ubound(s%u, 1)
  u => s%u(i)
  u%design = u%model
  u%base = u%design
  u%design%tao_branch = u%model%tao_branch
  u%base%tao_branch   = u%design%tao_branch
  u%data%design_value = u%data%model_value
  u%data%base_value   = u%data%model_value
  u%data%good_design  = u%data%good_model
  u%data%good_base    = u%data%good_model
  u%design%tao_branch%modes = u%design%tao_branch%modes_rf_on
enddo

call tao_var_repoint()

! Normally you will want to use tao_hook_init1. However, tao_hook_init2 can be used, for example, 
! to set model variable values different from design variable values.

call tao_hook_init2 ()     

! Draw everything

call tao_plot_setup ()     ! transfer data to the plotting structures
call tao_draw_plots ()     ! Update the plotting window

! Print bad data

do i = lbound(s%u, 1), ubound(s%u, 1)
  do j = 1, size(s%u(i)%data)
    data => s%u(i)%data(j)
    if (data%exists .and. .not. data%good_model) then
      call tao_evaluate_a_datum (data, s%u(i), s%u(i)%model, value, valid_value, why_invalid)
      call out_io(s_warn$, r_name, &
                  'DATUM EXISTS BUT CANNOT COMPUTE A MODEL VALUE: ' // tao_datum_name(data), &
                  '   INVALID SINCE: ' // why_invalid)
    endif
  enddo
enddo

! Look for a startup file

if (startup_file /= '') then
  call tao_open_file (startup_file, iu, file_name, -1)
  if (iu /= 0) then
    close (iu)
    call out_io (s_blank$, r_name, 'Using startup file: ' // file_name)
    s%com%cmd_from_cmd_file = .false.
    call tao_cmd_history_record ('call ' // startup_file)
    call tao_call_cmd (file_name)
  else if (startup_file /= 'tao.startup') then  ! If not default
    call out_io (s_error$, r_name, 'Tao startup file not found: ' // file_name)
  endif
endif

! Bookkeeping

call tao_set_data_useit_opt()
call tao_set_var_useit_opt()
err_flag = .false.

!------------------------------------------------------------------------------
contains

! every pointer and allocatable needs to be deallocated now before the universe
! is reallocated.

subroutine deallocate_everything ()

implicit none

type (tao_plot_struct), pointer :: plot
type (tao_curve_struct), pointer :: curve
type (tao_universe_struct), pointer :: u

integer i, j, k, istat

! building walls

if (allocated(s%building_wall%section)) deallocate (s%building_wall%section)

! Variables  

if (allocated (s%v1_var)) then
  deallocate(s%v1_var, stat=istat)
endif
  
if (allocated (s%var)) then
  do i = lbound(s%var,1), ubound(s%var,1)
    deallocate(s%var(i)%slave, stat=istat)
  enddo
  deallocate(s%var, stat=istat)
endif

! Keytable 

if (allocated(s%key)) deallocate(s%key, stat=istat)

! Plotting  

if (allocated(s%plot_page%region)) deallocate (s%plot_page%region)

do i = 1, size(s%plot_page%template)
  plot => s%plot_page%template(i)
  if (.not. allocated (plot%graph)) cycle
  deallocate(plot%graph, stat=istat)
enddo
deallocate (s%plot_page%template)

if (allocated(s%plot_page%lat_layout%ele_shape)) deallocate (s%plot_page%lat_layout%ele_shape)
if (allocated(s%plot_page%floor_plan%ele_shape)) deallocate (s%plot_page%floor_plan%ele_shape)
if (allocated(s%com%covar))                deallocate (s%com%covar, s%com%alpha)

! Universes 

if (allocated (s%u)) then
  do i = lbound(s%u, 1), ubound(s%u, 1)

    u => s%u(i)
    ! radiation integrals cache
    do ib = 0, ubound(u%model%lat%branch, 1)
      if (u%model%tao_branch(ib)%ix_rad_int_cache /= 0) call release_rad_int_cache(u%model%tao_branch(ib)%ix_rad_int_cache)
      if (u%design%tao_branch(ib)%ix_rad_int_cache /= 0) call release_rad_int_cache(u%design%tao_branch(ib)%ix_rad_int_cache)
      if (u%base%tao_branch(ib)%ix_rad_int_cache /= 0) call release_rad_int_cache(u%base%tao_branch(ib)%ix_rad_int_cache)
    enddo

    ! Orbits

    deallocate(u%model%tao_branch, stat=istat)
    deallocate(u%design%tao_branch, stat=istat)
    deallocate(u%base%tao_branch, stat=istat)

    ! Beams: All s%u(i)%ele point to the same place with common_lattice.

    if (i == 0 .or. .not. s%com%common_lattice) then
      do ib = 0, ubound(u%uni_branch, 1)
        call reallocate_beam(u%uni_branch(ib)%ele(0)%beam, 0, 0)
        deallocate (u%uni_branch(ib)%ele)
      enddo
      deallocate (u%uni_branch)
    endif

    call reallocate_beam(u%beam%start, 0, 0)

    ! Lattices

    call deallocate_lat_pointers (u%model%lat)
    call deallocate_lat_pointers (u%design%lat)
    call deallocate_lat_pointers (u%base%lat)

  enddo

  deallocate (s%u)

endif
    
end subroutine deallocate_everything
    
!------------------------------------------------------------------------------
! contains

subroutine set_this_file_name (file_name, init_name, tao_com_name)

character(*) file_name, init_name, tao_com_name

! file_name may already have been set from the tao_init file. If not, it is 'NOT SET!'.
! tao_com_name comes from the command line.
! init_name is the default if not set.

if (tao_com_name /= '') then
  file_name    = tao_com_name
elseif (file_name == 'NOT SET!') then
  file_name = init_name
elseif (file_name_is_relative(file_name)) then
  file_name = trim(s%com%init_tao_file_path) // trim(file_name)
endif

end subroutine set_this_file_name

end subroutine tao_init


