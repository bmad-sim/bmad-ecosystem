!+
! Subroutine tao_init (err_flag)
!
! Subroutine to initialize the tao structures.
!
! Output:
!   err_flag  -- logical: Set Treu if there is an error. False otherwise.
!-

subroutine tao_init (err_flag)

use tao_init_mod, dummy2 => tao_init
use tao_init_data_mod, dummy3 => tao_init
use tao_init_variables_mod, dummy4 => tao_init
use tao_lattice_calc_mod, dummy5 => tao_init
use tao_command_mod, dummy8 => tao_init
use tao_set_mod, dummy9 => tao_init
use tao_plot_mod, only: tao_draw_plots
use tao_plot_window_mod, only: tao_destroy_plot_window
use random_mod

!$ use omp_lib

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

real(rp) value, sigma(6,6), merit
real(rp), pointer :: ptr_attrib

character(100) line, line2
character(200) plot_file, data_file, var_file, file_name, startup_file, hook_init_file
character(200) building_wall_file, beam_file, why_invalid, init_tao_file
character(40) name1, name2
character(16) :: r_name = 'tao_init'
character(16) init_name

integer i, i0, j, i2, j2, n_universes, iu, ix, ib, ip, ios, ie
integer iu_log

logical err_flag
logical err, calc_ok, valid_value, this_calc_ok, using_default

namelist / tao_start / startup_file, building_wall_file, hook_init_file, &
               data_file, var_file, plot_file, n_universes, init_name, beam_file

! Put all informational messages in the tao_init.log file.
! Only print error messages. Not standard ones.

err_flag = .true.
iu_log = -1

if (s%init%log_startup_arg /= '') then
  iu_log = lunget()
  open (iu_log, file = 'tao_init.log', action = 'write', iostat = ios)
  if (ios == 0) then
    call out_io (s_dinfo$, r_name, 'Opening initialization logging file: tao_init.log')
    call output_direct (iu_log, .true., s_blank$, s_abort$)
    call output_direct (iu_log, .false., s_blank$, s_success$) ! Do not print 
  else
    iu_log = -1
    call out_io (s_error$, r_name, 'NOTE: Cannot open a file for logging initialization messages')
  endif
else
  call output_direct (-1, .false., s_blank$, s_success$) ! Do not print 
endif

! OpenMP info

!$ s%global%n_threads = omp_get_max_threads()
!$ call out_io (s_important$, r_name, 'OpenMP active with number of threads: ' // int_str(s%global%n_threads))

! Open the init file.
! If the init file name is *not* the default (that is, it has been set by
! the user) then an open failure is considered fatal.
! Additionally, if there is an open failure and no lattice file has been specified
! by the user, then there is nothing to do and is considered fatal.

init_tao_file = 'tao.init'
if (s%init%init_file_arg /= '') init_tao_file = s%init%init_file_arg
if (s%init%noinit_arg /= '') init_tao_file = ''

iu = 0

if (init_tao_file /= '') then
  call tao_open_file (init_tao_file, iu, file_name, s_blank$)
  if (iu == 0) then ! If open failure
    if (s%init%init_file_arg == '') then
      call out_io (s_info$, r_name, 'Tao initialization file not found.')
      init_tao_file = ''
    else
      call output_direct (-1, print_and_capture=s%com%print_to_terminal)
      call out_io (s_error$, r_name, 'TAO INITIALIZATION FILE NOT FOUND: ' // init_tao_file)
      call tao_init_global('')
      return
    endif
  endif
endif

if (iu == 0 .and. (s%init%lattice_file_arg == '' .and. s%init%hook_lat_file == '')) then
  call output_direct (-1, print_and_capture=s%com%print_to_terminal)
  call out_io (s_abort$, r_name, &
          'Note: To run Tao, you either need a Tao initialization file or', &
          '  use a lattice file using the syntax "tao -lat <lat_file_name>".', &
          '  See the Tao manual for more details...')
  call tao_print_command_line_info
  return
endif

! Set defaults.
! n_universes is present to accomodate files with the old syntax.

n_universes        = 1                  ! set default
init_name          = "Tao"              ! set default
plot_file          = 'NOT SET!'         ! set default
data_file          = 'NOT SET!'         ! set default
var_file           = 'NOT SET!'         ! set default
beam_file          = 'NOT SET!'         ! set default
building_wall_file = 'NOT SET!'       
startup_file       = 'NOT SET!'       
hook_init_file     = 'NOT SET!'

! Quit the plot window if it exists so it will be recreated    

call tao_destroy_plot_window
s%com%init_plot_needed = .true.

! Read the global parameters

if (iu /= 0) then

  ! Read symbolic constants

  i = 0
  do
    i=i+1; read (iu, '(a)', iostat = ios) line
    if (ios /= 0) exit
    call string_trim (line, line, ix)
    if (line(1:ix) /= '&symbolic_number') cycle
    ! Found a symbolic number
    call string_trim (line(ix+1:), line, ix)
    i0 = i
    do
      ix = max(1, len_trim(line))  ! Prevent ix = 0 which will bomb next line.
      if (line(ix:ix) == '/') exit
      i=i+1; read (iu, '(a)', iostat = ios) line2
      if (ios /= 0) then
        call out_io (s_error$, r_name, 'CANNOT FIND ENDING "/" FOR "&sybolic_number" NAMELIST STARTING AT LINE \i0\ ', &
                      i_array = [i0])
        exit
      endif
      line = trim(line) // line2
    enddo
    j = index(line, '=')
    if (j == 0) then
      call out_io (s_error$, r_name, 'CANNOT FIND ENDING EQUAL SIGN "=" IN "&sybolic_number" NAMELIST STARTING AT LINE \i0\ ', &
                    i_array = [i0])
      exit
    endif
    call tao_set_symbolic_number_cmd(line(1:j-1), line(j+1:ix-1))
  enddo

  rewind(iu)

  ! Read tao_start namelist.

  read (iu, nml = tao_start, iostat = ios)

  if (ios < 0) then
    call out_io (s_info$, r_name, 'Cannot read "tao_start" namelist in file: ' // file_name)
  endif

  if (ios > 0) then
    call out_io (s_error$, r_name, 'ERROR IN READING "TAO_START" NAMELIST IN FILE: ' // file_name)
    rewind (iu)
    read (iu, nml = tao_start)  ! And generate error message.    
  endif

  close (iu)
  s%init%init_name = init_name
  s%com%n_universes = n_universes
endif

! Set

call set_this_file_name (plot_file, init_tao_file, s%init%plot_file_arg, s%init%hook_plot_file)
call set_this_file_name (data_file, init_tao_file, s%init%data_file_arg, s%init%hook_data_file)
call set_this_file_name (var_file,  init_tao_file, s%init%var_file_arg, s%init%hook_var_file)
call set_this_file_name (beam_file, init_tao_file, s%init%beam_file_arg, s%init%hook_beam_file)
call set_this_file_name (building_wall_file, '',   s%init%building_wall_file_arg, s%init%hook_building_wall_file)
call set_this_file_name (startup_file, 'tao.startup', s%init%startup_file_arg, s%init%hook_startup_file)
call set_this_file_name (hook_init_file, 'tao_hook.init', s%init%hook_init_file_arg, '')
s%init%hook_init_file = hook_init_file

! Tao inits.
! Data can have variable info so init vars first.

if (allocated(s%u)) call deallocate_everything ()

global_com%exit_on_error = .false.

call tao_init_global(init_tao_file) ! The first global read is just to look for %debug and %stop_on_error values.
call tao_init_lattice (init_tao_file, err); if (err) return
call tao_init_global(init_tao_file)
call tao_init_dynamic_aperture (init_tao_file)
call tao_init_beams (beam_file)
call tao_init_variables (var_file)
call tao_init_data (data_file)
call tao_init_building_wall (building_wall_file)

if (associated(tao_hook_init1_ptr)) call tao_hook_init1_ptr (init_tao_file)

! Seed random number generator

if (s%global%random_seed /= -1) call ran_seed_put (s%global%random_seed)
if (s%global%random_engine /= '') call ran_engine (s%global%random_engine)
call ran_gauss_converter (s%global%random_gauss_converter)
if (s%global%random_sigma_cutoff > 0) call ran_gauss_converter (set_sigma_cut = s%global%random_sigma_cutoff)

! check variables
! check if vars are good

call tao_set_var_useit_opt

do i = 1, s%n_var_used
  var => s%var(i)
  if (.not. var%exists) cycle
  do j = 1, size(var%slave)
    var_slave => var%slave(j)
    u => s%u(var_slave%ix_uni)
    if (var_slave%ix_ele < 0) cycle  ! Do not check EG "particle_start".
    if (.not. attribute_free (var_slave%ix_ele, var_slave%ix_branch, var%attrib_name, &
                                                                u%model%lat, dependent_attribs_free = .true.)) then
      call out_io (s_error$, r_name, &
                '       VARIABLE TRYING TO CONTROL AN ATTRIBUTE THAT IS NOT FREE TO VARY.', &
                '       VARIABLE:  ' // tao_var1_name(var), &
                '       ELEMENT:   ' // var%ele_name, &
                '       ATTRIBUTE: ' // var%attrib_name)
      var%exists = .false.
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
        if (associated (s%var(i)%slave(j)%model_value, &
                          s%var(i2)%slave(j2)%model_value)) then
          write (name1, '(2a, i0, a)') trim(s%var(i)%v1%name), '[', s%var(i)%ix_v1, ']'  
          write (name2, '(2a, i0, a)') trim(s%var(i2)%v1%name), '[', s%var(i2)%ix_v1, ']'  
          call out_io (s_important$, r_name, &
               '       VARIABLE:     ' // name1, &
               '       AND VARIABLE: ' // name2, &
               '       CONTROL THE SAME EXACT THING!', &
               '       THIS CAN CAUSE STRANGE BEHAVIOR. YOU HAVE BEEN WARNED!!!')            
        endif
      enddo
    enddo
  enddo
enddo

! plotting

call tao_init_plotting (plot_file)


! Route all messages back to the terminal.
! Need to do this before calling tao_lattice_calc since we don't want to supress these messages.

call output_direct (-1, print_and_capture=s%com%print_to_terminal)

! Set up model and base lattices.
! Must first transfer to model lattice for tao_lattice_calc to run.

if (bmad_com%radiation_fluctuations_on .and. s%global%track_type == 'single') then
  call out_io (s_info$, r_name, 'Note: Radiation fluctuations (but not necessarily damping) are always turned off for single particle tracking...')
endif

! Calculate radiation integrals.

if (iu_log > 0) write (iu_log, '(a)') '*Init: Starting radiation and chrom calc.'

do i = lbound(s%u, 1), ubound(s%u, 1)
  u => s%u(i)
  if (u%design_same_as_previous) cycle

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
    call tao_single_track (tao_lat, this_calc_ok, ib)

    if (branch%param%particle /= photon$ .and. this_calc_ok) then
      if (s%global%rad_int_user_calc_on) call radiation_integrals (tao_lat%lat, tao_branch%orbit, &
                                  tao_branch%modes_ri, tao_branch%ix_rad_int_cache, ib, tao_lat%rad_int_by_ele_ri)

      if (tao_branch%track_state == moving_forward$) then
        if (s%global%rad_int_user_calc_on) call emit_6d(branch%ele(0), .false., tao_branch%modes_6d, sigma, tao_branch%orbit, tao_lat%rad_int_by_ele_6d)
        if (s%global%rad_int_user_calc_on) call emit_6d(branch%ele(0), .true., tao_branch%modes_6d, sigma, tao_branch%orbit, tao_lat%rad_int_by_ele_6d)
        call chrom_calc (tao_lat%lat, s%global%delta_e_chrom, tao_branch%a%chrom, tao_branch%b%chrom, err, tao_branch%orbit(0)%vec(6), &
                              tao_lat%low_E_lat, tao_lat%high_E_lat, tao_branch%low_E_orb, tao_branch%high_E_orb, ib, tao_branch%orbit(0))
      endif

      if (branch%param%geometry == closed$ .and. tao_branch%track_state == moving_forward$) then
        tao_branch%modes_6d%momentum_compaction = momentum_compaction(branch)
        if (tao_branch%modes_6d%a%j_damp < 0 .or. tao_branch%modes_6d%b%j_damp < 0 .or. &
                                                   (tao_branch%modes_6d%z%j_damp < 0 .and. rf_is_on(branch))) then
        call out_io (s_info$, r_name, &
          'Negative damping partition number detected therefore the lattice is unstable with radiation excitations.', &
          'Note1: This may not be a problem if the amount of radiation generated is low (like for protons).', &
          'Note2: Instability with respect to radiation excitations does not affect such things as the closed orbit calculation.')
        endif
      endif
    endif

  enddo
enddo

! Turn off RF if needed. But first calculate the synchrotron tune.

do i = lbound(s%u, 1), ubound(s%u, 1)
  u => s%u(i)

  do ib = 0, ubound(u%model%lat%branch, 1)
    branch => u%model%lat%branch(ib)
    if (branch%param%geometry == closed$) then
      call calc_z_tune(u%model%lat%branch(ib))
      if (.not. s%global%rf_on) call set_on_off (rfcavity$, u%model%lat, off$, ix_branch = ib)
    endif
    do ie = 1, branch%n_ele_max
      ! Make sure cache is calculated with respect to particle orbit.
      if (associated(branch%ele(ie)%rad_map)) branch%ele(ie)%rad_map%stale = .true.
    enddo
  enddo
enddo

!

if (iu_log > 0) write (iu_log, '(a)') '*Init: Starting lattice calc.'

s%u%calc%lattice = .true.
call tao_lattice_calc (calc_ok)

do i = lbound(s%u, 1), ubound(s%u, 1)
  u => s%u(i)

  u%design = u%model
  u%base = u%design

  u%design%name = 'design'
  u%model%name  = 'model'
  u%base%name   = 'base'

  u%data%design_value = u%data%model_value
  u%data%base_value   = u%data%model_value
  u%data%good_design  = u%data%good_model
  u%data%good_base    = u%data%good_model
enddo

if (iu_log > 0) write (iu_log, '(a)') '*Init: End lattice calc.'

call tao_var_repoint()

! Normally you will want to use tao_hook_init1. However, tao_hook_init2 can be used, for example, 
! to set model variable values different from design variable values.

if (associated(tao_hook_init2_ptr)) call tao_hook_init2_ptr ()

! Draw everything

if (iu_log > 0) write (iu_log, '(a)') '*Init: Draw plots.'

call tao_plot_setup ()     ! transfer data to the plotting structures
call tao_draw_plots ()     ! Update the plotting window

! Print bad data

if (iu_log > 0) write (iu_log, '(a)') '*Init: Print bad data.'

do i = lbound(s%u, 1), ubound(s%u, 1)
  if (.not. s%u(i)%is_on .or. .not. s%global%lattice_calc_on) cycle
  do j = 1, size(s%u(i)%data)
    data => s%u(i)%data(j)
    if (data%exists .and. data%data_type /= 'null' .and. .not. data%good_model) then
      call tao_evaluate_a_datum (data, s%u(i), s%u(i)%model, value, valid_value, why_invalid)
      call tao_set_invalid (data, why_invalid)
    endif
  enddo
enddo

! Look for a startup file

if (iu_log > 0) write (iu_log, '(a)') '*Init: Execute startup file.'

if (startup_file /= '' .and. s%init%nostartup_arg == '') then
  using_default = (startup_file == 'tao.startup')
  call tao_open_file (startup_file, iu, file_name, -1)
  if (iu == 0 .and. using_default) then ! If default
    startup_file = trim(s%init%init_file_arg_path) // 'tao.startup'
    call tao_open_file (startup_file, iu, file_name, -1)
  endif

  if (iu == 0 .and. .not. using_default) then
    call out_io (s_error$, r_name, 'Tao startup file not found: ' // file_name)

  elseif (iu /= 0) then
    close (iu)
    call out_io (s_blank$, r_name, 'Using startup file: ' // file_name)
    s%com%cmd_from_cmd_file = .false.
    call tao_cmd_history_record ('call ' // startup_file)
    call tao_call_cmd (file_name)
  endif
endif

! Bookkeeping

if (iu_log > 0) write (iu_log, '(a)') '*Init: End bookkeeping.'

call tao_set_data_useit_opt()
call tao_set_var_useit_opt()
err_flag = .false.

if (iu_log > 0) write (iu_log, '(a)') '*Init: And done.'
if (iu_log > 0) close (iu_log)

merit = tao_merit()  ! To calc initial merit contributions to the merit function.

!------------------------------------------------------------------------------
contains

! every pointer and allocatable needs to be deallocated now before the universe
! is reallocated.

subroutine deallocate_everything ()

use radiation_mod, only: release_rad_int_cache

implicit none

type (tao_plot_struct), pointer :: plot
type (tao_curve_struct), pointer :: curve
type (tao_universe_struct), pointer :: u

integer i, j, k, ib, istat

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

if (allocated(s%plot_page%template)) then
  do i = 1, size(s%plot_page%template)
    plot => s%plot_page%template(i)
    if (.not. allocated (plot%graph)) cycle
    deallocate(plot%graph, stat=istat)
  enddo
  deallocate (s%plot_page%template)
endif

if (allocated(s%plot_page%lat_layout%ele_shape)) deallocate (s%plot_page%lat_layout%ele_shape)
if (allocated(s%plot_page%floor_plan%ele_shape)) deallocate (s%plot_page%floor_plan%ele_shape)
if (allocated(s%plot_page%pattern))              deallocate (s%plot_page%pattern)
if (allocated(s%com%covar))                      deallocate (s%com%covar, s%com%alpha)

! Universes 

if (allocated (s%u)) then
  do i = lbound(s%u, 1), ubound(s%u, 1)

    u => s%u(i)
    ! radiation integrals cache
    if (allocated(u%model%tao_branch)) then
      do ib = 0, ubound(u%model%tao_branch, 1)
        if (u%model%tao_branch(ib)%ix_rad_int_cache /= 0) call release_rad_int_cache(u%model%tao_branch(ib)%ix_rad_int_cache)
        if (u%design%tao_branch(ib)%ix_rad_int_cache /= 0) call release_rad_int_cache(u%design%tao_branch(ib)%ix_rad_int_cache)
        if (u%base%tao_branch(ib)%ix_rad_int_cache /= 0) call release_rad_int_cache(u%base%tao_branch(ib)%ix_rad_int_cache)
      enddo

      ! Orbits

      deallocate(u%model%tao_branch, stat=istat)
      deallocate(u%design%tao_branch, stat=istat)
      deallocate(u%base%tao_branch, stat=istat)
    endif

    ! Beams

    if (associated(u%model_branch)) then
      do ib = 0, ubound(u%model_branch, 1)
        call reallocate_beam(u%model_branch(ib)%ele(0)%beam, 0, 0)
        deallocate (u%model_branch(ib)%ele)
        call reallocate_beam(u%model_branch(ib)%beam%beam_at_start, 0, 0)
      enddo
      deallocate (u%model_branch)
    endif

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

!+
! Subroutine set_this_file_name (file_name, default_name, arg_name, hook_name)
!
! Routine to set the name of the file based on the file name set from various sources.
!
! Input:
!   file_name     -- character(*): The file name as set in the tao init file.
!                       If it has not been set then file_name = 'NOT SET!'.
!   default_name  -- character(*): Default name if no other name is present.
!   arg_name      -- character(*): This name is from the startup command line.
!   hook_name     -- character(*): Name as set by the tao_hook_parse_command_args routine.
!                       Essentially the hook name overrides the default name.
!
! Output:
!   file_name     -- character(*): File name.
!-

subroutine set_this_file_name (file_name, default_name, arg_name, hook_name)

character(*) file_name, default_name, arg_name, hook_name
character(200) name

! Order of preference. Highest used first:
!   1) Name has been set on the command line.
!   2) Name has been set in the Tao init file.
!   3) Name has been set via a hook routine.
!   4) Default_name.

if (arg_name /= '') then
  name  = arg_name
elseif (file_name /= 'NOT SET!') then
  name = file_name
  if (file_name_is_relative(name)) name = trim(s%init%init_file_arg_path) // trim(name)
elseif (hook_name /= '') then
  name = hook_name
elseif (default_name /= '') then
  name = default_name
else
  name = default_name
endif

file_name = name


end subroutine set_this_file_name

end subroutine tao_init


