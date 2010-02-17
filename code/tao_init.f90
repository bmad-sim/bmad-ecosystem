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

use tao_mod
use tao_lattice_calc_mod
use tao_command_mod
use tao_plot_mod
use tao_init_mod

implicit none

type (tao_universe_struct), pointer :: u
type (tao_var_struct), pointer :: var
type (tao_this_var_struct), pointer :: this
type (tao_plot_struct), pointer :: p
type (tao_data_struct), pointer :: data
type (beam_struct), pointer :: beam

real(rp), pointer :: ptr_attrib

character(80) arg, arg2, startup_file
character(100) lattice_file, plot_file, data_file, var_file, file_name
character(100) wall_file
character(40) name1, name2
character(16) :: r_name = 'tao_init'
character(16) init_name

integer i, j, i2, j2, n_universes, iu, ix, n_arg, ib, ip, ios
integer iu_log

logical err, calc_ok  
logical, optional :: err_flag

namelist / tao_start / lattice_file, startup_file, wall_file, &
               data_file, var_file, plot_file, n_universes, init_name

! global inits

tao_com%n_alias = 0
tao_com%ix_key_bank = 0             ! For single mode.
tao_com%use_saved_beam_in_tracking = .false.
if (.not. allocated(tao_com%cmd_file)) allocate (tao_com%cmd_file(0:0))

! Put all informational messages in the tao_init.log file.
! Only print error messages. Not standard ones.

iu_log = lunget()
open (iu_log, file = 'tao_init.log', action = 'write', iostat = ios)
if (ios == 0) then
  call out_io (s_dinfo$, r_name, 'Opening initialization logging file: tao_init.log')
  call output_direct (iu_log, .true., s_blank$, s_abort$)
  call output_direct (iu_log, .false., s_blank$, s_success$) ! Do not print 
else
 call out_io (s_error$, r_name, &
                'NOTE: Cannot open a file for logging initialization information')
endif

! Open the init file.
! If the init file name is *not* the default (that is, it has been set by
! the user) then an open failure is considered fatal.
! Additionally, if there is an open failure and no lattice file has been specified
! by the user, then there is nothing to do and is considered fatal.

if (present(err_flag)) err_flag = .true.

call tao_open_file ('TAO_INIT_DIR', tao_com%init_tao_file, iu, file_name)
if (iu == 0) then ! If open failure
  call out_io (s_warn$, r_name, 'CANNOT OPEN TAO INITIALIZATION FILE!')
  if (tao_com%init_tao_file /= tao_com%default_init_tao_file .or. &
                                    tao_com%init_lat_file == '') stop
  tao_com%init_tao_file = ''
endif

! Set defaults.
! n_universes is present to accomodate files with the old syntax.

lattice_file       = tao_com%init_tao_file      ! set default
plot_file          = tao_com%init_tao_file      ! set default
data_file          = tao_com%init_tao_file      ! set default
var_file           = tao_com%init_tao_file      ! set default
wall_file          = ''
n_universes        = 1              ! set default
init_name          = "Tao"          ! set default
startup_file       = "tao.startup"

! Read the info

if (iu /= 0) then
  read (iu, nml = tao_start, iostat = ios)
  close (iu)

  if (ios /= 0) then
    call out_io (s_info$, r_name, &
                    'Cannot read "tao_start" namelist in file: ' // file_name)
  endif

  tao_com%init_name = init_name
  tao_com%n_universes = n_universes
endif

! Tao inits.
! Data can have variable info so init vars first.

if (allocated(s%u)) call deallocate_everything ()

bmad_status%exit_on_error = .false.

call tao_init_global(tao_com%init_tao_file)
call tao_init_lattice (lattice_file) 
call tao_init_beams_and_uni_connections (tao_com%init_tao_file)
call tao_init_variables (var_file)
call tao_init_data (data_file)
call tao_init_wall (wall_file)

call tao_hook_init1 (tao_com%init_tao_file)

! check variables
! check if vars are good

call tao_set_var_useit_opt

do i = 1, size(s%var)
  var => s%var(i)
  if (.not. var%exists) cycle
  do j = 1, size(var%this)
    this => var%this(j)
    u => s%u(this%ix_uni)
    if (.not. attribute_free (this%ix_ele, this%ix_branch, var%attrib_name, u%model%lat)) then
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

do i = 1, size(s%var)
  if (.not. allocated(s%var(i)%this)) cycle
  do j = 1, size(s%var(i)%this)
    do i2 = i, size(s%var)
      if (.not. allocated(s%var(i2)%this)) cycle
      do j2 = 1, size(s%var(i2)%this)
        if (i == i2 .and. j == j2) cycle
        if (tao_com%common_lattice .and. &
                          s%var(i)%this(j)%ix_uni /= s%var(i2)%this(j2)%ix_uni) cycle
        if (associated (s%var(i)%this(j)%model_value, &
                          s%var(i2)%this(j2)%model_value)) then
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
! Need to do this before calling tao_lattice_calc since we don't want to supress messages.

close (iu_log)
call output_direct (0, .true.)

! Init radiation damping and fluctuation constants if needed.

if (bmad_com%radiation_damping_on .or. bmad_com%radiation_fluctuations_on) then
  do i = 1, size(s%u)
    call setup_radiation_tracking (s%u(i)%model%lat)
  enddo
endif

! Set up model and base lattices.
! Must first transfer to model lattice for tao_lattice_calc to run.

s%u%lattice_recalc = .true.
call tao_lattice_calc (calc_ok) 

do i = lbound(s%u, 1), ubound(s%u, 1)
  s%u(i)%design = s%u(i)%model
  s%u(i)%base = s%u(i)%design
  s%u(i)%design%lat_branch = s%u(i)%model%lat_branch
  s%u(i)%base%lat_branch   = s%u(i)%design%lat_branch
  s%u(i)%data%design_value = s%u(i)%data%model_value
  s%u(i)%data%base_value   = s%u(i)%data%model_value
  s%u(i)%data%good_design  = s%u(i)%data%good_model
  s%u(i)%data%good_base    = s%u(i)%data%good_base
enddo

! tao_hook_init2 is for custom setup after the regular setup

call tao_hook_init2 ()     

! Draw everything

call tao_plot_setup ()     ! transfer data to the plotting structures
call tao_draw_plots ()     ! Update the plotting window

! Print bad data

do i = lbound(s%u, 1), ubound(s%u, 1)
  do j = 1, size(s%u(i)%data)
    data => s%u(i)%data(j)
    if (data%exists .and. .not. data%good_model) then
      call out_io(s_warn$, r_name, &
                    'DATUM EXISTS BUT CANNOT COMPUTE A MODEL VALUE: ' // tao_datum_name(data))
      cycle
    endif
  enddo
enddo

do i = lbound(s%u, 1), ubound(s%u, 1)
  do j = 1, size(s%u(i)%data)
    data => s%u(i)%data(j)
    if (data%exists .and. .not. data%good_model) call out_io(s_blank$, r_name, &
                                              '         ' // tao_datum_name(data))
  enddo
enddo

! Look for a startup file

call tao_open_file ('TAO_INIT_DIR', startup_file, iu, file_name, .false.)
if (iu /= 0) then
  call out_io (s_blank$, r_name, 'Using startup file: ' // file_name)
  tao_com%cmd_from_cmd_file = .false.
  call tao_cmd_history_record ('call ' // startup_file)
  call tao_call_cmd (file_name)
endif

! Bookkeeping

call tao_set_data_useit_opt()
call tao_set_var_useit_opt()
if (present(err_flag)) err_flag = .false.

contains
!------------------------------------------------------------------------------
! every pointer and allocatable needs to be deallocated now before the universe
! is reallocated.

subroutine deallocate_everything ()

implicit none

type (tao_plot_struct), pointer :: plot
type (tao_curve_struct), pointer :: curve
type (tao_universe_struct), pointer :: u

integer i, j, k, istat

! Tunnel walls

if (allocated(s%wall)) deallocate (s%wall)

! Variables  

if (allocated (s%v1_var)) then
  deallocate(s%v1_var, stat=istat)
endif
  
if (allocated (s%var)) then
  do i = lbound(s%var,1), ubound(s%var,1)
    deallocate(s%var(i)%this, stat=istat)
  enddo
  deallocate(s%var, stat=istat)
endif

! Keytable 

if (allocated(s%key)) deallocate(s%key, stat=istat)

! Plotting  

if (allocated(s%plot_region)) deallocate (s%plot_region)

do i = 1, size(s%template_plot)
  plot => s%template_plot(i)
  if (.not. allocated (plot%graph)) cycle
  deallocate(plot%graph, stat=istat)
enddo
deallocate (s%template_plot)

if (allocated(tao_com%ele_shape_lat_layout)) deallocate (tao_com%ele_shape_lat_layout)
if (allocated(tao_com%ele_shape_floor_plan)) deallocate (tao_com%ele_shape_floor_plan)
if (allocated(tao_com%covar))                deallocate (tao_com%covar, tao_com%alpha)

! Universes 

if (allocated (s%u)) then
  do i = lbound(s%u, 1), ubound(s%u, 1)

    u => s%u(i)
    ! radiation integrals cache
    if (u%ix_rad_int_cache /= 0) call release_rad_int_cache(u%ix_rad_int_cache)

    ! Orbits

    deallocate(u%model%lat_branch, stat=istat)
    deallocate(u%design%lat_branch, stat=istat)
    deallocate(u%base%lat_branch, stat=istat)
    
    deallocate(u%model%bunch_params2, stat=istat)
    deallocate(u%design%bunch_params2, stat=istat)
    deallocate(u%base%bunch_params2, stat=istat)
    
    ! Beams: All s%u(i)%ele point to the same place with common_lattice.

    if (i == 0 .or. .not. tao_com%common_lattice) then
      do ib = 0, ubound(u%uni_branch, 1)
        call reallocate_beam(u%uni_branch(ib)%ele(0)%beam, 0, 0)
        deallocate (u%uni_branch(ib)%ele)
      enddo
      deallocate (u%uni_branch)
    endif

    call reallocate_beam(u%current_beam, 0, 0)

    ! Connected universes
    call deallocate_ele_pointers (u%connect%match_ele)
    call reallocate_beam (u%connect%injecting_beam, 0, 0)
    
    ! d2_data
    do j = 1, size(u%d2_data)
      do k = 1, size(u%d2_data(j)%d1)
        nullify(u%d2_data(j)%d1(k)%d2)
        nullify(u%d2_data(j)%d1(k)%d)
      enddo
      deallocate(u%d2_data(j)%d1, stat=istat)
    enddo
    deallocate(u%d2_data, stat=istat)
 
    ! Data
    do j = lbound(u%data,1), ubound(u%data,1)
      nullify(u%data(j)%d1)
    enddo
    deallocate(u%data, stat=istat)
     
    ! Lattices
    call deallocate_lat_pointers (u%model%lat)
    call deallocate_lat_pointers (u%design%lat)
    call deallocate_lat_pointers (u%base%lat)

  enddo
endif

deallocate (s%u)
    
end subroutine deallocate_everything
    
end subroutine tao_init


