module tao_set_mod

use tao_mod
use tao_data_and_eval_mod
use quick_plot
use tao_lattice_calc_mod

contains

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! Subroutine tao_set_lattice_cmd (dest_lat, source_lat)
!
! Sets a lattice equal to another. This will also update the data structs
! If the 
!
! Input:
!   dest_lat -- Character(*): Maybe: 'model', 'design', or 'base' with 
!                     optional '@n' at beginning to indicate the universe
!   source_lat  -- Character(*): Maybe: 'model', 'design', or 'base' 
!
!  Output:
!    s%u(n) -- lat_struct: changes specified lattice in specified universe 
!-

subroutine tao_set_lattice_cmd (dest_lat, source_lat)

implicit none

character(*) dest_lat, source_lat
character(16) dest1_name
character(20) :: r_name = 'tao_set_lattice_cmd'

integer i, ib

logical, allocatable, save :: this_u(:)
logical err

call tao_pick_universe (dest_lat, dest1_name, this_u, err)
if (err) return

do i = lbound(s%u, 1), ubound(s%u, 1)
  if (.not. this_u(i)) cycle
  call set_lat (s%u(i))
  if (err) return
enddo

!-------------------------------------------
contains

subroutine set_lat (u)

implicit none

type (tao_universe_struct), target :: u
type (tao_lattice_struct), pointer :: dest1_lat
type (tao_lattice_struct), pointer :: source1_lat
real(rp), pointer :: dest_data(:)
real(rp), pointer :: source_data(:)
logical, pointer :: dest_good(:), source_good(:)
logical calc_ok

integer j

!

err = .false.

select case (dest1_name)
  case ('model')
    dest1_lat => u%model
    dest_data => u%data%model_value
    dest_good => u%data%good_model
  case ('base')
    dest1_lat => u%base
    dest_data => u%data%base_value
    dest_good => u%data%good_base
  case default
    call out_io (s_error$, r_name, 'BAD LATTICE: ' // dest_lat)
    err = .true.
    return
end select

select case (source_lat)
  case ('model')
    ! make sure model data is up to date
    tao_com%lattice_recalc = .true.
    call tao_lattice_calc (calc_ok)
    source1_lat => u%model
    source_data => u%data%model_value
    source_good => u%data%good_model
  case ('base')
    source1_lat => u%base
    source_data => u%data%base_value
    source_good => u%data%good_base
  case ('design')
    source1_lat => u%design
    source_data => u%data%design_value
    source_good => u%data%good_design
  case default
    call out_io (s_error$, r_name, 'BAD LATTICE: ' // source_lat)
    err = .true.
    return
end select

dest1_lat%lat          = source1_lat%lat
dest1_lat%lat_branch   = source1_lat%lat_branch
dest1_lat%modes        = source1_lat%modes
dest1_lat%a            = source1_lat%a
dest1_lat%b            = source1_lat%b
dest1_lat%rad_int      = source1_lat%rad_int

do ib = 0, ubound(dest1_lat%lat_branch, 1)
  do j = lbound(dest1_lat%lat_branch(ib)%bunch_params, 1), ubound(dest1_lat%lat_branch(ib)%bunch_params, 1)
    dest1_lat%lat_branch(ib)%bunch_params(j) = source1_lat%lat_branch(ib)%bunch_params(j)
  enddo
enddo

if (allocated(source1_lat%bunch_params2)) then
  if (allocated(dest1_lat%bunch_params2)) then
    if (size(dest1_lat%bunch_params2) < size(source1_lat%bunch_params2)) &
                                        deallocate (dest1_lat%bunch_params2)
  endif
  if (.not. allocated(dest1_lat%bunch_params2)) &
                        allocate (dest1_lat%bunch_params2(source1_lat%n_bunch_params2))
  dest1_lat%n_bunch_params2 = source1_lat%n_bunch_params2
else
  if (allocated(dest1_lat%bunch_params2)) then
    deallocate(dest1_lat%bunch_params2)
    dest1_lat%n_bunch_params2 = 0
  endif
endif

! Transfer the data

dest_data = source_data
dest_good = source_good

end subroutine set_lat

end subroutine tao_set_lattice_cmd

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_set_global_cmd (who, set_value)
!
! Routine to set global variables
! 
! Input:
!   who       -- Character(*): which global variable to set
!   set_value -- Character(*): Value to set to.
!
! Output:
!    s%global  -- Global variables structure.
!-

subroutine tao_set_global_cmd (who, set_value)

implicit none

type (tao_global_struct) global

character(*) who, set_value
character(20) :: r_name = 'tao_set_global_cmd'

integer iu, ios
logical err

namelist / params / global

! open a scratch file for a namelist read

iu = lunget()
open (iu, status = 'scratch', iostat = ios)
if (ios /= 0) then
  call out_io (s_error$, r_name, 'CANNOT OPEN A SCRATCH FILE!')
  return
endif

write (iu, *) '&params'
write (iu, *) ' global%' // trim(who) // ' = ' // trim(set_value)
write (iu, *) '/'
rewind (iu)
global = s%global  ! set defaults
read (iu, nml = params, iostat = ios)
close (iu)

call tao_data_check (err)
if (err) return

if (ios == 0) then
  s%global = global
  if (trim(who) == 'track_type') tao_com%lattice_recalc = .true.
else
  call out_io (s_error$, r_name, 'BAD COMPONENT OR NUMBER')
endif

end subroutine tao_set_global_cmd

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_set_wave_cmd (who, set_value, err)
!
! Routine to set wave variables
! 
! Input:
!   who       -- Character(*): which wave variable to set
!   set_value -- Character(*): Value to set to.
!
! Output:
!    s%wave  -- Wave variables structure.
!-

subroutine tao_set_wave_cmd (who, set_value, err)

implicit none

type (tao_wave_struct) wave

character(*) who, set_value
character(20) :: r_name = 'tao_set_wave_cmd'

real(rp) ix_a(2), ix_b(2)

integer iu, ios
logical err

namelist / params / ix_a, ix_b

! open a scratch file for a namelist read

err = .true.

iu = lunget()
open (iu, status = 'scratch', iostat = ios)
if (ios /= 0) then
  call out_io (s_error$, r_name, 'CANNOT OPEN A SCRATCH FILE!')
  return
endif

ix_a = [s%wave%ix_a1, s%wave%ix_a2]
ix_b = [s%wave%ix_b1, s%wave%ix_b2]

write (iu, *) '&params'
write (iu, *) trim(who) // ' = ' // trim(set_value)
write (iu, *) '/'
rewind (iu)
wave = s%wave  ! set defaults
read (iu, nml = params, iostat = ios)
close (iu)

if (ios /= 0) then
  call out_io (s_error$, r_name, 'BAD COMPONENT OR NUMBER')
  return
endif

s%wave%ix_a1 = ix_a(1)
s%wave%ix_a2 = ix_a(2)
s%wave%ix_b1 = ix_b(1)
s%wave%ix_b2 = ix_b(2)

err = .false.

end subroutine tao_set_wave_cmd

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_set_beam_init_cmd (who, set_value)
!
! Routine to set beam_init variables
! 
! Input:
!   who       -- Character(*): which beam_init variable to set
!   set_value -- Character(*): Value to set to.
!
! Output:
!    s%beam_init  -- Beam_init variables structure.
!-

subroutine tao_set_beam_init_cmd (who, set_value, ix_branch)

implicit none

type (beam_init_struct) beam_init
type (tao_universe_struct), pointer :: u
character(*) who, set_value
character(40) who2
character(20) :: r_name = 'tao_set_beam_init_cmd'

integer ix_branch
integer i, iu, ios
logical err
logical, allocatable :: picked_uni(:)

namelist / params / beam_init

! get universe

call tao_pick_universe (who, who2, picked_uni, err)

! open a scratch file for a namelist read

iu = lunget()
do i = lbound(s%u, 1), ubound(s%u, 1)

  open (iu, status = 'scratch', iostat = ios)
  if (ios /= 0) then
    call out_io (s_error$, r_name, 'CANNOT OPEN A SCRATCH FILE!')
    return
  endif

  write (iu, *) '&params'
  write (iu, *) ' beam_init%' // trim(who2) // ' = ' // trim(set_value)
  write (iu, *) '/'
  rewind (iu)
  u => s%u(i)
  beam_init = u%uni_branch(ix_branch)%beam_init  ! set defaults
  read (iu, nml = params, iostat = ios)
  close (iu)

  call tao_data_check (err)
  if (err) return

  if (ios == 0) then
    u%uni_branch(ix_branch)%beam_init = beam_init
    tao_com%lattice_recalc = .true.
  else
    call out_io (s_error$, r_name, 'BAD COMPONENT OR NUMBER')
    return
  endif

enddo

end subroutine tao_set_beam_init_cmd

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Subroutine tao_set_plot_page_cmd (component, set_value)
!
!  Set various aspects of the plotting window
!
! Input:
!   component     -- Character(*): Which component to set.
!   set_value     -- Character(*): What value to set to.
!
!  Output:
!    s%plot_page  -- tao_plot_page_struct:
!-

subroutine tao_set_plot_page_cmd (component, set_value1, set_value2)

implicit none

character(*) component, set_value1
character(*), optional :: set_value2
character(24) :: r_name = 'tao_set_plot_page_cmd'

real(rp) x, y
integer ix
logical error

select case (component)

case ('title')
  s%plot_page%title(1)%string = trim(set_value1)

case ('subtitle')
  s%plot_page%title(2)%string = trim(set_value1)
  s%plot_page%title(2)%draw_it = .true.

case ('subtitle_loc')

  if (.not. present(set_value2)) then
    call out_io(s_info$, r_name, "subtitle_loc requires two numbers.")
    return
  endif

  read(set_value1, '(f15.10)') x
  read(set_value2, '(f15.10)') y
  s%plot_page%title(2)%x = x
  s%plot_page%title(2)%y = y

case ('shape_height_max')
  call tao_real_set_value (s%plot_page%shape_height_max, component, set_value1, error)

case ('text_height')
  call tao_real_set_value (s%plot_page%text_height, component, set_value1, error)

case ('main_title_text_scale')
  call tao_real_set_value (s%plot_page%main_title_text_scale, component, set_value1, error)

case ('graph_title_text_scale')
  call tao_real_set_value (s%plot_page%graph_title_text_scale, component, set_value1, error)

case ('axis_number_text_scale')
  call tao_real_set_value (s%plot_page%axis_number_text_scale, component, set_value1, error)

case ('axis_label_text_scale')
  call tao_real_set_value (s%plot_page%axis_label_text_scale, component, set_value1, error)

case ('legend_text_scale')
  call tao_real_set_value (s%plot_page%legend_text_scale, component, set_value1, error)

case ('key_table_text_scale')
  call tao_real_set_value (s%plot_page%key_table_text_scale, component, set_value1, error)


case default
  call out_io (s_error$, r_name, 'PLOT COMPONENT NOT RECOGNIZED: ' // component)

end select

end subroutine tao_set_plot_page_cmd

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_set_curve_cmd (curve_name, component, set_value)
!
! Routine to set var values.
!
! Input:
!   curve_name -- Character(*): Which curve to set.
!   component  -- Character(*): Which component to set.
!   set_value  -- Character(*): What value to set it to.
!-

subroutine tao_set_curve_cmd (curve_name, component, set_value)

implicit none

type (tao_curve_array_struct), allocatable, save :: curve(:)
type (tao_graph_array_struct), allocatable, save :: graph(:)
type (lat_struct), pointer :: lat

integer i, j, ios, i_uni
integer, allocatable, save :: ix_ele(:)

character(*) curve_name, component, set_value
character(20) :: r_name = 'tao_set_curve_cmd'

logical err

!

call tao_find_plots (err, curve_name, 'REGION', curve = curve, always_allocate = .true.)
if (err) return

if (.not. allocated(curve) .or. size(curve) == 0) then
  call out_io (s_error$, r_name, 'CURVE OR GRAPH NOT SPECIFIED')
  return
else
  do i = 1, size(curve)
    call set_this_curve (curve(i)%c)
  enddo
endif

!---------------------------------------------
contains

subroutine set_this_curve (this_curve)

type (tao_curve_struct) this_curve
type (tao_graph_struct), pointer :: this_graph
type (tao_universe_struct), pointer :: u
type (tao_universe_branch_struct), pointer :: uni_branch
type (ele_pointer_struct), allocatable :: eles(:)

integer ix, ib, i_branch
logical error

!

i_branch = this_curve%ix_branch
i_uni = tao_universe_number(this_curve%ix_universe)

this_graph => this_curve%g

! if the universe is changed then need to check ele_ref

select case (component)

case ('ele_ref_name')
  this_curve%ele_ref_name = set_value
  call tao_locate_elements (this_curve%ele_ref_name, i_uni, eles, error, .true.)
  if (size(eles) == 0) return
  this_curve%ix_ele_ref = eles(1)%ele%ix_ele
  this_curve%ix_branch  = eles(1)%ele%ix_branch
  call tao_ele_to_ele_track (i_uni, i_branch, this_curve%ix_ele_ref, this_curve%ix_ele_ref_track)
  
case ('ix_ele_ref')
  call tao_integer_set_value (this_curve%ix_ele_ref, component, &
                    set_value, error, 0, s%u(i_uni)%model%lat%branch(i_branch)%n_ele_max)
  this_curve%ele_ref_name = s%u(i_uni)%model%lat%ele(this_curve%ix_ele_ref)%name
  call tao_ele_to_ele_track (this_curve%ix_universe, ib, &
                                this_curve%ix_ele_ref, this_curve%ix_ele_ref_track)

case ('ix_universe')
  call tao_integer_set_value (this_curve%ix_universe, component, &
                                            set_value, error, 0, ubound(s%u, 1))
  if (error) return
  call tao_locate_elements (this_curve%ele_ref_name, this_curve%ix_universe, eles, error, .true.)
  if (size(eles) == 0) return
  this_curve%ix_ele_ref = eles(1)%ele%ix_ele
  this_curve%ix_branch  = eles(1)%ele%ix_branch
  call tao_ele_to_ele_track (this_curve%ix_universe, this_curve%ix_branch, &
                                     this_curve%ix_ele_ref, this_curve%ix_ele_ref_track)

case ('ix_branch') 
  call tao_integer_set_value (this_curve%ix_branch, component, set_value, error)

case ('ix_bunch')
  u => tao_pointer_to_universe (this_curve%ix_universe)
  if (.not. associated(u)) return
  call tao_integer_set_value (this_curve%ix_bunch, component, &
                        set_value, error, -1, u%uni_branch(i_branch)%beam_init%n_bunch)

case ('symbol_every')
  call tao_integer_set_value (this_curve%symbol_every, component, &
                                            set_value, error, 0, size(this_curve%x_symb))

case ('draw_line')
  call tao_logical_set_value (this_curve%draw_line, component, set_value, error)

case ('draw_symbols')
  call tao_logical_set_value (this_curve%draw_symbols, component, set_value, error)

case ('draw_symbol_index')
  call tao_logical_set_value (this_curve%draw_symbol_index, component, set_value, error)

case default
  call out_io (s_error$, r_name, "BAD CURVE COMPONENT")
  return
    
end select

! Set ix_ele_ref_track if necessary

select case (component)
case ('ele_ref_name', 'ix_ele_ref', 'ix_universe')

end select

! Enable

if (this_graph%type == 'phase_space') then
  uni_branch => s%u(i_uni)%uni_branch(i_branch)
  if (.not. uni_branch%ele(this_curve%ix_ele_ref)%save_beam) then
    tao_com%lattice_recalc = .true.
    uni_branch%ele(this_curve%ix_ele_ref)%save_beam = .true.
  endif
endif

end subroutine

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_set_graph_cmd (graph_name, component, set_value)
!
! Routine to set var values.
!
! Input:
!   graph_name -- Character(*): Which graph to set.
!   component  -- Character(*): Which component to set.
!   set_value  -- Character(*): What value to set it to.
!
!  Output:
!-

subroutine tao_set_graph_cmd (graph_name, component, set_value)

implicit none

type (tao_plot_array_struct), allocatable, save :: plot(:)
type (tao_graph_array_struct), allocatable, save :: graph(:)

character(*) graph_name, component, set_value
character(20) :: r_name = 'tao_set_graph_cmd'

integer i, j, ios, i_uni
logical err

!

call tao_find_plots (err, graph_name, 'REGION', plot = plot, graph = graph)
if (err) return

if (allocated(graph)) then
  do i = 1, size(graph)
    call set_this_graph (graph(i)%g)
  enddo
elseif (allocated(plot)) then
  do i = 1, size(plot)
    do j = 1, size(plot(i)%p%graph)
      call set_this_graph (plot(i)%p%graph(j))
    enddo
  enddo
else
  call out_io (s_error$, r_name, 'GRAPH OR PLOT NOT SPECIFIED')
  return
endif

!---------------------------------------------
contains

subroutine set_this_graph (this_graph)

type (tao_graph_struct) this_graph
character(40) comp
integer iset, iw, ix
logical logic, error

!

i_uni = this_graph%ix_universe
if (i_uni == 0) i_uni = s%global%u_view
comp = component

select case (comp)

  case ('component')
    this_graph%component = set_value

  case ('clip')
    read (set_value, '(l)', iostat = ios) logic
    if (ios /= 0) then
      call out_io (s_error$, r_name, 'BAD CLIP VALUE.')
      return
    endif
    this_graph%clip = logic

  case ('draw_axes')
    read (set_value, '(l)', iostat = ios) logic
    if (ios /= 0) then
      call out_io (s_error$, r_name, 'BAD CLIP VALUE.')
      return
    endif
    this_graph%draw_axes = logic

  case ('ix_universe')
    call tao_integer_set_value (this_graph%ix_universe, comp, set_value, error, 1, ubound(s%u, 1))

  case ('margin%x1')
    call tao_real_set_value(this_graph%margin%x1, comp, set_value, error)

  case ('margin%x2')
    call tao_real_set_value(this_graph%margin%x2, comp, set_value, error)

  case ('margin%y1')
    call tao_real_set_value(this_graph%margin%y1, comp, set_value, error)

  case ('margin%y2')
    call tao_real_set_value(this_graph%margin%y2, comp, set_value, error)

  case default
    call out_io (s_error$, r_name, "BAD GRAPH COMPONENT: " // component)
    return
    
end select

tao_com%lattice_recalc = .true.

end subroutine
end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_set_var_cmd (var_str, value_str)
!
! Routine to set var values.
!
! Input:
!   var_str  -- Character(*): Which var name to set.
!   value_str  -- Character(*): What value to set it to.
!
!  Output:
!-

subroutine tao_set_var_cmd (var_str, value_str)

implicit none

type (tao_v1_var_struct), pointer :: v1_ptr
type (tao_real_pointer_struct), allocatable, save    :: r_var(:), r_set(:)
type (tao_logical_array_struct), allocatable, save :: l_var(:), l_set(:)
type (tao_var_array_struct), allocatable, save     :: v_var(:)
type (tao_string_array_struct), allocatable, save :: s_var(:), s_set(:)

real(rp), allocatable, save :: r_value(:)
real(rp) value
integer i, j

character(*) var_str, value_str
character(20) :: r_name = 'tao_set_var_cmd'
character(20) set_is, component
character(40) :: merit_type_names(2) = (/ 'target ', 'limit  ' /)

logical err, l_value, err_flag
logical, allocatable, save :: good(:)

! Decode variable component to set.

call tao_find_var (err, var_str, v_array = v_var, re_array=r_var, &
                   log_array=l_var, str_array = s_var, component = component)

! A logical value_str is either a logical or an array of datum values.

if (size(l_var) > 0) then
  if (is_logical(value_str)) then
    read (value_str, *) l_value
    do i = 1, size(l_var)
      l_var(i)%l = l_value
    enddo

  else
    call tao_find_var (err, value_str, log_array=l_set)
    if (size(l_set) /= size(l_var)) then
      call out_io (s_error$, r_name, 'ARRAY SIZES ARE NOT THE SAME')
      return
    endif
    do i = 1, size(l_var)
      l_var(i)%l = l_set(i)%l
    enddo
  endif

! Must be merit_type for a string.
! If value_string has "|" then it must be a datum array

elseif (size(s_var) /= 0) then
  if (index(value_str, '|') == 0) then
    if (all (value_str /= merit_type_names)) then
      call out_io (s_error$, r_name, 'BAD MERIT_TYPE NAME:' // value_str)
      return
    endif
    do i = 1, size(s_var)
      s_var(i)%s = value_str
    enddo

  else
    call tao_find_var (err, value_str, str_array=s_set)
    if (size(l_set) /= size(l_var)) then
      call out_io (s_error$, r_name, 'ARRAY SIZES ARE NOT THE SAME')
      return
    endif
    do i = 1, size(s_var)
      s_var(i)%s = s_set(i)%s
    enddo
  endif

! Only possibility left is real/ The value_str might be a number or it might 
! be a mathematical expression involving datum values or array of values.

elseif (size(r_var) /= 0) then
  call tao_evaluate_expression (value_str, size(r_var),  .false., r_value, good, err)
  if (err) then
    call out_io (s_error$, r_name, 'BAD SET VALUE ' // value_str)
    return
  endif

  do i = 1, size(r_var)
    if (component == 'model') then
      call tao_set_var_model_value (v_var(i)%v, r_value(i))
    else
      r_var(i)%r = r_value(i)
    endif
  enddo

! Else must be an error

else
  call out_io (s_error$, r_name, 'NOTHING TO SET!')
  return

endif

call tao_set_var_useit_opt()

end subroutine tao_set_var_cmd

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_set_data_cmd (who_str, value_str)
!
! Routine to set data values.
!
! Input:
!   who_str   -- Character(*): Which data component(s) to set.
!   value_str -- Character(*): What value to set it to.
!
!  Output:
!-

subroutine tao_set_data_cmd (who_str, value_str)

implicit none

type (tao_real_pointer_struct), allocatable, save    :: r_dat(:), r_set(:)
type (tao_data_array_struct), allocatable, save    :: d_dat(:)
type (tao_integer_array_struct), allocatable, save :: int_dat(:), int_set(:)
type (tao_logical_array_struct), allocatable, save :: l_dat(:), l_set(:)
type (tao_string_array_struct), allocatable, save :: s_dat(:), s_set(:)
type (tao_universe_struct), pointer :: u
type (branch_struct), pointer :: branch

real(rp), allocatable, save :: r_value(:)
integer i, ix, int_value

integer, allocatable :: i_save(:)

character(*) who_str, value_str
character(20) component
character(20) :: r_name = 'tao_set_data_cmd'
character(40) :: merit_type_names(5) = &
              (/ 'target ', 'min    ', 'max    ', 'abs_min', 'abs_max' /)
logical err, l_value, valid_value
logical, allocatable :: good(:)

! Decode data component to set.

call tao_find_data (err, who_str, d_array = d_dat, re_array=r_dat, &
          log_array=l_dat, str_array = s_dat, int_array = int_dat, component = component)
if (err) return

! A logical value_str is either a logical or an array of datum values.

if (size(l_dat) /= 0) then
  if (is_logical(value_str)) then
    read (value_str, *) l_value
    do i = 1, size(l_dat)
      l_dat(i)%l = l_value
    enddo

  else
    call tao_find_data (err, value_str, log_array=l_set)
    if (size(l_set) /= size(l_dat)) then
      call out_io (s_error$, r_name, 'ARRAY SIZES ARE NOT THE SAME')
      return
    endif
    do i = 1, size(l_dat)
      l_dat(i)%l = l_set(i)%l
    enddo
  endif

! An integer value_str is either an integer or an array of datum values.

elseif (size(int_dat) /= 0) then

  allocate (i_save(size(int_dat)))

  if (is_integer(value_str)) then
    read (value_str, *) int_value
    do i = 1, size(int_dat)
      i_save(i) = int_dat(i)%i
      int_dat(i)%i = int_value
    enddo

  else
    call tao_find_data (err, value_str, int_array=int_set)
    if (size(int_set) /= size(int_dat)) then
      call out_io (s_error$, r_name, 'ARRAY SIZES ARE NOT THE SAME')
      return
    endif
    do i = 1, size(int_dat)
      i_save(i) = int_dat(i)%i
      int_dat(i)%i = int_set(i)%i
    enddo
  endif

  if (component == 'ix_ele' .or. component == 'ix_ele_start' .or. component == 'ix_ele_ref') then
    do i = 1, size(int_dat)
      u => s%u(d_dat(i)%d%d1%d2%ix_uni)
      branch => u%design%lat%branch(d_dat(i)%d%ix_branch)
      if (int_dat(i)%i < 0 .or. int_dat(i)%i > branch%n_ele_max) then
        int_dat(i)%i = i_save(i)
      else
        if (component == 'ix_ele') then
          d_dat(i)%d%ele_name = branch%ele(int_dat(i)%i)%name
        elseif (component == 'ix_ele_start') then
          d_dat(i)%d%ele_start_name = branch%ele(int_dat(i)%i)%name
        else
          d_dat(i)%d%ele_ref_name = branch%ele(int_dat(i)%i)%name
        endif
      endif
    enddo
  endif

  deallocate (i_save)

! A string:

elseif (size(s_dat) /= 0) then

  ! If value_string has "|" then it must be a datum array

  if (index(value_str, '|') == 0) then
    if (all (component == 'merit_type' .and. value_str /= merit_type_names)) then
      call out_io (s_error$, r_name, 'BAD MERIT_TYPE NAME:' // value_str)
      return
    endif
    do i = 1, size(s_dat)
      s_dat(i)%s = value_str
    enddo

  else
    call tao_find_data (err, value_str, str_array=s_set)
    if (size(s_set) /= size(s_dat) .and. size(s_set) /= 1) then
      call out_io (s_error$, r_name, 'ARRAY SIZES ARE NOT THE SAME')
      return
    endif

    do i = 1, size(s_dat)
      if (size(s_set) == 1) then
        s_dat(i)%s = s_set(1)%s
      else
        s_dat(i)%s = s_set(i)%s
      endif
    enddo

  endif

  if (component == 'ele_name' .or. component == 'ele_start_name' .or. component == 'ele_ref_name') then
    do i = 1, size(d_dat)
      u => s%u(d_dat(i)%d%d1%d2%ix_uni)
      call element_locator (s_dat(i)%s, u%design%lat, ix)
      if (ix < 0) then
        call out_io (s_error$, r_name, 'ELEMENT NOT LOCATED: ' // s_dat(i)%s)
        cycle
      endif
      if (component == 'ele_name') then
        d_dat(i)%d%ix_ele = ix
      elseif (component == 'ele_start_name') then
        d_dat(i)%d%ix_ele_start = ix
      else
        d_dat(i)%d%ix_ele_ref = ix
      endif
    enddo
  endif

! Only possibility left is real. The value_str might be a number or it might 
! be a mathematical expression involving datum values or array of values.

elseif (size(r_dat) /= 0) then
  call tao_evaluate_expression (value_str, size(r_dat), .false., r_value, good, err)
  if (err) then
    call out_io (s_error$, r_name, 'BAD SET VALUE ' // value_str)
    return
  endif

  do i = 1, size(r_dat)
    r_dat(i)%r = r_value(i)
    if (component == 'meas') d_dat(i)%d%good_meas = .true.
    if (component == 'ref')  d_dat(i)%d%good_ref = .true.
  enddo

endif

do i = 1, size(d_dat)
  u => s%u(d_dat(i)%d%d1%d2%ix_uni)  
  call tao_evaluate_a_datum (d_dat(i)%d, u, u%model, d_dat(i)%d%model_value, valid_value)
enddo

call tao_set_data_useit_opt()

end subroutine tao_set_data_cmd

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_set_uni_cmd (uni, who, what)
!
! Sets a universe on or off, or sets the recalculate or mat6_recalc logicals
!
! Input:
!   uni     -- Character(*): which universe; 0 => current viewed universe
!   who     -- Character(*): "on", "off", "recalculate", or "mat6_recalc"
!   what    -- Character(*): "on" or "off" for who = "mat6_recalc".
!
!-

subroutine tao_set_uni_cmd (uni, who, what)

implicit none

integer i, n_uni

character(*) uni, who, what
character(20) :: r_name = "tao_set_universe_cmd"

logical is_on, err, recalc, mat6_toggle


! Pick universe

if (uni /= '*') then
  call tao_to_int (uni, n_uni, err)
  if (err) return
  if (n_uni < -1 .or. n_uni > ubound(s%u, 1)) then
    call out_io (s_warn$, r_name, "Invalid Universe specifier")
    return 
  endif
  n_uni = tao_universe_number (n_uni)
endif

!

if (index('mat6_recalc', trim(who)) == 1) then
  if (what == 'on') then
    is_on = .true.
  elseif (what == 'off') then
    is_on = .false.
  else
    call out_io (s_error$, r_name, 'Syntax is: "set universe <uni_num> mat6_recalc on/off"')
    return
  endif
  if (uni == '*') then
    s%u(:)%mat6_recalc_on = is_on
  else
    s%u(n_uni)%mat6_recalc_on = is_on
  endif
  tao_com%lattice_recalc = .true.
  return
endif
  
!

if (what /= '') then
  call out_io (s_error$, r_name, 'Extra stuff on line. Nothing done.')
  return
endif


recalc = .false.

if (who == 'on') then
  is_on = .true.
  recalc = .true.
elseif (who == 'off') then
  is_on = .false.
elseif (index('recalculate', trim(who)) == 1) then
  recalc = .true.
else
  call out_io (s_error$, r_name, "Choices are: 'on', 'off', 'recalculate', or 'mat6_recalc")
  return
endif

!

if (uni == '*') then
  call out_io (s_blank$, r_name, "Setting all universes to: " // on_off_logic(is_on))
  s%u(:)%is_on = is_on
else
  s%u(n_uni)%is_on = is_on
  call out_io (s_blank$, r_name, "Setting universe \i0\ to: " // on_off_logic(is_on), n_uni)
endif

! make sure lattice calculation is up to date if turning lattice on

call tao_set_data_useit_opt()
if (recalc) tao_com%lattice_recalc = .true.
  
end subroutine tao_set_uni_cmd

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_logical_set_value (var, var_str, value_str, error)
!
! Subroutine to read and set the value of an logical varialbe.
!
! If the value is out of the range [min_val, max_val] then an error message will
! be generated and the variable will not be set.
!
! Input:
!   var_str   -- Character(*): Used for error messages.
!   value_str -- Character(*): String with encoded value.
!
! Output:
!   var   -- Logical: Variable to set.
!   error -- Logical: Set True on an error. False otherwise.
!-

subroutine tao_logical_set_value (var, var_str, value_str, error)

implicit none

logical var, ix
integer ios

character(*) var_str, value_str
character(24) :: r_name = 'tao_logical_set_value'
logical error

!

error = .true.
read (value_str, '(l)', iostat = ios) ix

if (ios /= 0 .or. len_trim(value_str) == 0) then
  call out_io (s_error$, r_name, 'BAD ' // trim(var_str) // ' VALUE.')
  return
endif

var = ix      
error = .false.

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_integer_set_value (var, var_str, value_str, error, min_val, max_val)
!
! Subroutine to read and set the value of an integer varialbe.
!
! If the value is out of the range [min_val, max_val] then an error message will
! be generated and the variable will not be set.
!
! Input:
!   var_str   -- Character(*): Used for error messages.
!   value_str -- Character(*): String with encoded value.
!   min_val   -- Integer, optional: Minimum value. 
!   max_val   -- Integer, optional: Maximum value.
!
! Output:
!   var   -- Integer: Variable to set.
!   error -- Logical: Set True on an error. False otherwise.
!-

subroutine tao_integer_set_value (var, var_str, value_str, error, min_val, max_val)

implicit none

integer var
integer, optional :: min_val, max_val
integer ios, ix

character(*) var_str, value_str
character(20) :: r_name = 'tao_integer_set_value'
logical error

!

error = .true.
read (value_str, *, iostat = ios) ix

if (ios /= 0 .or. len_trim(value_str) == 0) then
  call out_io (s_error$, r_name, 'BAD ' // trim(var_str) // ' VALUE.')
  return
endif

if (ix < min_val .or. ix > max_val) then
  call out_io (s_error$, r_name, var_str // ' VALUE OUT OF RANGE.')
  return
endif

var = ix      
error = .false.

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_real_set_value (var, var_str, value_str, error, min_val, max_val)
!
! Subroutine to read and set the value of a real varialbe.
!
! If the value is out of the range [min_val, max_val] then an error message will
! be generated and the variable will not be set.
!
! Input:
!   var_str   -- Character(*): Used for error messages.
!   value_str -- Character(*): String with encoded value.
!   min_val   -- real(rp), optional: Minimum value. 
!   max_val   -- real(rp), optional: Maximum value.
!
! Output:
!   var   -- real(rp): Variable to set.
!   error -- Logical: Set True on an error. False otherwise.
!-

subroutine tao_real_set_value (var, var_str, value_str, error, min_val, max_val)

implicit none

real(rp) var, var_value
real(rp), optional :: min_val, max_val
integer ios

character(*) var_str, value_str
character(20) :: r_name = 'tao_real_set_value'
logical error

!

error = .true.
read (value_str, *, iostat = ios) var_value

if (ios /= 0 .or. len_trim(var_str) == 0) then
  call out_io (s_error$, r_name, 'BAD ' // trim(var_str) // ' VALUE.')
  return
endif

if (var_value < min_val .or. var_value > max_val) then
  call out_io (s_error$, r_name, var_str // ' VALUE OUT OF RANGE.')
  return
endif

var = var_value
error = .false.

end subroutine

end module tao_set_mod
