!+
! Module tao_utils
!
! helper subroutines available for communal usage.
!-

module tao_utils

use tao_struct
use tao_interface
use bmad

! used for parsing expressions
integer, private :: plus$ = 1, minus$ = 2, times$ = 3, divide$ = 4
integer, private :: l_parens$ = 5, r_parens$ = 6, power$ = 7
integer, private :: unary_minus$ = 8, unary_plus$ = 9, no_delim$ = 10
integer, private :: sin$ = 11, cos$ = 12, tan$ = 13
integer, private :: asin$ = 14, acos$ = 15, atan$ = 16, abs$ = 17, sqrt$ = 18
integer, private :: log$ = 19, exp$ = 20, ran$ = 21, ran_gauss$ = 22
integer, private :: numeric$ = 100

integer, private :: eval_level(22) = (/ 1, 1, 2, 2, 0, 0, 4, 3, 3, -1, &
                            9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9 /)

type eval_stack_struct
  integer type
  real(rp) value
end type

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_pick_universe (data_type_in, data_type_out, picked, err)
!
! Subroutine to pick what universe the data name is comming from.
! If data_type_in ends in ";*" or ";0" choose all universes.
! If data_type_in ends in ";n" then choose universe n.
! If not then choose universe s%global%u_view.
! data_type_out is data_type_in without any ";n"
!
! Input:
!   data_type_in -- Character(*): data name.
!
! Output:
!   data_type_out -- Character(*): data_type_in without any ";n" ending.
!   picked(:)     -- Logica: Array showing picked universes.
!   err           -- Logical: Set True if an error is detected.
!-

subroutine tao_pick_universe (data_type_in, data_type_out, picked, err)

implicit none

character(*) data_type_in, data_type_out
character(20) :: r_name = 'tao_pick_universe'
character(1) char

integer ix, n

logical picked(:)
logical err

! Init

err = .false.
picked = .false.

! No ";" then simply choose s%global%u_view

ix = index (data_type_in, ';')
if (ix == 0) then
  picked (s%global%u_view) = .true.
  data_type_out = data_type_in
  return
endif

! Here whn ";" is found...

if (data_type_in(ix+2:) /= ' ') then
  call out_io (s_error$, r_name, 'BAD DATA NAME: ' // data_type_in)
  err = .true.
  return
endif

char = data_type_in(ix+1:ix+1)

if (char == '*' .or. char == '0') then
  picked = .true.
  data_type_out = data_type_in(:ix-1)
  return
endif

n = index ('123456789', char)
if (n /= 0) then
  picked(n) = .true.
  data_type_out = data_type_in(:ix-1)
  return
endif

err = .true.
call out_io (s_error$, r_name, 'BAD UNIVERSE ENCODING IN DATA NAME: ' // data_type_in)

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_point_v1_to_var (v1, var, n, n_var)
!
! used for arbitrary variable pointer indexing
!
! v1       -- tao_v1_var_struct: Contains the pointer
! var(n:)  -- tao_var_struct: the variable
! n        -- integer: starting index for the var array.
!-

subroutine tao_point_v1_to_var (v1, var, n, n_var)

implicit none

integer n, i, n_var

type (tao_v1_var_struct), target :: v1
type (tao_var_struct), target :: var(n:)

v1%v => var

do i = lbound(var, 1), ubound(var, 1)
  var(i)%ix_v1 = i
  var(i)%ix_var = n_var + i - n
  var(i)%v1 => v1
enddo

end subroutine 

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_point_d1_to_data (ip, ii, n, n_data)
!
! Routine used for arbitrary data pointer indexing
!
! ip     -- tao_data_struct: the pointer
! ii     -- tao_data_struct: the data
! n      -- integer: starting index for the pointer
! n_data -- integer: starting index for the next data point in the big data array
!-

subroutine tao_point_d1_to_data (ip, ii, n, n_data)

implicit none

integer n, i, n0, n_data

type (tao_data_struct), pointer :: ip(:)
type (tao_data_struct), target :: ii(n:)

ip => ii

forall (i = lbound(ii, 1):ubound(ii, 1)) 
  ii(i)%ix_d1 = i
  ii(i)%ix_data = n_data + i - n
end forall

end subroutine 

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_locate_element (string, ix_universe, ix_ele, ignore_blank) 
!
! Subroutine to find a lattice element.
!
! Input:
!   string       -- Character(*): String with element name or index
!   ix_universe  -- Integer: Universe to search. 0 => search s%global%u_view.
!   ignore_blank -- Logical, optional: If present and true then do nothing if
!     string is blank. otherwise treated as an error.
!
! Output:
!   ix_ele  -- Integer: Index of element. Set to -1 if element not found.
!-

subroutine tao_locate_element (string, ix_universe, ix_ele, ignore_blank)

implicit none

integer ix_ele, ios, ix, ix_universe

character(*) string
character(16) ele_name
character(20) :: r_name = 'tao_locate_element'

logical, optional :: ignore_blank

! If it is a number translate it:

call str_upcase (ele_name, string)
call string_trim (ele_name, ele_name, ix)

if (ix == 0 .and. logic_option(.false., ignore_blank)) return

if (ix == 0) then
  ix_ele = -1
  call out_io (s_error$, r_name, 'ELEMENT NAME IS BLANK')
  return
endif

ix = ix_universe
if (ix == 0) ix = s%global%u_view

read (ele_name, *, iostat = ios) ix_ele
if (ios .eq. 0) then
  !it's a number
  if (ix_ele < 0 .or. ix_ele > s%u(ix)%model%n_ele_max) then
    ix_ele = -1
    call out_io (s_error$, r_name, 'ELEMENT INDEX OUT OF RANGE: ' // ele_name)
  endif
  return
endif

call element_locator (ele_name, s%u(ix)%model, ix_ele)

if (ix_ele < 0) call out_io (s_error$, r_name, 'ELEMENT NOT FOUND: ' // string)

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_pointer_to_var_in_lattice (var, this, ix_uni, is_ele, err)
! 
! Routine to set a pointer to the appropriate variable in a lattice
!
! Input:
!   var    -- Tao_var_struct: Structure has the info of where to point.
!   this   -- Tao_this_var_struct: the variables pointers to point with
!   ix_uni -- Integer: the universe to use
!   ix_ele -- (Optional) Integer: Point to this element
!
! Output:
!   err   -- Logical: Set True if there is an error. False otherwise.
!-

subroutine tao_pointer_to_var_in_lattice (var, this, ix_uni, ix_ele, err)

implicit none

type (tao_var_struct) var
type (tao_universe_struct), pointer :: u
type (tao_this_var_struct) this

integer, optional :: ix_ele
integer ix, ie, ix_uni
logical, optional :: err
logical error
character(30) :: r_name = 'tao_pointer_to_var_in_lattice'

! locate element

  if (present(err)) err = .true.

  u => s%u(ix_uni)
  if (present(ix_ele)) then
    ie = ix_ele
  else
    call element_locator (var%ele_name, u%model, ie)
  endif

  if (ie < 0) then
    call out_io (s_error$, r_name, 'ELEMENT NAME NOT FOUND: ' // var%ele_name)
    if (present(err)) return
    call err_exit
  endif

  ! locate attribute

  call pointer_to_attribute (u%model%ele_(ie), var%attrib_name, .true., this%model_ptr, ix, error)
  call pointer_to_attribute (u%base%ele_(ie),  var%attrib_name, .true., this%base_ptr,  ix, error)
  if (error) then
    if (present(err)) return
    call err_exit
  endif

  if (present(err)) err = .false.

  this%ix_ele = ie
  this%ix_uni = ix_uni

end subroutine tao_pointer_to_var_in_lattice

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_find_plot_by_region (err, where, plot, graph, curve, region, print_flag)
!
! Routine to find a plot using the region name.
! Optionally find a graph of the plot.
! A region name is something like: where = "top"
! A graph name is something like: where  = "top:x"
! A curve name is something like: where  = "top:x:1"

! Input:
!   where      -- Character(*): Name to match to.
!   print_flag -- Logical, optional: If present and False then surpress error
!                   messages. Default is True.
!
! Output:
!   err      -- logical: Set True on error. False otherwise.
!   plot     -- Tao_plot_struct, pointer, optional: Pointer to the appropriate plot.
!   graph    -- Tao_graph_struct, pointer, optional: Pointer to the appropriate graph.
!   curve    -- Tao_curve_struct, pointer, optional: Pointer to the appropriate curve.
!   region   -- Tao_plot_region_struct, pointer, optional: Region found.
!-

subroutine tao_find_plot_by_region (err, where, plot, graph, curve, region, print_flag)

implicit none

type (tao_plot_region_struct), pointer, optional :: region
type (tao_plot_struct), pointer, optional :: plot
type (tao_graph_struct), pointer, optional :: graph
type (tao_curve_struct), pointer, optional :: curve
type (tao_plot_region_struct), pointer :: this_region

integer i, ix

character(*) where
character(32) plot_name, graph_name
character(28) :: r_name = 'tao_find_plot_by_region'

logical, optional :: print_flag
logical err

! Parse where argument

err = .false.

ix = index(where, ':')
if (ix == 0) then
  plot_name = where
  graph_name = ' '
else
  plot_name = where(1:ix-1)
  graph_name = where(ix+1:)
endif

! Match plot name to region

do i = 1, size(s%plot_page%region)

  this_region => s%plot_page%region(i)
  if (plot_name == this_region%name) exit

  if (i == size(s%plot_page%region)) then
    if (logic_option(.true., print_flag)) call out_io (s_error$, r_name, &
                                             'PLOT LOCATION NOT FOUND: ' // plot_name)
    err = .true.
    return
  endif

enddo

if (present(region)) region => this_region
if (present(plot)) plot => this_region%plot

! Find the graph and/or curve

call tao_find_graph_or_curve_from_plot (err, graph_name, this_region%plot, &
                                                            graph, curve, print_flag)

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_find_template_plot (err, where, plot, graph, curve, print_flag)
!
! Routine to find a template plot using the template name.
! A plot name is something like: where = "top"
!
! Input:
!   where    -- Character(*): Name to match to.
!   print_flag -- Logical, optional: If present and False then surpress error
!                   messages. Default is True.
!
! Output:
!   plot     -- Tao_plot_struct, pointer, optional: Pointer to the appropriate plot.
!   graph    -- Tao_graph_struct, pointer, optional: Pointer to the appropriate graph.
!   curve    -- Tao_curve_struct, pointer, optional: Pointer to the appropriate curve.
!-

subroutine tao_find_template_plot (err, where, plot, graph, curve, print_flag)

implicit none

type (tao_plot_struct), pointer, optional :: plot
type (tao_graph_struct), pointer, optional :: graph
type (tao_curve_struct), pointer, optional :: curve

integer i, j, ix

character(*) where
character(32) plot_name, graph_name
character(28) :: r_name = 'tao_find_template_plot'

logical, optional :: print_flag
logical err

! Find plot

err = .false.

ix = index(where, ':')
if (ix == 0) then
  plot_name = where
  graph_name = ' '
else
  plot_name = where(1:ix-1)
  graph_name = where(ix+1:)
endif

do i = 1, size(s%template_plot)
  if (plot_name == s%template_plot(i)%name) then
    if (present(plot)) plot => s%template_plot(i)
    exit
  endif
enddo

if (i == size(s%template_plot) + 1) then
  if (logic_option(.true., print_flag)) call out_io (s_error$, r_name, &
                                                  'PLOT LOCATION NOT FOUND: ' // where)
  err = .true.
  return
endif

call tao_find_graph_or_curve_from_plot (err, graph_name, s%template_plot(i), &
                                                               graph, curve, print_flag)

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_find_graph_or_curve_from_plot (err, where, plot, graph, curve, print_flag)
!
! Subroutine to find a graph or curve of a given plot.
!
! Input:
!   where    -- Character(*): Name to match to.
!   plot     -- Tao_plot_struct: Plot containing graph and/or curve.
!   print_flag -- Logical, optional: If present and False then surpress error
!                   messages. Default is True.
!
! Output:
!   err      -- logical: Set True on error. False otherwise.
!   graph    -- Tao_graph_struct, pointer, optional: Pointer to the appropriate graph.
!   curve    -- Tao_curve_struct, pointer, optional: Pointer to the appropriate curve.
!-

subroutine tao_find_graph_or_curve_from_plot (err, where, plot, graph, curve, print_flag)

type (tao_plot_struct) :: plot
type (tao_graph_struct), pointer, optional :: graph
type (tao_curve_struct), pointer, optional :: curve
type (tao_graph_struct), pointer :: this_graph

integer i, j, ix

character(*) where
character(32) plot_name, graph_name, curve_name
character(28) :: r_name = 'tao_find_graph_or_curve_from_plot'

logical, optional :: print_flag
logical err

!

if (present(graph)) nullify(graph)
if (present(curve)) nullify(curve)

ix = index(where, ':')
if (ix == 0) then
  graph_name = where
  curve_name = ' '
else
  graph_name = where(1:ix-1)
  curve_name = where(ix+1:)
endif

if (graph_name == ' ') return

do j = 1, size(plot%graph)
  this_graph => plot%graph(j)
  if (graph_name == this_graph%name) then
    if (present(graph)) graph => this_graph
    exit
  endif
enddo

if (j == size(plot%graph) + 1) then
  if (logic_option(.true., print_flag)) call out_io (s_error$, r_name, &
                                                  'GRAPH NOT FOUND: ' // where)
  err = .true.
endif

! Find the curve

if (curve_name == ' ') return
if (.not. present(curve)) return

do j = 1, size(this_graph%curve)
  if (curve_name == this_graph%curve(j)%name) then
    curve => this_graph%curve(j)
    return
  endif
enddo

if (logic_option(.true., print_flag)) call out_io (s_error$, r_name, &
                                                  'CURVE NOT FOUND: ' // where)
err = .true.

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_data_useit_plot_calc (graph, data)
!
! Subroutine to set the data for plotting.
!
! Input:
!
! Output:
!   data     -- Tao_data_struct:
!     %useit_plot -- True if good for plotting.
!                  = %exists & %good_plot (w/o measured & reference data)
!                  = %exists & %good_plot & %good_user & %good_meas (w/ meas data)
!                  = %exists & %good_plot & %good_user & %good_ref (w/ reference data)
!                  = %exists & %good_plot & %good_user & %good_meas & %good_ref 
!                                                        (w/ measured & reference data)
!-

subroutine tao_data_useit_plot_calc (graph, data)

implicit none

type (tao_graph_struct) graph
type (tao_data_struct) data(:)

!

data%useit_plot = data%exists .and. data%good_plot
if (any(graph%who%name == 'meas')) &
         data%useit_plot = data%useit_plot .and. data%good_user .and. data%good_meas
if (any(graph%who%name == 'ref'))  &
         data%useit_plot = data%useit_plot .and. data%good_user .and. data%good_ref

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_var_useit_plot_calc (graph, var)
!
! Subroutine to set the variables for plotting.
!
! Input:
!
! Output:
!   var     -- Tao_var_struct:
!     %useit_plot -- True if good for plotting.
!-

subroutine tao_var_useit_plot_calc (graph, var)

implicit none

type (tao_graph_struct) graph
type (tao_var_struct) var(:)

!

var%useit_plot = var%exists .and. var%good_user .and. var%good_plot &
                                                .and. var%good_var
if (any(graph%who%name == 'meas')) var%useit_plot = var%useit_plot
if (any(graph%who%name == 'ref'))  var%useit_plot = var%useit_plot

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_find_data (err, u, data_type, d2_ptr, d1_ptr, &
!                                        data_number, d_ptr, print_err)
!
! Routine to set data pointers to the correct data.
! Note: if, say, data_type = 'orbit' then d1_ptr will be nullified unless there
! if only one d1_data in which case d1_ptr will point to this.
!
! Input:
!   u            -- Tao_universe_struct
!   data_type    -- Character(*): the data name type. Eg: "orbit:x"
!   data_number  -- Character(*), optional: the data point index.
!                     If data_number = 'null' then d_ptr will be nullified.
!   print_err    -- Logical, optional: Print error message if data is 
!                     not found? Default is True.
!
! Output:
!   err     -- logical: err condition
!   d2_ptr  -- tao_d2_data_struct, optional, pointer: to the d2 data
!   d1_ptr  -- tao_d1_data_struct, optional: pointer to the d1 data
!   d_ptr   -- tao_data_struct, optional: pointer to the data point
!-

subroutine tao_find_data (err, u, data_type, d2_ptr, d1_ptr, &
                                    data_number, d_ptr, print_err)

implicit none

type (tao_universe_struct), target           :: u
type (tao_d2_data_struct), pointer, optional :: d2_ptr
type (tao_d1_data_struct), pointer, optional :: d1_ptr
type (tao_data_struct), pointer, optional    :: d_ptr
type (tao_d2_data_struct), pointer :: d2_pointer
type (tao_d1_data_struct), pointer :: d1_pointer

character(*)                                :: data_type
character(*), optional                      :: data_number
character(20) :: r_name = 'tao_find_data'
character(16) name, d2_name, d1_name

integer :: data_num, ios
integer i, ix, ix_plane

logical err
logical, optional :: print_err

! init

err = .false.

if (present(d2_ptr)) nullify(d2_ptr)
if (present(d1_ptr)) nullify(d1_ptr)
if (present(d_ptr)) nullify(d_ptr)

ix = index(data_type, ':')
if (ix == 0) then
  name = data_type
  d1_name = ' '
else
  name = data_type(1:ix-1)
  d1_name = data_type(ix+1:)
endif

! Point to the correct d2 data type 

call string_trim (name, name, ix) ! Strip off all whitespace

do i = 1, size(u%d2_data)
  if (name == u%d2_data(i)%name) then
    d2_pointer => u%d2_data(i) 
    if (present(d2_ptr)) d2_ptr => d2_pointer
    exit
  endif
  if (i == size(u%d2_data)) then
    if (logic_option (.true., print_err)) &
          call out_io (s_error$, r_name, "Couldn't find d2_data name: " // name)
    err = .true.
    return
  endif
enddo

! strip off all whitespace

call string_trim (d1_name, name, ix)
  
! point to the correct d1 data type

do i = 1, size(d2_pointer%d1)
  if (d1_name == d2_pointer%d1(i)%name) then
    d1_pointer => d2_pointer%d1(i)
    if (present(d1_ptr)) d1_ptr => d1_pointer
    exit
  endif
  if (i .eq. size(d2_pointer%d1)) then
    if (d1_name == ' ') then
      if (size(d2_pointer%d1) .eq. 1) then
        d1_pointer => d2_pointer%d1(1)
        if (present(d1_ptr)) d1_ptr => d1_pointer
      endif
      return
    endif
    if (logic_option (.true., print_err)) call out_io (s_error$, r_name, &
                                    "Couldn't find d1_data name: " // d1_name)
    err = .true.
    return  
  endif
enddo

! point to data point

if (.not. present(data_number)) return
if (data_number == 'null') return

read (data_number, '(i)', iostat = ios) data_num
if (ios /= 0) then
  call out_io (s_error$, r_name, "BAD DATA_NUMBER: " // data_number)
  err = .true.
  return  
endif
if (data_num < lbound(d1_ptr%d, 1) .or. data_num > ubound(d1_ptr%d, 1)) then
  call out_io (s_error$, r_name, "DATA NUMBER OUT OF RANGE: " // data_number)
  err = .true.
  return  
endif

if (present(d_ptr)) d_ptr => d1_pointer%d(data_num)

end subroutine tao_find_data

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine: tao_find_var (err, var_name, v1_ptr, var_number, v_ptr)
!
! find a v1 variable type, and variable data then point to it
!
!Input:
! u    -- tao_universe_struct
! var_name     -- character(*): the v1_var name
! var_number    -- integer: (optional) the variable data point.
!                     If var_number = 'null' then d_ptr will be nullified.
!
!Output:
! err    -- logical: err condition
! v1_ptr:  -- tao_v1_var_struct: pointer to the v1 variable
! v_ptr:        -- tao_var_struct: (optional) pointer to the variable data point
!-

subroutine tao_find_var (err, var_name, v1_ptr, var_number, v_ptr)

implicit none

logical err
type (tao_v1_var_struct), pointer, optional  :: v1_ptr
type (tao_v1_var_struct), pointer :: v1
type (tao_var_struct), pointer, optional     :: v_ptr

character(*)                                 :: var_name
character(*), optional                       :: var_number
character(20) :: r_name = 'tao_find_var'
character(32) v_name

integer i, ix, n_var, ios

err = .false.

! Strip off all whitespace

if (present(v1_ptr)) nullify (v1_ptr)
if (present(v_ptr)) nullify (v_ptr)

call string_trim(var_name, v_name, ix)
if (ix == 0) then
  err = .true.
  call out_io (s_error$, r_name, 'VARIABLE NAME NAME IS BLANK')
  return
endif

! Point to the correct v1 data type 

  do i = 1, size(s%v1_var)
    if (v_name == s%v1_var(i)%name) then
      v1 => s%v1_var(i) 
      if (present(v1_ptr)) v1_ptr => s%v1_var(i) 
      exit
    endif
    if (i .eq. size(s%v1_var)) then
      call out_io (s_error$, r_name, &
                        "COULD NOT FIND VARIABLE NAME: " // v_name)
      err = .true.
      return
    endif
  enddo

! point to variable data point

  if (.not. present(var_number)) return
  if (var_number == 'null') return

  read (var_number, '(i)', iostat = ios) n_var
  if (ios /= 0) then
    call out_io (s_error$, r_name, "BAD VAR_NUMBER: " // var_number)
    err = .true.
    return  
  endif
  if (n_var < lbound(v1%v, 1) .or. n_var > ubound(v1%v, 1)) then
    call out_io (s_error$, r_name, &
                                "VAR_NUMBER OUT OF RANGE: " // var_number)
    err = .true.
    return  
  endif

  if (present(v_ptr)) v_ptr => v1%v(n_var)

end subroutine tao_find_var

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_var_target_calc ()
! 
! Subroutine to calculate the variable target values (the values that they need
! to be set to to do a correction of the orbit, phase, etc.
!
! Input:
!
! Output:
!-

subroutine tao_var_target_calc ()

implicit none

type (tao_var_struct), pointer :: var

integer i, j

!

do j = 1, size(s%var)
  var => s%var(j)
  var%correction_value = var%meas_value + (var%design_value - var%model_value)
enddo

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

subroutine tao_set_var_model_value (var, value)

implicit none

type (tao_var_struct) var

real(rp) value
integer i

!

if (.not. (var%exists .and. var%good_var)) return

var%model_value = value
do i = 1, size(var%this)
  var%this(i)%model_ptr = value
enddo

s%global%lattice_recalc = .true.

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! subroutine tao_count_strings (string, pattern, num)
! 
! Subroutine to count the number of a specific pattern in the string
!
! Input:    
!  string    -- character(*): the string to look at
!  pattern   -- character(*): the search pattern
!
! Output:
!  num       -- integer: number of occurances
!-

subroutine tao_count_strings (string, pattern, num)

implicit none

character(*) string, pattern

integer num, len_string, len_pattern, i

num = 0
len_pattern = len(pattern)
len_string  = len(string)

do i = 1, len(string)
  if (i+len_pattern-1 .gt. len_string) return
  if (string(i:i+len_pattern-1) .eq. pattern) num = num + 1
enddo

end subroutine tao_count_strings

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_lat_bookkeeper (u, lattice)
!
! This will make sure all bookkeeping is up to date
!
! Input:
!  u                -- tao_universe_struct
!  lat              -- ring_string: lattice to bookkeep
!
! Output:
!  lat              -- ring_struct
!-

subroutine tao_lat_bookkeeper (u, lat)

implicit none

type (tao_universe_struct) :: u
type (ring_struct) :: lat

character(20) :: r_name = "tao_lat_bookkeeper"

  call lattice_bookkeeper (lat)

end subroutine tao_lat_bookkeeper

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine tao_to_real (word, value, err_flag)
!
! mathematically evaluates a character expression.
!
! Input:
!  expression   -- character(*): arithmetic expression
!  
! Output:
!  value        -- real(rp): value of arithmetic expression
!  err_flag     -- Logical: TRUE on error
!-

subroutine tao_to_real (expression, value, err_flag)

use random_mod

implicit none

type (eval_stack_struct) stk(200)

integer i_lev, i_op, i, ios

integer op_(200), ix_word, i_delim, i2, ix_word2

real(rp) value

character(*), intent(in) :: expression
character(100) phrase
character(1) delim
character(40) word, word2
character(16) :: r_name = "parse_expression"

logical delim_found, split, ran_function_pending
logical err_flag

! Don't destroy the input expression
  if (len(expression) .gt. len(phrase)) then
    call out_io (s_warn$, r_name, &
      "Expression cannot be longer than /I3/ characters", len(phrase))
    err_flag = .true.
    return
  endif
  read (expression, '(a)') phrase

! if phrase is blank then return 0.0
  call string_trim (phrase, phrase, ios)
  if (ios .eq. 0) then
    value = 0.0
    return
  endif
 
! The general idea is to rewrite the phrase on a stack in reverse polish.
! Reverse polish means that the operand goes last so that 2 * 3 is writen 
! on the stack as: [2, 3, *]

! The stack is called: stk
! Since operations move towards the end of the stack we need a separate
! stack called op_ which keeps track of what operations have not yet
! been put on stk.

! init

  err_flag = .false.
  i_lev = 0
  i_op = 0
  ran_function_pending = .false.

! parsing loop to build up the stack.

  parsing_loop: do

! get a word

    call word_read (phrase, '+-*/()^,:} ', word, ix_word, delim, &
                    delim_found, phrase)
    call str_upcase (word, word)
!   call get_next_word (word, ix_word, '+-*/()^,:} ', delim, delim_found)

    if (delim == '*' .and. word(1:1) == '*') then
      call out_io (s_warn$, r_name, 'EXPONENTIATION SYMBOL IS "^" AS OPPOSED TO "**" for!')
      err_flag = .true.
      return
    endif

    if (ran_function_pending .and. (ix_word /= 0 .or. delim /= ')')) then
          call out_io (s_warn$, r_name, &
             'RAN AND RAN_GAUSS DO NOT TAKE AN ARGUMENT')
      err_flag = .true.
      return
    endif

!--------------------------
! Preliminary: If we have split up something that should have not been split
! then put it back together again...

! just make sure we are not chopping a number in two, e.g. "3.5e-7" should not
! get split at the "-" even though "-" is a delimiter

    split = .true.         ! assume initially that we have a split number
    if (ix_word == 0) then
      split = .false.
    elseif (word(ix_word:ix_word) /= 'E') then
      split = .false.
    endif
    if (delim(1:1) /= '-' .and. delim(1:1) /= '+') split = .false.
    do i = 1, ix_word-1
      if (index('.0123456789', word(i:i)) == 0) split = .false.
    enddo

! If still SPLIT = .TRUE. then we need to unsplit

    if (split) then
      word = word(:ix_word) // delim
      call word_read (phrase, '+-*/()^,:}', word2, ix_word2, delim, &
                    delim_found, phrase)
      call str_upcase (word2, word2)
      word = word(:ix_word+1) // word2
      ix_word = ix_word + ix_word2
    endif

! Something like "lcav[lr(2).freq]" will get split on the "("

    if (delim == '(' .and. index(word, '[LR') /= 0) then
      call word_read (phrase, '+-*/(^,:}', word2, ix_word2, delim, &
                    delim_found, phrase)
      call str_upcase (word2, word2)
      word = word(:ix_word) // '(' // word2
      ix_word = ix_word + ix_word2
    endif

!---------------------------
! Now see what we got...

! For a "(" delim we must have a function

    if (delim == '(') then

      ran_function_pending = .false.
      if (ix_word /= 0) then
        select case (word)
        case ('SIN') 
          call pushit (op_, i_op, sin$)
        case ('COS') 
          call pushit (op_, i_op, cos$)
        case ('TAN') 
          call pushit (op_, i_op, tan$)
        case ('ASIN') 
          call pushit (op_, i_op, asin$)
        case ('ACOS') 
          call pushit (op_, i_op, acos$)
        case ('ATAN') 
          call pushit (op_, i_op, atan$)
        case ('ABS') 
          call pushit (op_, i_op, abs$)
        case ('SQRT') 
          call pushit (op_, i_op, sqrt$)
        case ('LOG') 
          call pushit (op_, i_op, log$)
        case ('EXP') 
          call pushit (op_, i_op, exp$)
        case ('RAN') 
          call pushit (op_, i_op, ran$)
          ran_function_pending = .true.
        case ('RAN_GAUSS') 
          call pushit (op_, i_op, ran_gauss$)
          ran_function_pending = .true.
        case default
          call out_io (s_warn$, r_name, &
               'UNEXPECTED CHARACTERS ON RHS BEFORE "(": ')
          err_flag = .true.
          return
        end select
      endif

      call pushit (op_, i_op, l_parens$)
      cycle parsing_loop

! for a unary "-"

    elseif (delim == '-' .and. ix_word == 0) then
      call pushit (op_, i_op, unary_minus$)
      cycle parsing_loop

! for a unary "+"

      call pushit (op_, i_op, unary_plus$)
      cycle parsing_loop

! for a ")" delim

    elseif (delim == ')') then
      if (ix_word == 0) then
        if (.not. ran_function_pending) call out_io (s_warn$, r_name, &
              'CONSTANT OR VARIABLE MISSING BEFORE ")"')
        err_flag = .true.
        return
      else
        read (word, *, iostat = ios) value
        if (ios .ne. 0) then
          call out_io (s_warn$, r_name, &
                "This doesn't seem to be a number: /a/", word)
          err_flag = .true.
          return
        endif
        call pushit (stk%type, i_lev, numeric$)
        stk(i_lev)%value = value
      endif

      do
        do i = i_op, 1, -1     ! release pending ops
          if (op_(i) == l_parens$) exit          ! break do loop
          call pushit (stk%type, i_lev, op_(i))
        enddo

        if (i == 0) then
          call out_io (s_warn$, r_name, 'UNMATCHED ")" ON RHS')
          err_flag = .true.
          return
        endif

        i_op = i - 1

        call word_read (phrase, '+-*/()^,:}', word, ix_word, delim, &
                    delim_found, phrase)
        call str_upcase (word, word)
        if (ix_word /= 0) then
          call out_io (s_warn$, r_name, &
                   'UNEXPECTED CHARACTERS ON RHS AFTER ")"')
          err_flag = .true.
          return
        endif

        if (delim /= ')') exit  ! if no more ')' then no need to release more
      enddo


      if (delim == '(') then
        call out_io (s_warn$, r_name,  &
                    '")(" CONSTRUCT DOES NOT MAKE SENSE')
        err_flag = .true.
        return
      endif

! For binary "+-/*^" delims

    else
      if (ix_word == 0) then
        call out_io (s_warn$, r_name, 'CONSTANT OR VARIABLE MISSING')
        err_flag = .true.
        return
      endif
      read (word, *, iostat = ios) value
      if (ios .ne. 0) then
        call out_io (s_warn$, r_name, &
              "This doesn't seem to be a number: /a/", word)
        err_flag = .true.
        return
      endif
      call pushit (stk%type, i_lev, numeric$)
      stk(i_lev)%value = value
    endif

! If we are here then we have an operation that is waiting to be identified

    if (.not. delim_found) delim = ':'

    select case (delim)
    case ('+')
      i_delim = plus$
    case ('-')
      i_delim = minus$
    case ('*')
      i_delim = times$
    case ('/')
      i_delim = divide$
    case (')')
      i_delim = r_parens$
    case ('^')
      i_delim = power$
    case (',', '}', ':')
      i_delim = no_delim$
    case default
        call out_io (s_error$, r_name, 'INTERNAL ERROR')
        call err_exit
    end select

! now see if there are operations on the OP_ stack that need to be transferred
! to the STK_ stack

    do i = i_op, 1, -1
      if (eval_level(op_(i)) >= eval_level(i_delim)) then
        if (op_(i) == l_parens$) then
          call out_io (s_warn$, r_name, 'UNMATCHED "("')
          err_flag = .true.
          return
        endif
        call pushit (stk%type, i_lev, op_(i))
      else
        exit
      endif
    enddo

! put the pending operation on the OP_ stack

    i_op = i
    if (i_delim == no_delim$) then
      exit parsing_loop
    else
      call pushit (op_, i_op, i_delim)
    endif

  enddo parsing_loop

!------------------------------------------------------------------
! now go through the stack and perform the operations

  if (i_op /= 0) then
    call out_io (s_warn$, r_name, 'UNMATCHED "("')
    err_flag = .true.
    return
  endif

  if (i_lev == 0) then
    call out_io (s_warn$, r_name, 'NO VALUE FOUND')
    err_flag = .true.
    return
  endif

  i2 = 0
  do i = 1, i_lev
    if (stk(i)%type == numeric$) then
      i2 = i2 + 1
      stk(i2)%value = stk(i)%value
    elseif (stk(i)%type == unary_minus$) then
      stk(i2)%value = -stk(i2)%value
    elseif (stk(i)%type == unary_plus$) then
      stk(i2)%value = stk(i2)%value
    elseif (stk(i)%type == plus$) then
      stk(i2-1)%value = stk(i2-1)%value + stk(i2)%value
      i2 = i2 - 1
    elseif (stk(i)%type == minus$) then
      stk(i2-1)%value = stk(i2-1)%value - stk(i2)%value
      i2 = i2 - 1
    elseif (stk(i)%type == times$) then
      stk(i2-1)%value = stk(i2-1)%value * stk(i2)%value
      i2 = i2 - 1
    elseif (stk(i)%type == divide$) then
      if (stk(i2)%value == 0) then
        call out_io  (s_warn$, r_name, 'DIVIDE BY 0 ON RHS')
        err_flag = .true.
        return
      endif
      stk(i2-1)%value= stk(i2-1)%value / stk(i2)%value
      i2 = i2 - 1
    elseif (stk(i)%type == power$) then
      stk(i2-1)%value = stk(i2-1)%value**stk(i2)%value
      i2 = i2 - 1
    elseif (stk(i)%type == sin$) then
      stk(i2)%value = sin(stk(i2)%value)
    elseif (stk(i)%type == cos$) then
      stk(i2)%value = cos(stk(i2)%value)
    elseif (stk(i)%type == tan$) then
      stk(i2)%value = tan(stk(i2)%value)
    elseif (stk(i)%type == asin$) then
      stk(i2)%value = asin(stk(i2)%value)
    elseif (stk(i)%type == acos$) then
      stk(i2)%value = acos(stk(i2)%value)
    elseif (stk(i)%type == atan$) then
      stk(i2)%value = atan(stk(i2)%value)
    elseif (stk(i)%type == abs$) then
      stk(i2)%value = abs(stk(i2)%value)
    elseif (stk(i)%type == sqrt$) then
      stk(i2)%value = sqrt(stk(i2)%value)
    elseif (stk(i)%type == log$) then
      stk(i2)%value = log(stk(i2)%value)
    elseif (stk(i)%type == exp$) then
      stk(i2)%value = exp(stk(i2)%value)
    elseif (stk(i)%type == ran$) then
      i2 = i2 + 1
      call ran_uniform(stk(i2)%value)
    elseif (stk(i)%type == ran_gauss$) then
      i2 = i2 + 1
      call ran_gauss(stk(i2)%value)
    else
      call out_io (s_warn$, r_name, 'INTERNAL ERROR')
      err_flag = .true.
      return
    endif
  enddo


  if (i2 /= 1) call out_io (s_warn$, r_name, 'INTERNAL ERROR')

  value = stk(1)%value


contains

!-------------------------------------------------------------------------

subroutine pushit (stack, i_lev, value)

  implicit none

  integer stack(:), i_lev, value

  character(6) :: r_name = "pushit"

!

  i_lev = i_lev + 1

  if (i_lev > size(stack)) then
    call out_io (s_warn$, r_name, 'STACK OVERFLOW.')
    call err_exit
  endif

  stack(i_lev) = value

end subroutine pushit
                       
end subroutine tao_to_real

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_to_int (str, i_int, err)
! 
! Converts a string to an integer
!
! If the string str is blank then i_int = 0
!-

subroutine tao_to_int (str, i_int, err)

character(*) str
integer ios, i_int
logical err
character(12) :: r_name = "tao_to_int"

!

  call string_trim (str, str, ios)
  if (ios .eq. 0) then
    i_int = 0
    return
  endif
 
  err = .false.
  read (str, *, iostat = ios) i_int

  if (ios /= 0) then
    call out_io (s_error$, r_name, 'EXPECTING INTEGER: ' // str)
    err = .true.
    return
  endif

end subroutine

end module tao_utils
