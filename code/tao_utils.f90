!+
! Module tao_utils
!
! helper subroutines available for communal usage.
!-

module tao_utils

use tao_struct
use tao_interface
use bmad

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
! Subroutine tao_point_v1_to_var (ip, ii, n, n_var)
!
! used for arbitrary variable pointer indexing
!
! ip       -- tao_var_struct: the pointer
! ii:      -- tao_var_struct: the variable
! n:       -- integer: starting index for the pointer
!-

subroutine tao_point_v1_to_var (ip, ii, n, n_var)

implicit none

integer n, i, n_var

type (tao_var_struct), pointer :: ip(:)
type (tao_var_struct), target :: ii(n:)

ip => ii

forall (i = lbound(ii, 1):ubound(ii, 1)) 
  ii(i)%ix_v1 = i
  ii(i)%ix_var = n_var + i - n
end forall

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
! Subroutine tao_locate_element (string, lattice, ix_ele) 
!
! Subroutine to find a lattice element.
!
! Input:
!   string  -- Character(*): String with element name
!   lattice -- Ring_struct: Lattice to search
!
! Output:
!   ix_ele  -- Integer: Index of element. Set to -1 if element not found.
!-

subroutine tao_locate_element (string, lattice, ix_ele)

implicit none

type (ring_struct) lattice

integer ix_ele, ios, ix

character(*) string
character(16) ele_name
character(20) :: r_name = 'tao_find_element'

! If it is a number translate it:

call str_upcase (ele_name, string)
call string_trim (ele_name, ele_name, ix)

if (ix == 0) then
  ix_ele = -1
  call out_io (s_error$, r_name, 'ELEMENT NAME IS BLANK')
  return
endif

read (ele_name, *, iostat = ios) ix_ele
if (ios .eq. 0) then
  !it's a number
  if (ix_ele < 0 .or. ix_ele > lattice%n_ele_max) then
    ix_ele = -1
    call out_io (s_error$, r_name, 'ELEMENT INDEX OUT OF RANGE: ' // ele_name)
  endif
  return
endif

call element_locator (ele_name, lattice, ix_ele)
if (ix_ele < 0) call out_io (s_error$, r_name, 'ELEMENT NOT FOUND: ' // string)

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_pointer_to_var_in_lattice (var, this, ix_uni, err)
! 
! Routine to set a pointer to the appropriate variable in a lattice
!
! Input:
!   var   -- Tao_var_struct: Structure has the info of where to point.
!
! Output:
!   err   -- Logical: Set True if there is an error. False otherwise.
!-

subroutine tao_pointer_to_var_in_lattice (var, this, ix_uni, err)

implicit none

type (tao_var_struct) var
type (tao_universe_struct), pointer :: u
type (tao_this_var_struct) this

integer ix, ie, ix_uni
logical, optional :: err
logical error
character(30) :: r_name = 'tao_pointer_to_var_in_lattice'

! locate element

if (present(err)) err = .true.

u => s%u(ix_uni)
call element_locator (var%ele_name, u%model, ie)
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

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_find_plot (err, plots, by_who, where, plot, graph)
!
! Routine to find a plot or graph by name.
! A plot name is something like: where = "top"
! A graph name is something like: where = "top:x"
!
! Input:
!   plots(:) -- Tao_plot_struct: Array of plots to look at.
!     %name        -- Match name if by_who = "BY_TYPE"
!     %region%name -- Match name if by_who = "BY_REGION"
!   by_who   -- Character(*): Either: "BY_REGION" or "BY_TYPE".
!   where    -- Character(*): Name to match to.
!
! Output:
!   plot     -- Tao_plot_struct, pointer, optional: Pointer to the appropriate plot.
!   graph    -- Tao_graph_struct, pointer, optional: Pointer to the appropriate graph.
!-

subroutine tao_find_plot (err, plots, by_who, where, plot, graph)

implicit none

type (tao_plot_struct), target :: plots(:)
type (tao_plot_struct), pointer, optional :: plot
type (tao_graph_struct), pointer, optional :: graph

integer i, j, ix

character(*) by_who, where
character(16) plot_name, graph_name
character(20) :: r_name = 'tao_find_plot'

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

do i = 1, size(plots)

  if (by_who == 'BY_REGION') then
    if (plot_name == plots(i)%region%name) exit
  elseif (by_who == 'BY_TYPE') then
    if (plot_name == plots(i)%name) exit
  else
    call out_io (s_error$, r_name, 'BAD "BY_WHO" ARGUMENT: ' // by_who)
  endif

  if (i == size(plots)) then
    call out_io (s_error$, r_name, 'PLOT LOCATION NOT FOUND: ' // plot_name)
    err = .true.
    return
  endif

enddo

if (present(plot)) plot => plots(i)

! Find graph

if (.not. present(graph)) return 
if (graph_name == ' ')  then
  nullify(graph)
  return
endif

do j = 1, size(plots(i)%graph)
  graph => plots(i)%graph(j)
  if (graph_name == graph%name) return
enddo

call out_io (s_error$, r_name, 'GRAPH NOT FOUND: ' // where)
err = .true.
nullify(graph)

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_useit_plot_calc (plot, data)
!
! Subroutine to set the data for plotting.
!
! Input:
!
! Output:
!   data     -- Tao_data_struct:
!     %useit_plot -- True if good for plotting.
!-

subroutine tao_useit_plot_calc (plot, data)

implicit none

type (tao_plot_struct) plot
type (tao_data_struct) data(:)

!

data%useit_plot = data%exists .and. data%good_user .and. data%good_plot
if (any(plot%who%name == 'meas')) data%useit_plot = data%useit_plot .and. data%good_data
if (any(plot%who%name == 'ref'))  data%useit_plot = data%useit_plot .and. data%good_ref

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_find_data (err, u, data_type, d2_ptr, d1_ptr, &
!                                        data_number, d_ptr, print_err)
!
! Routine to set data pointers to the correct data.
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
character(20) :: r_name = 'TAO_FIND_DATA'
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

! If d1_name == ' ' nothing more to be done.
! Except if there is only 1 d2_pointer%d1 then point to that.

if (d1_name == ' ') then
  if (size(d2_pointer%d1) == 1 .and. present(d1_ptr)) then
    d1_ptr => d2_pointer%d1(1) 
  else
    return
  endif
endif

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

var%model_value = value
do i = 1, size(var%this)
  var%this(i)%model_ptr = value
enddo

end subroutine

end module tao_utils
