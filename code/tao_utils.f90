!+
! Module tao_utils
!
! helper subroutines available for communal usage.
!-

module tao_utils

use tao_struct
use tao_interface
use bmad
use output_mod
use lat_ele_loc_mod

integer, parameter :: apparent_emit$ = 1, projected_emit$ = 2

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

subroutine tao_key_info_to_str (ix_key, ix_min_key, ix_max_key, key_str, header_str)

implicit none

type (tao_var_struct), pointer :: var

integer ix_key, ix_min_key, ix_max_key
integer i, n, m, p, ix_var
integer :: j_var1, j_att

real(rp) :: y_here, norm, v, x1, x2, y1, y2

character(*) key_str
character(*) header_str
character(60) fmt, fmt2
character(11) model_str, val0_str, delta_str
character(4) exp_str
character(24) :: r_name = 'tao_key_info_to_str'

! Compute widths of var1 and attrib fields.

j_var1 = 4
j_att = 5

do i = ix_min_key, ix_max_key
  if (i > ubound(s%key, 1)) cycle
  ix_var = s%key(i)
  if (ix_var == 0) cycle
  j_var1 = max(j_var1, len_trim(tao_var1_name(s%var(ix_var))))
  j_att  = max(j_att,  len_trim(tao_var_attrib_name(s%var(ix_var))))
enddo

! Write header 

write (fmt, '(a, i5, a, i2, a)') '(a, ', j_var1-2, 'x, a, ', j_att, 'x, a)'
write (header_str, fmt) 'Name', 'Attrib', 'Value     Value0      Delta Opt'

! Write key info

key_str = ''
if (ix_key > ubound(s%key, 1)) return
  
ix_var = s%key(ix_key)
if (ix_var == 0) return

var => s%var(ix_var)
v = max(abs(var%model_value), abs(var%key_val0), abs(var%key_delta))
if (v == 0) then
  n = 0
  m = 2
else
  m = 1.001 * log10(v)
  n = 3 * floor(m/3.0)
  p = 3 - (m - n)
endif

if (m >= -1 .and. m <= 1) then
  fmt2 = '(f11.4, a0)'
  n = 0
elseif (m == 2) then
  fmt2 = '(f11.2, a0)'
  n = 0
else
  write (fmt2, '(a, i1, a)') '(f7.', p, ', a)'
  write (exp_str, '(a, i3.2)') 'E', n
  if (exp_str(2:2) == ' ') exp_str(2:2) = '+'
endif

write (model_str, fmt2) var%model_value / 10.0**n, exp_str
write (val0_str,  fmt2) var%key_val0 / 10.0**n, exp_str
write (delta_str, fmt2) var%key_delta / 10.0**n, exp_str

write (fmt, '(3(a, i2.2))') '(a', j_var1, ', 2x, a', j_att, ', 3a11, 3x, l1)'
write (key_str, fmt) tao_var1_name(var), tao_var_attrib_name(var), model_str, &
                        val0_str, delta_str, var%useit_opt

end subroutine tao_key_info_to_str

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_data_check (err)
!
! Routine to do some checking of data.
!-

subroutine tao_data_check (err)

implicit none

type (tao_data_struct), pointer :: datum
integer iu, id
logical err
character(16) :: r_name = 'tao_data_check'

!

err = .false.

do iu = lbound(s%u, 1), ubound(s%u, 1)
  do id = 1, size(s%u(iu)%data)
    datum => s%u(iu)%data(id)
    if (datum%merit_type(1:4) == 'int_' .and. &
            (s%global%opt_with_ref .or. s%global%opt_with_base)) then
      call out_io (s_error$, r_name, &
                        'BAD DATUM INTEGRATION FOR: ' // tao_datum_name(datum))
      err = .true.
    endif
  enddo
enddo

end subroutine tao_data_check

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_pick_universe (name_in, name_out, picked, err, ix_uni)
!
! Subroutine to pick what universe the data name is comming from.
! Examples:
!   "*@..."           -- Choose all universes.
!   "3@..."           -- Choose universe 3. 
!   "[1:30,34]@..."   -- Choose universes 1 through 30 and 34
!   No "@" in name    -- Choose universe s%com%default_universe.
!
! Also see:
!   tao_pointer_to_universe
!
! Input:
!   name_in    -- Character(*): data name.
!
! Output:
!   name_out   -- Character(*): name_in without any "n@" beginning.
!   picked(:)  -- Logical, allocatable: Array showing picked universes.
!                   The array will be resized if necessary.
!   err        -- Logical: Set True if an error is detected.
!   ix_uni     -- Integer, optional: Set to the picked universe with the highest index.
!-

subroutine tao_pick_universe (name_in, name_out, picked, err, ix_uni)

implicit none

character(*) name_in, name_out
character(20) :: r_name = 'tao_pick_universe'
character(40) uni

integer, optional :: ix_uni
integer i, ix, n, ios, iu, num, ic

logical, allocatable :: picked(:)
logical, allocatable :: p(:)
logical err

! Init

call re_allocate2 (picked, lbound(s%u, 1), ubound(s%u, 1))
call re_allocate2 (p, -1, ubound(s%u, 1))

err = .false.
picked = .false.
p = .false.
if (present(ix_uni)) ix_uni = -1

! No "@" then simply choose s%com%default_universe.

ix = index (name_in, '@')
ic = index (name_in, '::')
if (ix == 0 .or. (ic /= 0 .and. ix > ic)) then
  picked (s%com%default_universe) = .true.
  name_out = name_in
  if (present(ix_uni)) ix_uni = s%com%default_universe
  return
endif

! Here whn "@" is found...

uni = name_in(:ix-1)
name_out = name_in(ix+1:)

! Strip off '[' and ']'

if (uni(1:1) == '[' .and. uni(ix-1:ix-1) == ']') uni = uni(2:ix-2)

if (uni == '*') then
  picked = .true.
  if (present(ix_uni)) ix_uni = lbound(s%u, 1)
  return
endif


call location_decode (uni, p, lbound(p, 1), num)
if (num == -1) then
  call out_io (s_error$, r_name, "BAD UNIVERSE NUMBER: " // uni)
  err = .true.
  return
endif

do i = lbound(p, 1), ubound(p, 1)
  if (.not. p(i)) cycle
  iu = tao_universe_number (i)
  if (iu < lbound(s%u, 1) .or. iu > ubound(s%u, 1)) then
    call out_io (s_error$, r_name, "NUMBER DOES NOT CORRESPOND TO A UNIVERSE: " // uni)
    err = .true.
    return
  endif
  picked(iu) = .true.
  if (present(ix_uni)) ix_uni = iu
enddo

end subroutine tao_pick_universe

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_locate_all_elements (ele_list, eles, err, ignore_blank) 
!
! Subroutine to find the lattice elements in the lattice corresponding to the ele_list argument. 
!
! See also tao_locate_elements for a routine that only searches one given universe.
!
! Input:
!   ele_list     -- Character(*): String with element names using element list format.
!   ignore_blank -- Logical, optional: If present and true then do nothing if
!     ele_list is blank. otherwise treated as an error.
!
! Output:
!   eles  -- ele_pointer_struct(:), allocatable: Array of elements in the model lat. 
!     %id  -- Set to universe number.
!   err   -- Logical: Set true on error.
!-

subroutine tao_locate_all_elements (ele_list, eles, err, ignore_blank)

implicit none

type (tao_universe_struct), pointer :: u
type (ele_pointer_struct), allocatable :: eles(:)
type (ele_pointer_struct), allocatable :: this_eles(:)

integer i, ix, ix_word, n_eles, n0

character(*) ele_list
character(100) ele_name, word
character(20) :: r_name = 'tao_locate_all_elements'
character(1) delim

logical err, delim_found
logical, optional :: ignore_blank
logical, allocatable :: picked(:)

!

err = .true.

call re_allocate_eles (eles, 0, exact = .true.)

call str_upcase (ele_name, ele_list)
call string_trim (ele_name, ele_name, ix)

if (ix == 0 .and. logic_option(.false., ignore_blank)) return

if (ix == 0) then
  call out_io (s_error$, r_name, 'ELEMENT NAME IS BLANK')
  return
endif

! Loop over all items in the element list

do 
  call word_read (ele_list, ', ', word, ix_word, delim, delim_found, ele_list)
  if (ix_word == 0) exit
  call tao_pick_universe (word, word, picked, err)
  if (err) return
  do i = lbound(s%u, 1), ubound(s%u, 1)
    if (.not. picked(i)) cycle
    call lat_ele_locator (word, s%u(i)%model%lat, this_eles, n_eles, err)
    if (err) return
    if (n_eles == 0) cycle
    n0 = size(eles)
    call re_allocate_eles (eles, n0+n_eles, .true., .true.)
    eles(n0+1:n0+n_eles) = this_eles(1:n_eles)
    eles(n0+1:n0+n_eles)%id = i
  enddo
enddo

if (size(eles) == 0) then
  call out_io (s_error$, r_name, 'ELEMENT(S) NOT FOUND: ' // ele_list)
  err = .true.
  return
endif

end subroutine tao_locate_all_elements

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_locate_elements (ele_list, ix_universe, eles, err, lat_type, ignore_blank, print_err, ix_dflt_branch) 
!
! Subroutine to find the lattice elements in the lattice
! corresponding to the ele_list argument. 
!
! Note: If ele_list can contain elements from different universes, use
! the routine:
!   tao_locate_all_elements 
!
! Input:
!   ele_list       -- Character(*): String with element names using element list format.
!   ix_universe    -- Integer: Universe to search. 0 => search s%com%default_universe.
!   lat_type       -- Integer, optional: model$ (default), design$, or base$.
!   ignore_blank   -- Logical, optional: If present and true then do nothing if
!                     ele_list is blank. otherwise treated as an error.
!   print_err      -- Logical, optional: If present and False then do not print error messages.
!   ix_dflt_branch -- Integer, optional: If present and positive then use this as the branch index 
!                       for elements specified using an integer index (EG: "43").
!                       If not present or -1 the default branch is branch 0.
!
! Output:
!   eles  -- ele_pointer_struct(:), allocatable: Array of elements in the model lat. 
!   err   -- Logical: Set true on error.
!-

subroutine tao_locate_elements (ele_list, ix_universe, eles, err, lat_type, ignore_blank, &
                                                                 print_err, above_ubound_is_err, ix_dflt_branch)

implicit none

type (tao_universe_struct), pointer :: u
type (tao_lattice_struct), pointer :: tao_lat
type (ele_pointer_struct), allocatable :: eles(:)

integer, optional :: lat_type, ix_dflt_branch
integer ios, ix, ix_universe, num, i, i_ix_ele, n_loc

character(*) ele_list
character(200) ele_name
character(20) :: r_name = 'tao_locate_elements'

logical err, printit
logical, optional :: ignore_blank, print_err, above_ubound_is_err

! 

err = .true.
printit = logic_option (.true., print_err)

call re_allocate_eles (eles, 0, exact = .true.)

call str_upcase (ele_name, ele_list)
call string_trim (ele_name, ele_name, ix)

if (ix == 0 .and. logic_option(.false., ignore_blank)) return

if (ix == 0) then
  if (printit) call out_io (s_error$, r_name, 'ELEMENT NAME IS BLANK')
  return
endif

u => tao_pointer_to_universe (ix_universe)
if (.not. associated(u)) return

tao_lat => tao_pointer_to_tao_lat (u, lat_type)

call lat_ele_locator (ele_name, tao_lat%lat, eles, n_loc, err, above_ubound_is_err, ix_dflt_branch)
if (err) return

if (n_loc == 0) then
  if (printit) call out_io (s_error$, r_name, 'ELEMENT(S) NOT FOUND: ' // ele_list)
  err = .true.
  return
endif

call re_allocate_eles (eles, n_loc, .true., .true.)

end subroutine tao_locate_elements

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Function tao_pointer_to_tao_lat (u, lat_type) result (tao_lat)
!
! Routine to set a pointer to a tao_lat.
!
! Also see:
!   tao_pointer_to_universe
!
! Input:
!   u         -- Tao_universe_struct: Universe to work with
!   lat_type  -- Integer, optional: model$ (default), design$, or base$.
!
! Output:
!   tao_lat   -- tao_lattice_struct, pointer: Tao_lat pointer. 
!                   Points to u%model, u%design, or u%base
!-

function tao_pointer_to_tao_lat (u, lat_type) result (tao_lat)

implicit none

type (tao_universe_struct), target :: u
type (tao_lattice_struct), pointer :: tao_lat
integer, optional :: lat_type
character(28) :: r_name = 'tao_pointer_to_tao_lat'

!

if (.not. present(lat_type)) then
  tao_lat => u%model
  return
endif

select case (lat_type)
case (design$)
  tao_lat => u%design
case (model$)
  tao_lat => u%model
case (base$)
  tao_lat => u%base
case default
  call err_exit
end select

end function tao_pointer_to_tao_lat

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_find_plot_region (err, where, region, print_flag)
!
! Routine to find a region using the region name.
!
! Input:
!   where      -- Character(*): Region name.
!   print_flag -- Logical, optional: If present and False then surpress error
!                   messages. Default is True.
!
! Output:
!   err      -- logical: Set True on error. False otherwise.
!   region   -- Tao_plot_region_struct, pointer: Region found.
!-

subroutine tao_find_plot_region (err, where, region, print_flag)

implicit none

type (tao_plot_region_struct), pointer :: region

integer i, ix

character(*) where
character(40) plot_name, graph_name
character(28) :: r_name = 'tao_find_plot_region'

logical, optional :: print_flag
logical err

! Parse where argument

ix = index(where, '.')
if (ix == 0) then
  plot_name = where
else
  plot_name = where(1:ix-1)
endif

! Match plot name to region

err = .false.

call match_word (plot_name, s%plot_page%region%name, ix, .true.)

if (ix < 1) then
  if (logic_option(.true., print_flag)) call out_io (s_error$, r_name, &
                                    'PLOT LOCATION NOT FOUND: ' // plot_name)
  err = .true.
else
  region => s%plot_page%region(ix)  
endif

end subroutine tao_find_plot_region

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_find_plots (err, name, where, plot, graph, curve, print_flag, always_allocate)
!
! Routine to find a plot using the region or plot name.
! A region or plot name is something like: name = "top"
! A graph name is something like: name  = "top.x"
! A curve name is something like: name  = "top.x.c1"
! The wild card "*" can be used so name = "top.*.c1" could 
!   return "top.x.c1" and "top.y.c1".
!
! Unless always_allocate = T, The graph(:) array will only be allocated if 
! the graph portion of name is not blank.
!   For example: name = "top" will leave graph(:) unallocated.
! Unless always_allocate = T, The curve(:) array will only be allocated if 
! the curve portion of name is not blank.
!   For example: name = "top.x" will leave curve(:) unallocated.
!
! Input:
!   name       -- Character(*): Name of plot or region.
!   where      -- Character(*): Where to look: 'TEMPLATE', 'REGION', or 'BOTH'
!                   For where = 'BOTH', if something is found in a plot region,
!                   then the templates will not be searched
!   print_flag -- Logical, optional: If present and False then surpress error
!                   messages. Default is True.
!   always_allocate 
!              -- Logical, optional: If present and True then always allocate 
!                   graph(:) and curve(:) arrays except if there is an error. 
!
! Output:
!   err      -- logical: Set True on error. False otherwise.
!   plot(:)  -- Tao_plot_array_struct, allocatable, optional: Array of plots.
!   graph(:) -- Tao_graph_array_struct, allocatable, optional: Array of graphs.
!   curve(:) -- Tao_curve_array_struct, allocatable, optional: Array of curves.
!-

subroutine tao_find_plots (err, name, where, plot, graph, curve, print_flag, always_allocate)

implicit none

type (tao_plot_array_struct), allocatable, optional :: plot(:)
type (tao_graph_array_struct), allocatable, optional :: graph(:)
type (tao_curve_array_struct), allocatable, optional :: curve(:)
type (tao_plot_array_struct), allocatable :: p(:)
type (tao_graph_array_struct), allocatable :: g(:)
type (tao_curve_array_struct), allocatable :: c(:)

integer i, j, k, ix, np, ng, nc

character(*) name, where
character(40) plot_name, graph_name, curve_name
character(28) :: r_name = 'tao_find_plots'

logical, optional :: print_flag, always_allocate
logical err

! Init

if (present(plot)) then
  if (allocated(plot)) deallocate(plot)
endif

if (present(graph)) then
  if (allocated(graph)) deallocate(graph)
endif

if (present(curve)) then
  if (allocated(curve)) deallocate(curve)
endif

if (allocated(p)) deallocate(p)
if (allocated(g)) deallocate(g)
if (allocated(c)) deallocate(c)

! Error check

if (where /= 'REGION' .and. where /= 'BOTH' .and. where /= 'TEMPLATE') then
  if (logic_option(.true., print_flag)) call out_io (s_fatal$, r_name, &
                                             'BAD "WHERE" LOCATION: ' // where)
  call err_exit
endif

! Parse name argument

err = .false.

if (name == "") then
  if (logic_option(.true., print_flag)) &
                call out_io (s_error$, r_name, 'BLANK "WHERE" LOCATION')
  err = .true.
  return
endif

ix = index(name, '.')
if (ix == 0) then
  plot_name = name
  graph_name = ' '
else
  plot_name = name(1:ix-1)
  graph_name = name(ix+1:)
endif

! Match name to region or plot and count how many there are.

np = 0

if (where == 'REGION' .or. where == 'BOTH') then
  do i = 1, size(s%plot_page%region)
    if (index(s%plot_page%region(i)%name, trim(plot_name)) == 1 .or. &
        index(s%plot_page%region(i)%plot%name, trim(plot_name)) == 1 .or. plot_name == '*') np = np + 1
  enddo
endif

if (where == 'TEMPLATE' .or. (where == 'BOTH' .and. np == 0)) then
  do i = 1, size(s%plot_page%template)
    if (index(s%plot_page%template(i)%name, trim(plot_name)) == 1 .or. plot_name == '*') np = np + 1
  enddo
endif

! Allocate space

if (np == 0) then
  if (logic_option(.true., print_flag)) call out_io (s_error$, r_name, &
                                             'PLOT NOT FOUND: ' // plot_name)
  err = .true.
  return
endif

allocate (p(np))
if (present(plot)) allocate(plot(np))

! Now that we have counted and allocated the matches set the plot pointers.

np = 0

if (where == 'REGION' .or. where == 'BOTH') then
  do i = 1, size(s%plot_page%region)
    if (index(s%plot_page%region(i)%name, trim(plot_name)) == 1 .or. &
        index(s%plot_page%region(i)%plot%name, trim(plot_name)) == 1 .or. plot_name == '*') then
      np = np + 1
      p(np)%p => s%plot_page%region(i)%plot
    endif
  enddo
endif

if (where == 'TEMPLATE' .or. (where == 'BOTH' .and. np == 0)) then
  do i = 1, size(s%plot_page%template)
    if (index(s%plot_page%template(i)%name, trim(plot_name)) == 1 .or. plot_name == '*') then
      np = np + 1
      p(np)%p => s%plot_page%template(i)
    endif
  enddo
endif

if (present(plot)) plot = p

! Find the number of graphs and allocate.

if (.not. present(graph) .and. .not. present(curve)) return

ix = index(graph_name, '.')
if (ix == 0) then
  curve_name = ' '
else
  curve_name = graph_name(ix+1:)
  graph_name = graph_name(1:ix-1)
endif

if (logic_option(.false., always_allocate) .and. graph_name == ' ') graph_name = '*'
if (graph_name == ' ') return

ng = 0
do i = 1, np
  do j = 1, size(p(i)%p%graph)
    if (p(i)%p%graph(j)%name == graph_name .or. graph_name == '*') ng = ng + 1
  enddo
enddo

if (ng == 0) then
  if (logic_option(.true., print_flag)) call out_io (s_error$, r_name, &
                  'GRAPH NOT FOUND: ' // trim(plot_name) // '.' // graph_name)
  err = .true.
  return
endif

allocate (g(ng))
if (present(graph)) allocate (graph(ng))

! Now that we have counted and allocated the matches set the graph pointers

ng = 0
do i = 1, np
  do j = 1, size(p(i)%p%graph)
    if (p(i)%p%graph(j)%name == graph_name .or. graph_name == '*') then
      ng = ng + 1
      g(ng)%g => p(i)%p%graph(j)
    endif
  enddo
enddo

if (present(graph)) graph = g

! Find number of curves that match

if (.not. present(curve)) return

if (logic_option(.false., always_allocate) .and. curve_name == ' ') curve_name = '*'
if (curve_name == ' ') return

nc = 0
do j = 1, ng
  do k = 1, size(g(j)%g%curve)
    if (g(j)%g%curve(k)%name == curve_name .or. curve_name == '*') nc = nc + 1
  enddo
enddo

if (nc == 0) then
  if (logic_option(.true., print_flag)) call out_io (s_error$, r_name, &
                  'CURVE NOT FOUND: ' // name)
  err = .true.
  return
endif

allocate (c(nc))
if (present(curve)) allocate (curve(nc))

! Now that we have counted and allocated the matches set the curve pointers

nc = 0
do j = 1, ng
  do k = 1, size(g(j)%g%curve)
    if (g(j)%g%curve(k)%name == curve_name .or. curve_name == '*') then
      nc = nc + 1
      c(nc)%c => g(j)%g%curve(k)
    endif
  enddo
enddo

if (present(curve)) curve = c

end subroutine tao_find_plots

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

var%useit_plot = var%exists .and. var%good_plot .and. var%good_var .and. &
                 (var%good_user .or. .not. graph%draw_only_good_user_data_or_vars)

end subroutine tao_var_useit_plot_calc

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_evaluate_element_parameters (err, param_name, values, print_err, dflt_source, dflt_component)
!
! Routine to evaluate a lattice element parameter of the form 
!     <universe>@ele::{<class>}::<ele_name_or_num>[<parameter>]{|<component>}
! or to evaluate at the middle of the element
!     <universe>@ele_mid::{<class>}::<ele_name_or_num>[<parameter>]{|<component>}
! Note: size(values) can be zero without an error
! 
! Input:
!   param_name      -- Character(*): parameter name.
!   print_err       -- Logical: Print error message? 
!   dflt_source     -- Character(*): Default source
!   dflt_component  -- Character(*): Default component
!
! Output:
!   err       -- Logical: True if there is an error in syntax. False otherwise
!   values(:) -- Real(rp), allocatable: Array of datum valuse.
!-

subroutine tao_evaluate_element_parameters (err, param_name, values, print_err, dflt_source, dflt_component)

implicit none

type (tao_universe_struct), pointer :: u
type (lat_struct), pointer :: lat
type (ele_struct) ele3
type (coord_struct), pointer :: this_orb(:)
type (coord_struct) orb
type (branch_struct), pointer :: branch

character(*) param_name
character(*) dflt_source
character(*), optional :: dflt_component
character(60) name, class_ele, parameter, component
character(40) :: r_name = 'tao_evaluate_element_parameters'

real(rp), allocatable :: values(:)
real(rp) :: real_val

integer i, j, ix, num, ixe, ix1, ios, n_tot

logical err, valid, middle
logical :: print_err

!

call tao_pick_universe (param_name, name, scratch%this_u, err)
if (err) return

err = .true.

if (name(1:5) == 'ele::') then
  name = name(6:)  ! Strip off 'ele::'
  middle = .false.
elseif (name(1:9) == 'ele_mid::') then   
  name = name(10:)  ! Strip off 'ele_mid::'
  middle = .true.
elseif (dflt_source /= 'element') then
  return
endif

! Get component

ix = index(name, '|')
if (ix == 0) then
  component = 'model'
  if (present(dflt_component)) then
    if (dflt_component /= '') component = dflt_component
  endif
else
  component = name(ix+1:)
  name = name(1:ix-1)
endif

! Get class:name

ix1 = index(name, '[');  
if (ix1 == 0) then
  return

else
  ix1 = index(name, '[');  if (ix1 == 0) return
  class_ele = name(1:ix1-1)
  name = name(ix1+1:)
  if (class_ele(1:2) == '::') class_ele = class_ele(3:)
  ix1 = index(name, ']');  if (ix1 == 0) return
  parameter = name(1:ix1-1)
endif

! Evaluate

n_tot = 0
do i = lbound(s%u, 1), ubound(s%u, 1)
  if (.not. scratch%this_u(i)) cycle
  u => s%u(i)
  call tao_locate_elements (class_ele, u%ix_uni, scratch%eles, err)
  if (err) return
  call re_allocate (values, n_tot + size(scratch%eles))

  do j = 1, size(scratch%eles)
    ixe = scratch%eles(j)%ele%ix_ele

    if (parameter == 'index') then
      values(n_tot+j) = ixe
      cycle
    endif

    select case (component)
    case ('model')   
      lat => u%model%lat
      this_orb => u%model%lat_branch(scratch%eles(j)%ele%ix_branch)%orbit
      branch => u%model%lat%branch(scratch%eles(j)%ele%ix_branch)
    case ('base')  
      lat => u%base%lat
      this_orb => u%base%lat_branch(scratch%eles(j)%ele%ix_branch)%orbit
      branch => u%base%lat%branch(scratch%eles(j)%ele%ix_branch)
    case ('design')
      lat => u%design%lat
      this_orb => u%design%lat_branch(scratch%eles(j)%ele%ix_branch)%orbit
      branch => u%design%lat%branch(scratch%eles(j)%ele%ix_branch)
    case default
      call out_io (s_error$, r_name, 'BAD DATUM COMPONENT FOR: ' // param_name)
      return
    end select

    if (middle .and. ixe /= 0) then
      call twiss_and_track_intra_ele (branch%ele(ixe), lat%param, 0.0_rp, branch%ele(ixe)%value(l$)/2, &
                .true., .false., this_orb(ixe-1), orb, branch%ele(ixe-1), ele3, err)
      call tao_orbit_value (parameter, orb, values(n_tot+j), err)
      if (err) then
        call ele_attribute_value (ele3, parameter, .true., values(n_tot+j), err, print_err)
        if (err) return
      endif

    else
      call tao_orbit_value (parameter, this_orb(ixe), values(n_tot+j), err)
      if (err) then
        call ele_attribute_value (branch%ele(ixe), parameter, .true., values(n_tot+j), err, print_err)
        if (err) return
      endif
    endif

  enddo

  n_tot = n_tot + size(values)
enddo

err = .false.

end subroutine tao_evaluate_element_parameters

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_orbit_value (component, orbit, value, err)
!
! Routine to return the orbit component indicated by component
!
! Input:
!   component -- character(*): 'orbit.x', 'orbit.px', 'intensity.x', 'phase.y', 'energy', 'pc', etc.
!   orbit     -- coord_struct: Particle orbit.
!
! Output:
!   value -- real(rp): orbit component.
!   err   -- logical: Set True if component is not recognized. False otherwise.
!-

subroutine tao_orbit_value (component, orbit, value, err)

implicit none

type (coord_struct) orbit

real(rp) value
character(*) component
logical err

!

err = .true.

if (component == 'state') then
  value = orbit%state
  err = .false.
  return
endif

if (orbit%state /= alive$) return

select case (component)
case ('orbit_x', 'orbit.x')
  value = orbit%vec(1)
case ('orbit_px', 'orbit.px')
  value = orbit%vec(2)
case ('orbit_y', 'orbit.y')
  value = orbit%vec(3)
case ('orbit_py', 'orbit.py')
  value = orbit%vec(4)
case ('orbit_z', 'orbit.z')
  value = orbit%vec(5)
case ('orbit_pz', 'orbit.pz')
  value = orbit%vec(6)
case ('intensity')
  value = orbit%field(1)**2 + orbit%field(2)**2
case ('intensity_x', 'intensity.x')
  value = orbit%field(1)**2
case ('intensity_y', 'intensity.y')
  value = orbit%field(2)**2
case ('phase_x', 'phase.x')
  value = orbit%phase(1)
case ('phase_y', 'phase.y')
  value = orbit%phase(2)
case ('energy')
  if (orbit%species == photon$) then
    value = orbit%p0c
  else
    call convert_pc_to(orbit%p0c * (1 + orbit%vec(6)), orbit%species, e_tot = value)
  endif
case ('pc')
  if (orbit%species == photon$) then
    value = orbit%p0c
  else
    value = orbit%p0c * (1 + orbit%vec(6))
  endif
case default
  return
end select

err = .false.

end subroutine tao_orbit_value

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_find_data (err, data_name, d2_array, d1_array, d_array, re_array, 
!                       log_array, str_array, int_array, ix_uni, dflt_index, print_err, component)
!
! Routine to set data pointers to the correct data structures. 
!
! The re_array will be used if the component is one of:
!   model, base, design, meas, ref, old, fit, weight
! The l_array will be used if the component is one of:
!   exists, good_meas, good_ref, good_user, good_opt, good_plot
! 
! Example:
!   data_name = '*@orbit.x'
! re_array & l_array will have size = 0 since there is no data component specified.
!
! Example:
!   data_name = 'orbit'
! In this case the default universe will be used. The d1_array will have two components
! pointing to orbit.x and orbit.y.
! re_array & l_array will have size = 0 since there is no data component specified.
!
! Example:
!   data_name = '2@orbit.x[3,7:9]|meas'
! The measured values for the 3rd, 7th, 8th and 9th elements of orbit.x in universe #2.
! r_arrray will be allocated and l_array will have size = 0.
!
! Example:
!   data_name = 'orbit.x'
!   dflt_index = '4'
! This is equivalent to:
!   data_name = 'orbit.x[4]'
! Notice that if dflt_index is not present, or is negative, 'orbit.x' will evaluate
! to an array of numbers.
!
! Input:
!   data_name    -- Character(*): The data name type. Eg: "3@orbit.x[2:5,10]|meas"
!   ix_uni       -- Integer, optional: Index of default universe to use.
!                     If ix_uni = 0 then "viewed" universe will be used.
!                     Also, if not present then the "viewed" universe will be used.
!   dflt_index   -- character, optional: If present and non-negative, and if no index is specified
!                     by the data_name argument, this index is used in the evaluation.
!   print_err    -- Logical, optional: Print error message if data is 
!                     not found? Default is True.
!
! Output:
!   err          -- Logical: Err condition
!   d2_array(:)  -- Tao_d2_data_array_struct, allocatable, optional: Array of pointers to all 
!                     the matching d2_data structure. Size(d2_array) = 0 if no structures found.
!   d1_array(:)  -- Tao_d1_data_array_struct, allocatable, optional: Array of pointers to all 
!                     the matching d1_data structures. Size(d1_array) = 0 if no structures found.
!   d_array(:)   -- Tao_data_array_struct, allocatable, optional: Array of pointers to all 
!                     the matching tao_data_structs.  Size(d_array) = 0 if no structures found.
!   re_array(:)  -- Tao_real_pointer_struct, allocatable, optional: Array of pointers to real 
!                     component values.  Size(re_array) = 0 if no structures found.
!   log_array(:) -- Tao_logical_array_struct, allocatable, optional: Array of pointers to
!                     logical component values.  Size(log_array) = 0 if no structures found.
!   str_array(:) -- Tao_string_array_struct, allocatable, optional: Array of pointers to 
!                     character component values.  Size(str_array) = 0 if no structures found.
!   int_array(:) -- Tao_integer_array_struct, allocatable, optional: Array of pointers to
!                     integer component values.  Size(int_array) = 0 if no structures found.
!   component    -- Character(*), optional: Name of the component. E.G: 'good_user'
!                     set to ' ' if no component present.
!-

subroutine tao_find_data (err, data_name, d2_array, d1_array, d_array, re_array, &
                           log_array, str_array, int_array, ix_uni, dflt_index, print_err, component)

implicit none

type (tao_d2_data_array_struct), allocatable, optional :: d2_array(:)
type (tao_d1_data_array_struct), allocatable, optional :: d1_array(:)
type (tao_data_array_struct), allocatable, optional    :: d_array(:)
type (tao_real_pointer_struct), allocatable, optional    :: re_array(:)
type (tao_integer_array_struct), allocatable, optional :: int_array(:)
type (tao_logical_array_struct), allocatable, optional :: log_array(:)
type (tao_string_array_struct), allocatable, optional  :: str_array(:)
type (tao_universe_struct), pointer :: u

character(*) :: data_name
character(*), optional :: component
character(*), optional :: dflt_index

character(20) :: r_name = 'tao_find_data'
character(80) dat_name, component_name
character(16), parameter :: real_components(9) = [ &
             'model    ', 'base     ', 'design   ', 'meas     ', 'ref      ', &
             'old      ', 'fit      ', 'weight   ', 'invalid  ']
character(16), parameter :: logic_components(9) = [ &
             'exists    ', 'good_meas ', 'good_ref  ', 'good_user ', 'good_opt  ', &
             'good_plot ', 'good_base ', 'useit_opt ', 'useit_plot']
character(16), parameter :: integer_components(5) = [ &
             'ix_ele      ', 'ix_ele_start', 'ix_ele_ref  ', &
             'ix_d1       ', 'ix_uni      ']
character(16), parameter :: string_components(4) = ['merit_type    ', 'ele_name      ', &
                                                    'ele_start_name', 'ele_ref_name  ']

integer, optional :: ix_uni
integer :: data_num, ios, n_found
integer i, ix, iu

logical err, component_here, this_err, print_error, error, found_data
logical, optional :: print_err

! Init

print_error = logic_option(.true., print_err)

if (present(d2_array)) then
  if (allocated (d2_array)) then
    if (size(d2_array) /= 0) deallocate (d2_array)
  endif
  if (.not. allocated(d2_array)) allocate (d2_array(0))
endif

if (present(d1_array)) then
  if (allocated (d1_array)) then
    if (size(d1_array) /= 0) deallocate (d1_array)
  endif
  if (.not. allocated(d1_array)) allocate (d1_array(0))
endif

if (present(d_array)) then
  if (allocated (d_array)) then
    if (size(d_array) /= 0) deallocate (d_array)
  endif
  if (.not. allocated(d_array)) allocate (d_array(0))
endif

if (present(re_array)) then
  if (allocated (re_array)) then
    if (size(re_array) /= 0) deallocate (re_array)
  endif
  if (.not. allocated(re_array)) allocate (re_array(0))
endif

if (present(int_array)) then
  if (allocated (int_array)) then
    if (size(int_array) /= 0) deallocate (int_array)
  endif
  if (.not. allocated(int_array)) allocate (int_array(0))
endif

if (present(log_array)) then
  if (allocated (log_array)) then
    if (size(log_array) /= 0) deallocate (log_array)
  endif
  if (.not. allocated(log_array)) allocate (log_array(0))
endif

if (present(str_array)) then
  if (allocated (str_array)) then
    if (size(str_array) /= 0) deallocate (str_array)
  endif
  if (.not. allocated(str_array)) allocate (str_array(0))
endif

err = .true.
found_data = .false.

if (data_name == '') then
  if (print_error) call out_io (s_error$, r_name, 'DATA NAME IS BLANK')
  return
endif

! Select meas, ref, etc.

ix = index(data_name, '|')
if (ix == 0) then  ! not present
  component_here = .false.
  component_name = ' '
  dat_name = data_name
else
  component_here = .true.
  component_name = data_name(ix+1:)
  dat_name = data_name(:ix-1)
endif
if (present(component)) component = component_name

call string_trim (dat_name, dat_name, ix)

if (component_here) then
  call string_trim (component_name, component_name, ix)
  if (.not. any(component_name == real_components) .and. &
      .not. any(component_name == logic_components) .and. &
      .not. any(component_name == integer_components) .and. &
      .not. any(component_name == string_components)) then
    if (print_error) call out_io (s_error$, r_name, "BAD COMPONENT NAME: " // data_name)
    return            
  endif
endif

! Select universe

call tao_pick_universe (dat_name, dat_name, scratch%picked, this_err)
if (this_err) return

! Trim 'data::' suffix if present

if (index(dat_name, 'dat::') /= 0) then
  call out_io (s_error$, r_name, 'NAME USES OLD "dat::" SYNTAX. PLEASE CHANGE TO "data::": ' // dat_name)
  call err_exit
endif

if (dat_name(1:6) == 'data::') dat_name = dat_name(7:)

! Find the d2 data.

if (present(ix_uni)) then
  u => tao_pointer_to_universe (ix_uni)
  if (.not. associated(u)) return
  call find_this_d2 (u, dat_name, this_err)
else
  do i = lbound(s%u, 1), ubound(s%u, 1)
    if (.not. scratch%picked(i)) cycle
    u => tao_pointer_to_universe (i)
    call find_this_d2 (u, dat_name, this_err)
    if (this_err) return
  enddo
endif

! error check

if (this_err .or. .not. found_data) then
  if (print_error) call out_io (s_error$, r_name, "Couldn't find data: " // data_name)
  return
endif

err = .false.

!----------------------------------------------------------------------------
contains

subroutine find_this_d2 (uu, name, this_err)

type (tao_universe_struct), target :: uu
integer i, ix
character(*) name
character(80) d1_name, d2_name
logical this_err

! Everything before a period is the d2 name.
! if no period then must be something like name = "orbit" and everything is the d2 name.

this_err = .false.

ix = index(name, '.')
if (ix == 0) then
  ix = index(name, '[')
  if (ix /= 0) then
    d2_name = name(1:ix-1)
    d1_name = name(ix:)
  else
    d2_name = name
    d1_name = '*'
    if (present(dflt_index)) d1_name = dflt_index
  endif
else
  d2_name = name(1:ix-1)
  d1_name = name(ix+1:)
endif

! loop over matching d2 names

do i = 1, uu%n_d2_data_used
  if (.not. match_wild(uu%d2_data(i)%name, d2_name)) cycle
  call find_this_d1 (uu%d2_data(i), d1_name, .false., this_err)
  if (this_err) return
enddo

end subroutine find_this_d2

!----------------------------------------------------------------------------
! contains

subroutine find_this_d1 (d2, name, found_d1, this_err)

type (tao_d2_data_struct), target :: d2
type (tao_d2_data_array_struct), allocatable :: d2_temp(:)
integer i, ix, nd
character(*) name
character(80) d1_name, d_name
logical found_d1, this_err

! d2_array

if (present(d2_array)) then
  if (allocated(d2_array)) then
    nd = size(d2_array)
    call move_alloc(d2_array, d2_temp)
    allocate (d2_array(nd+1))
    d2_array(1:nd) = d2_temp
    deallocate(d2_temp)
    d2_array(nd+1)%d2 => d2
  else
    allocate (d2_array(1))
    d2_array(1)%d2 => d2
  endif
endif

! Everything before a '[' is the d1 name.

ix = index(name, '[')

if (ix == 0) then
  d1_name = name
  d_name = '*'
  if (present(dflt_index)) d_name = dflt_index
else
  d1_name = name(1:ix-1)
  d_name = name(ix+1:)
  ix = index(d_name, ']')
  if (ix == 0) then
    if (print_error) call out_io (s_error$, r_name, "NO MATCHING ']': " // data_name)
    this_err = .true.
    return
  endif
  if (d_name(ix+1:) /= ' ') then
    if (print_error) call out_io (s_error$, r_name, "GARBAGE AFTER ']': " // data_name)
    this_err = .true.
    return
  endif
  d_name = d_name(:ix-1)
endif

if (size(d2%d1) == 1 .and. d1_name == '') d1_name = '*'

do i = 1, size(d2%d1)
  if (.not. match_wild(d2%d1(i)%name, d1_name)) cycle
  call find_this_data (d2%d1(i), d_name, this_err)
  if (this_err) return
enddo

end subroutine find_this_d1

!----------------------------------------------------------------------------
! contains

subroutine find_this_data (d1, name, this_err)

type (tao_d1_data_struct), target :: d1
type (tao_d1_data_array_struct), allocatable :: d1_temp(:)
type (tao_data_array_struct), allocatable, save :: da(:)
type (tao_real_pointer_struct), allocatable, save :: ra(:)
type (tao_integer_array_struct), allocatable, save :: ia(:)
type (tao_logical_array_struct), allocatable, save :: la(:)
type (tao_string_array_struct), allocatable, save  :: sa(:)

integer i, j, nd, nl, i1, i2, num

character(*) name
character(80) d1_name, d_name

logical this_err
logical, allocatable, save :: list(:)

! d1_array

if (present(d1_array)) then
  if (allocated(d1_array)) then
    nd = size(d1_array)
    call move_alloc(d1_array, d1_temp)
    allocate (d1_array(nd+1))
    d1_array(1:nd) = d1_temp
    deallocate(d1_temp)
    d1_array(nd+1)%d1 => d1
  else
    allocate (d1_array(1))
    d1_array(1)%d1 => d1
  endif
endif

!

if (allocated(list)) deallocate(list)
i1 = lbound(d1%d, 1)
i2 = ubound(d1%d, 1)
allocate (list(i1:i2))

if (name == '*') then
  list = .true.

else
  call location_decode (name, list, i1, num)
  if (num <  0) then
    call out_io (s_error$, r_name, "BAD DATA INDEX NUMBER(S): " // name)
    this_err = .true.
    return  
  endif
endif

err = .false.
nl = count(list)
if (nl > 0) found_data = .true.

! data array

if (present(d_array)) then

  if (allocated(d_array)) then
    nd = size(d_array)
    allocate (da(nd))
    da = d_array
    deallocate(d_array)
    allocate (d_array(nl+nd))
    j = nd
    d_array(1:nd) = da
    deallocate(da)
  else
    allocate (d_array(nl))
    j = 0
  endif

  do i = i1, i2
    if (list(i)) then
      j = j + 1
      d_array(j)%d => d1%d(i)
    endif
  enddo

endif

! Integer component array

if (present(int_array) .and.  any(component_name == integer_components)) then

  if (allocated(int_array)) then
    nd = size(int_array)
    allocate (ia(nd))
    ia = int_array
    deallocate(int_array)
    allocate (int_array(nl+nd))
    j = nd
    int_array(1:nd) = ia
    deallocate(ia)
  else
    allocate (int_array(nl))
    j = 0
  endif

  do i = i1, i2
    if (list(i)) then
      j = j + 1
      select case (component_name)
      case ('ix_ele')
        int_array(j)%i => d1%d(i)%ix_ele
      case ('ix_ele_start')
        int_array(j)%i => d1%d(i)%ix_ele_start
      case ('ix_ele_ref')
        int_array(j)%i => d1%d(i)%ix_ele_ref
      case ('ix_d1')
        int_array(j)%i => d1%d(i)%ix_d1
      case ('ix_uni')
        int_array(j)%i => d1%d(i)%d1%d2%ix_uni
      case default
        call out_io (s_fatal$, r_name, "INTERNAL ERROR: INTEGER DATA")
        call err_exit
      end select
    endif
  enddo

endif

! real component array

if (present(re_array) .and.  any(component_name == real_components)) then

  if (allocated(re_array)) then
    nd = size(re_array)
    allocate (ra(nd))
    ra = re_array
    deallocate(re_array)
    allocate (re_array(nl+nd))
    j = nd
    re_array(1:nd) = ra
    deallocate(ra)
  else
    allocate (re_array(nl))
    j = 0
  endif

  do i = i1, i2
    if (list(i)) then
      j = j + 1
      select case (component_name)
      case ('model')
        re_array(j)%r => d1%d(i)%model_value
        re_array(j)%good_user  => d1%d(i)%good_user
        re_array(j)%good_value => d1%d(i)%good_model
      case ('base')
        re_array(j)%r => d1%d(i)%base_value
        re_array(j)%good_user  => d1%d(i)%good_user
        re_array(j)%good_value => d1%d(i)%good_base
      case ('design')
        re_array(j)%r => d1%d(i)%design_value
        re_array(j)%good_user  => d1%d(i)%good_user
        re_array(j)%good_value => d1%d(i)%good_model
      case ('meas')
        re_array(j)%r => d1%d(i)%meas_value
        re_array(j)%good_user  => d1%d(i)%good_user
        re_array(j)%good_value => d1%d(i)%good_meas
      case ('ref')
        re_array(j)%r => d1%d(i)%ref_value
        re_array(j)%good_user  => d1%d(i)%good_user
        re_array(j)%good_value => d1%d(i)%good_ref
      case ('old')
        re_array(j)%r => d1%d(i)%old_value
        re_array(j)%good_user  => d1%d(i)%good_user
        re_array(j)%good_value => forever_true$
      case ('invalid')
        re_array(j)%r => d1%d(i)%invalid_value
        re_array(j)%good_user  => d1%d(i)%good_user
        re_array(j)%good_value => forever_true$
      case ('weight')
        re_array(j)%r => d1%d(i)%weight
        re_array(j)%good_user  => d1%d(i)%good_user
        re_array(j)%good_value => forever_true$
      case default
        call out_io (s_fatal$, r_name, "INTERNAL ERROR: REAL DATA")
        call err_exit
      end select
    endif
  enddo

endif

! logical component array

if (present(log_array) .and. any(component_name == logic_components)) then

  if (allocated(log_array) .and. component_here) then
    nd = size(log_array)
    allocate (la(nd))
    la = log_array
    deallocate(log_array)
    allocate (log_array(nl+nd))
    j = nd
    log_array(1:nd) = la
    deallocate(la)
  else
    allocate (log_array(nl))
    j = 0
  endif

  do i = i1, i2
    if (list(i)) then
      j = j + 1
      select case (component_name)
      case ('exists')
        log_array(j)%l => d1%d(i)%exists
      case ('good_base')
        log_array(j)%l => d1%d(i)%good_base
      case ('good_meas')
        log_array(j)%l => d1%d(i)%good_meas
      case ('good_ref')
        log_array(j)%l => d1%d(i)%good_ref
      case ('good_user')
        log_array(j)%l => d1%d(i)%good_user
      case ('good_opt')
        log_array(j)%l => d1%d(i)%good_opt
      case ('good_plot')
        log_array(j)%l => d1%d(i)%good_plot
      case ('useit_opt')
        log_array(j)%l => d1%d(i)%useit_opt
      case ('useit_plot')
        log_array(j)%l => d1%d(i)%useit_plot
      case default
        call out_io (s_fatal$, r_name, "INTERNAL ERROR: LOGIC DATA")
        call err_exit
      end select
    endif
  enddo

endif

! string component array

if (present(str_array) .and. any(component_name == string_components)) then

  if (allocated(str_array) .and. component_here) then
    nd = size(str_array)
    allocate (sa(nd))
    sa = str_array
    deallocate(str_array)
    allocate (str_array(nl+nd))
    j = nd
    str_array(1:nd) = sa
    deallocate(sa)
  else
    allocate (str_array(nl))
    j = 0
  endif

  do i = i1, i2
    if (list(i)) then
      j = j + 1
      select case (component_name)
      case ('merit_type')
        str_array(j)%s => d1%d(i)%merit_type
      case ('ele_name')
        str_array(j)%s => d1%d(i)%ele_name
      case ('ele_start_name')
        str_array(j)%s => d1%d(i)%ele_start_name
      case ('ele_ref_name')
        str_array(j)%s => d1%d(i)%ele_ref_name
      case default
        call out_io (s_fatal$, r_name, "INTERNAL ERROR: STRING DATA")
        call err_exit
      end select
    endif
  enddo

endif

end subroutine find_this_data

end subroutine tao_find_data

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_find_var (err, var_name, v1_array, v_array, re_array, log_array,  
!                                                    str_array, print_err, component)
!
! Find a v1 variable type, and variable component then point to it.
!
! The re_array will be used if the component is one of:
!   model, base, design, meas, ref, old, step, weight, high_lim, low_lim 
! The log_array will be used if the component is one of:
!   exists, good_var, good_user, good_opt, good_plot
! 
! Note: Any of the output allocatable arrays will have size = 0 if not used.
! For example, if var_name is "model", then str_array, if present, have size = 0.
!
! Example:
!   var_name = 'quad_k1[3]|design'
!
! Input:
!   var_name     -- Character(*): Name of the variable.
!   print_err    -- Logical, optional: Print error message if data is 
!                     not found? Default is True.
!
! Output:
!   err          -- Logical: err condition
!   v1_array(:)  -- Tao_v1_var_array_struct, allocatable, optional: Array of pointers to 
!                     all the v1_var structures.
!   v_array(:)   -- Tao_var_array_struct, allocatable, optional: Array of pointers to the 
!                     variable data point.
!   re_array(:)  -- Tao_real_pointer_struct, allocatable, optional: Array of pointers to 
!                     the real component values.
!   log_array(:) -- Tao_logical_array_struct, allocatable, optional: Array of pointers to
!                     logical component values.
!   str_array(:) -- Tao_string_array_struct, allocatable, optional: Array of pointers to 
!                     character component values.
!   component    -- Character(*), optional: Name of the component. E.G: 'good_user'
!                   set to ' ' if no component present.
!-

subroutine tao_find_var (err, var_name, v1_array, v_array, re_array, log_array, &
                                                    str_array, print_err, component)

implicit none

type (tao_v1_var_array_struct), allocatable, optional  :: v1_array(:)
type (tao_var_array_struct), allocatable, optional     :: v_array(:)
type (tao_real_pointer_struct), allocatable, optional  :: re_array(:)
type (tao_logical_array_struct), allocatable, optional :: log_array(:)
type (tao_string_array_struct), allocatable, optional  :: str_array(:)

integer i, ix, n_var, ios

character(16), parameter :: real_components(11) = [ &
             'model    ', 'base     ', 'design   ', 'meas     ', 'ref      ', &
             'old      ', 'step     ', 'weight   ', 'high_lim ', 'low_lim  ', 'key_delta']
character(16), parameter :: logic_components(8) = [ &
             'exists    ', 'good_var  ', 'good_user ', 'good_opt  ', 'good_plot ', &
             'useit_opt ', 'useit_plot', 'key_bound ']
character(16), parameter :: string_components(1) = ['merit_type']

character(*) :: var_name
character(*), optional :: component
character(20) :: r_name = 'tao_find_var'
character(80) v1_name, v_name, component_name

logical, optional :: print_err
logical err, component_here, this_err, print_error

! Init

print_error = logic_option(.true., print_err)

if (present(v1_array)) then
  if (allocated (v1_array)) then
    if (size(v1_array) /= 0) deallocate (v1_array)
  endif
  if (.not. allocated(v1_array)) allocate (v1_array(0))
endif

if (present(v_array)) then
  if (allocated (v_array)) then
    if (size(v_array) /= 0) deallocate (v_array)
  endif
  if (.not. allocated(v_array)) allocate (v_array(0))
endif

if (present(re_array)) then
  if (allocated (re_array)) then
    if (size(re_array) /= 0) deallocate (re_array)
  endif
  if (.not. allocated(re_array)) allocate (re_array(0))
endif

if (present(log_array)) then
  if (allocated (log_array)) then
    if (size(log_array) /= 0) deallocate (log_array)
  endif
  if (.not. allocated(log_array)) allocate (log_array(0))
endif

if (present(str_array)) then
  if (allocated (str_array)) then
    if (size(str_array) /= 0) deallocate (str_array)
  endif
  if (.not. allocated(str_array)) allocate (str_array(0))
endif

err = .true.

! Error if no variables exist

if (s%n_var_used == 0) then
  if (print_error) call out_io (s_warn$, r_name, &
                        "NO VARIABLES HAVE BEEN DEFINED IN THE INPUT FILES!")
  return
endif

! Select meas, ref, etc.

ix = index(var_name, '|')
if (ix == 0) then  ! not present
  component_here = .false.
  component_name = ' '   ! garbage
  v1_name = var_name
else
  component_here = .true.
  component_name = var_name(ix+1:)
  v1_name = var_name(:ix-1)
endif
if (present(component)) component = component_name

call string_trim (v1_name, v1_name, ix)
call string_trim (component_name, component_name, ix)

if (component_here) then
  if (.not. any(component_name == real_components) .and. &
      .not. any(component_name == logic_components) .and. &
      .not. any(component_name == string_components)) then
    if (print_error) call out_io (s_error$, r_name, "BAD COMPONENT NAME: " // var_name)
    return            
  endif
endif

! Trim 'var::' suffix if present

if (v1_name(1:5) == 'var::') v1_name = v1_name(6:)

! split on '['

ix = index(v1_name, '[')
if (ix == 0) then
  v_name = '*'
else
  v_name  = v1_name(ix+1:)
  v1_name = v1_name(1:ix-1)
  ix = index(v_name, ']')
  if (ix == 0) then
    if (print_error) call out_io (s_error$, r_name, "NO MATCHING ']': " // var_name)
    return
  endif
  if (v_name(ix+1:) /= ' ') then
    if (print_error) call out_io (s_error$, r_name, "GARBAGE AFTER ']': " // var_name)
    return
  endif
  v_name = v_name(:ix-1)
endif

call string_trim(v1_name, v1_name, ix)
if (ix == 0) then
  if (print_error) call out_io (s_error$, r_name, 'VARIABLE NAME IS BLANK')
  return
endif

! Point to the correct v1 var type 

if (v1_name == '*') then
  if (present(v1_array)) then
    deallocate (v1_array)
    allocate (v1_array(s%n_v1_var_used))
  endif
  do i = 1, s%n_v1_var_used
    if (present(v1_array)) v1_array(i)%v1 => s%v1_var(i)
    call find_this_var (s%v1_var(i), v_name, this_err)
    if (this_err) return
  enddo

else
  do i = 1, s%n_v1_var_used
    if (v1_name == s%v1_var(i)%name) then
      if (present(v1_array)) then
        deallocate (v1_array)
        allocate (v1_array(1))
        v1_array(1)%v1 => s%v1_var(i)
      endif
      call find_this_var (s%v1_var(i), v_name, this_err)
      exit
    endif
  enddo
endif

! error check

if (err) then
  if (print_error) call out_io (s_error$, r_name, "COULDN'T FIND VARIABLE: " // var_name)
  return
endif

!----------------------------------------------------------------------------
contains

subroutine find_this_var (v1, name, this_err)

type (tao_v1_var_struct) :: v1
type (tao_var_array_struct), allocatable, save :: va(:)
type (tao_real_pointer_struct), allocatable, save :: ra(:)
type (tao_logical_array_struct), allocatable, save :: la(:)
type (tao_string_array_struct), allocatable, save  :: sa(:)

integer i, j, nd, nl, i1, i2, num

character(*) name
character(40), allocatable, save :: names(:)
character(80) v1_name, v_name

logical this_err
logical, allocatable, save :: list(:)

!

if (allocated(list)) deallocate(list, names)
i1 = lbound(v1%v, 1)
i2 = ubound(v1%v, 1)
allocate (list(i1:i2), names(i1:i2))
this_err = .false.

if (name == '*') then
  list = .true.

else
  do i = i1, i2
    names(i) = v1%v(i)%ele_name
  enddo
  call location_decode (name, list, i1, num, names)
  if (num <  0) then
    call out_io (s_error$, r_name, "BAD VAR INDEX NUMBER(S): " // name)
    this_err = .true.
    return  
  endif
endif

err = .false.
nl = count(list)

! real array

if (present(v_array)) then

  if (allocated(v_array)) then
    nd = size(v_array)
    allocate (va(nd))
    va = v_array
    deallocate(v_array)
    allocate (v_array(nl+nd))
    j = nd
    v_array(1:nd) = va
    deallocate(va)
  else
    allocate (v_array(nl))
    j = 0
  endif

  do i = i1, i2
    if (list(i)) then
      j = j + 1
      v_array(j)%v => v1%v(i)
    endif
  enddo

endif

! real component array

if (present(re_array) .and.  any(component_name == real_components)) then

  if (allocated(re_array)) then
    nd = size(re_array)
    allocate (ra(nd))
    ra = re_array
    deallocate(re_array)
    allocate (re_array(nl+nd))
    j = nd
    re_array(1:nd) = ra
    deallocate(ra)
  else
    allocate (re_array(nl))
    j = 0
  endif

  do i = i1, i2
    if (list(i)) then
      j = j + 1
      select case (component_name)
      case ('model')
        re_array(j)%r => v1%v(i)%model_value
      case ('base')
        re_array(j)%r => v1%v(i)%base_value
      case ('design')
        re_array(j)%r => v1%v(i)%design_value
      case ('meas')
        re_array(j)%r => v1%v(i)%meas_value
      case ('ref')
        re_array(j)%r => v1%v(i)%ref_value
      case ('old')
        re_array(j)%r => v1%v(i)%old_value
      case ('step')
        re_array(j)%r => v1%v(i)%step
      case ('weight')
        re_array(j)%r => v1%v(i)%weight
      case ('high_lim')
        re_array(j)%r => v1%v(i)%high_lim
      case ('low_lim')
        re_array(j)%r => v1%v(i)%low_lim
      case ('key_delta')
        re_array(j)%r => v1%v(i)%key_delta
      case default
        call out_io (s_fatal$, r_name, "INTERNAL ERROR: REAL VAR")
        call err_exit
      end select
    endif
  enddo

endif

! logical component array

if (present(log_array) .and. any(component_name == logic_components)) then

  if (allocated(log_array) .and. component_here) then
    nd = size(log_array)
    allocate (la(nd))
    la = log_array
    deallocate(log_array)
    allocate (log_array(nl+nd))
    j = nd
    log_array(1:nd) = la
    deallocate(la)
  else
    allocate (log_array(nl))
    j = 0
  endif

  do i = i1, i2
    if (list(i)) then
      j = j + 1
      select case (component_name)
      case ('exists')
        log_array(j)%l => v1%v(i)%exists
      case ('good_var')
        log_array(j)%l => v1%v(i)%good_var
      case ('good_user')
        log_array(j)%l => v1%v(i)%good_user
      case ('good_opt')
        log_array(j)%l => v1%v(i)%good_opt
      case ('good_plot')
        log_array(j)%l => v1%v(i)%good_plot
      case ('useit_opt')
        log_array(j)%l => v1%v(i)%useit_opt
      case ('useit_plot')
        log_array(j)%l => v1%v(i)%useit_plot
      case ('key_bound')
        log_array(j)%l => v1%v(i)%key_bound
      case default
        call out_io (s_fatal$, r_name, "INTERNAL ERROR: LOGIC VAR")
        call err_exit
      end select
    endif
  enddo

endif

! string component array

if (present(str_array) .and. any(component_name == string_components)) then

  if (allocated(str_array) .and. component_here) then
    nd = size(str_array)
    allocate (sa(nd))
    sa = str_array
    deallocate(str_array)
    allocate (str_array(nl+nd))
    j = nd
    str_array(1:nd) = sa
    deallocate(sa)
  else
    allocate (str_array(nl))
    j = 0
  endif

  do i = i1, i2
    if (list(i)) then
      j = j + 1
      select case (component_name)
      case ('merit_type')
        str_array(j)%s => v1%v(i)%merit_type
      case default
        call out_io (s_fatal$, r_name, "INTERNAL ERROR: STRING VAR")
        call err_exit
      end select
    endif
  enddo

endif

end subroutine find_this_var

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

do j = 1, s%n_var_used
  var => s%var(j)
  var%correction_value = var%meas_value + (var%design_value - var%model_value)
enddo

end subroutine tao_var_target_calc

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_set_var_model_value (var, value, print_limit_warning)
!
! Subroutine to set the value for a model variable and do the necessary bookkeeping.
! If value is past the variable's limit, and s%global%var_limits_on = True, the
! variable will be set to the limit value.
!
! Input:
!   var   -- Tao_var_struct: Variable to set
!   value -- Real(rp): Value to set to
!   print_limit_warning
!         -- Logical, optional: Print a warning if the value is past the variable's limits.
!             Default is True.
!-

subroutine tao_set_var_model_value (var, value, print_limit_warning)

implicit none

type (tao_var_struct), target :: var
type (tao_universe_struct), pointer :: u
type (lat_struct), pointer :: lat
type (tao_this_var_struct), pointer :: t
type (ele_struct), pointer :: ele

real(rp) value
integer i
logical, optional :: print_limit_warning

!

if (.not. var%exists) return

! check if hit variable limit
if (s%global%var_limits_on .and. (.not. s%global%var_limit_only_used .or. var%useit_opt)) then
  if (value < var%low_lim) then
    if (logic_option (.true., print_limit_warning)) &
          call out_io (s_blank$, ' ', "Hit lower limit of variable: " // tao_var1_name(var))
    value = var%low_lim
  elseif (value > var%high_lim) then
    if (logic_option (.true., print_limit_warning)) &
          call out_io (s_blank$, ' ', "Hit upper limit of variable: " // tao_var1_name(var))
    value = var%high_lim
  endif
endif

var%model_value = value
do i = 1, size(var%this)
  t => var%this(i)
  t%model_value = value
  u => s%u(t%ix_uni)

  if (s%com%common_lattice .and.  t%ix_uni == ix_common_uni$) then
    s%u(:)%calc%lattice = .true.
  else
    u%calc%lattice = .true.
  endif

  lat => u%model%lat
  if (var%ele_name == 'BEAM_START') then
    u%model%lat_branch(0)%orb0%vec = lat%beam_start%vec
    u%model%lat_branch(0)%orb0%t   = lat%beam_start%t
    u%model%lat_branch(0)%orb0%p0c = lat%beam_start%p0c
  else
    ele => lat%branch(t%ix_branch)%ele(t%ix_ele)
    call set_flags_for_changed_attribute (ele, t%model_value)
  endif
enddo

end subroutine tao_set_var_model_value

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
! Subroutine tao_lat_bookkeeper (u, tao_lat)
!
! This will make sure all bookkeeping is up to date.
!
! Input:
!  u            -- tao_universe_struct
!  lat_name     -- Integer: Which lattice
!
! Output:
!  tao_lat      -- lat_struct
!-

subroutine tao_lat_bookkeeper (u, tao_lat)

implicit none

type (tao_universe_struct), target :: u
type (tao_lattice_struct) :: tao_lat

integer i, j

character(20) :: r_name = "tao_lat_bookkeeper"

! Setup from common if it exists

if (associated(u%common)) then

  ! First put in the common values

  do i = 1, s%n_var_used
    do j = 1, size(s%var(i)%this)
      s%var(i)%this(j)%model_value = s%var(i)%common%model_value
      s%var(i)%this(j)%base_value  = s%var(i)%common%base_value
    enddo
  enddo

  ! Then put in the values for this universe

  do i = 1, s%n_var_used
    do j = 1, size(s%var(i)%this)
      if (s%var(i)%this(j)%ix_uni /= u%ix_uni) cycle
      s%var(i)%this(j)%model_value = s%var(i)%model_value
      s%var(i)%this(j)%base_value = s%var(i)%base_value
    enddo
  enddo

endif

! Do bookkeeping

call lattice_bookkeeper (tao_lat%lat)

end subroutine tao_lat_bookkeeper

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Function tao_read_this_index (name, ixc) result (ix)
!
! Returns the integer value in the array <name> at position <ixc>. This is used
! for finding a 6-dimensional index reference so any value at <ixc> greater than
! 6 returns an error.
!
! Input:
!  name     -- Character(*): character array holding the index
!  ixc      -- Integer: location within <name> to evaluate
!
! Output:
!  ix       -- Integer: Index at <name>(<ixc>:<ixc>)
!
! Example:
!      name = r:26
!      ixc  = 3
!
! Gives:
!      ix = 2
!
! Example:
!      name = mat_94
!      ixc  = 7
! Gives:
!      Error: "BAD INDEX CONSTRAINT: mat_94"
!-

function tao_read_this_index (name, ixc) result (ix)

character(*) name
integer ix, ixc
character(20) :: r_name = 'tao_read_this_index'

ix = index('123456', name(ixc:ixc))
if (ix == 0) then
  call out_io (s_abort$, r_name, 'BAD INDEX CONSTRAINT: ' // name)
  call err_exit
endif

end function tao_read_this_index


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Function tao_var1_name(var) result (var1_name)
!
! Function to return the variable name in the form:
!   var1_name[index]
! For example:
!   quad_k1[23]
!
! Input:
!   var -- Tao_var_struct: Variable
!
! Output:
!   var1_name -- Character(60): Appropriate name.
!-

function tao_var1_name(var) result (var1_name)

implicit none

type (tao_var_struct) var
character(60) var1_name

!

write (var1_name, '(2a, i0, a)') trim(var%v1%name), '[', var%ix_v1, ']'

end function tao_var1_name

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Function tao_var_attrib_name(var) result (var_attrib_name)
!
! Function to return the variable name in the form:
!   {universe@}element[attribute]
! For example:
!   Q03W[k1]
!   [2,3]@V11E[HKICK]
!
! Input:
!   var -- Tao_var_struct: Variable
!
! Output:
!   var_attrib_name -- Character(60): Appropriate name.
!-

function tao_var_attrib_name(var) result (var_attrib_name)

use location_encode_mod

implicit none

type (tao_var_struct) var

character(60) var_attrib_name
integer i

!

if (size(s%u) > 1 .and. size(var%this) > 0) then
  call location_encode (var_attrib_name, var%this%ix_uni, ',')
  if (index(var_attrib_name, ',') /= 0 .or. index(var_attrib_name, ':') /= 0) then  
    var_attrib_name = '[' // trim(var_attrib_name) // ']' 
  endif
  var_attrib_name = trim(var_attrib_name) // '@' // &
                    trim(var%ele_name) // '[' // trim(var%attrib_name) // ']'
else
  var_attrib_name = trim(var%ele_name) // '[' // trim(var%attrib_name) // ']'
endif

end function tao_var_attrib_name

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Function tao_constraint_type_name (datum) result (datum_name)
!
! Function to return the constraint type in the form:
!   data_type <merit_type>
! For example:
!   eta.x <target>
!
! Input:
!   datum      -- Tao_data_struct: Datum
!
! Output:
!   datum_name -- Character(60): Appropriate name.
!-

function tao_constraint_type_name(datum) result (datum_name)

implicit none

type (tao_data_struct) datum
character(200) datum_name
integer ix

! Expressions are too long so shorten the name

datum_name = trim(datum%data_type) // ' <' // trim(datum%merit_type) // '>'
if (datum_name(1:11) == 'expression:') call string_trim (datum_name(12:), datum_name, ix)

end function

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Function tao_datum_name (datum, show_universe) result (datum_name)
!
! Function to return the datum name in the form:
!   d2_name.d1_name[index]
! or (if show_universe is True and there is more than one universe):
!   universe@d2_name.d1_name[index]
! For example:
!   2@orbit.x[23]
!   
!
! Input:
!   datum         -- Tao_data_struct: Datum
!   show_universe -- Logical, optional: Show the datum's universe.
!                       Default is True.
!
! Output:
!   datum_name -- Character(60): Appropriate name.
!-

function tao_datum_name(datum, show_universe) result (datum_name)

implicit none

type (tao_data_struct) datum
character(60) datum_name
logical, optional :: show_universe

! If this datum is "isolated". That is, it does not have an associated d1_data 
! structure then just use it's data_type.
! This can happen if the datum is derived from a curve.

if (.not. associated(datum%d1)) then
  datum_name = datum%data_type
  return
endif

! Normal case

datum_name = tao_d2_d1_name (datum%d1, show_universe)
write (datum_name, '(2a, i0, a)') trim(datum_name), '[', datum%ix_d1, ']'

end function tao_datum_name

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Function tao_d2_d1_name (d1, show_universe) result (d2_d1_name)
!
! Function to return the datum name in the form:
!   d2_name.d1_name
! If there is only one d1_data array associated with the d2_data 
! array then the name is shortened to:
!   d2_name
! Additionally, if show_universe is True and there is more than one universe
! then "universe@" is prepended to the name.
! For example:
!   2@orbit.x
!   
!
! Input:
!   d1            -- Tao_d1_data_struct: Data array.
!   show_universe -- Logical, optional: Show the datum's universe.
!                       Default is True.
!
! Output:
!   d2_d1_name -- Character(60): Appropriate name.
!-

function tao_d2_d1_name(d1, show_universe) result (d2_d1_name)

implicit none

type (tao_d1_data_struct) d1
character(60) d2_d1_name
logical, optional :: show_universe

! If there is only one d1 array associated with the d2_data array then
! drop the d1 name.

if (size(d1%d2%d1) == 1) then
  write (d2_d1_name, '(a)') trim(d1%d2%name)
else
  write (d2_d1_name, '(3a)') trim(d1%d2%name), '.', trim(d1%name)
endif

!

if (size(s%u) > 1 .and. logic_option(.true., show_universe)) &
       write (d2_d1_name, '(i0, 2a)') d1%d2%ix_uni, '@', trim(d2_d1_name)

end function tao_d2_d1_name

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Function tao_curve_name (curve, use_region) result (curve_name)
!
! Function to return the curve name in the form:
!   plot_name.graph_name.curve_name
! For example:
!   orbit.x.c1
!
! Input:
!   curve      -- Tao_curve_struct: Curve
!   use_region -- Logical: If present and True then use the region 
!                  name instead of the plot name. Region name is
!                  'NULL_REGION' if there is no assocaited region.
!
! Output:
!   curve_name -- Character(60): Appropriate name.
!-

function tao_curve_name(curve, use_region) result (curve_name)

implicit none

type (tao_curve_struct) curve
character(60) curve_name
logical, optional :: use_region

!

curve_name = '.' // trim(curve%g%name) // '.' // trim(curve%name)

if (logic_option(.false., use_region)) then
  if (associated(curve%g%p%r)) then
    curve_name = trim(curve%g%p%r%name) // trim(curve_name)
  else
    curve_name = 'NULL_REGION' // trim(curve_name)
  endif
else
    curve_name = trim(curve%g%p%name) // trim(curve_name)
endif

end function

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Function tao_universe_number (i_uni) result (i_this_uni)
!
! Fnction to return the universe number.
! i_this_uni = i_uni except when i_uni is -1. 
! In this case i_this_uni = s%com%default_universe.
!
! Input:
!   i_uni -- Integer: Nominal universe number.
!
! Output:
!   i_this_uni -- Integer: Universe number. 
!-

function tao_universe_number (i_uni) result (i_this_uni)

implicit none

integer i_uni, i_this_uni

i_this_uni = i_uni
if (i_uni == -1) i_this_uni = s%com%default_universe

end function tao_universe_number

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine tao_parse_command_args (error, cmd_words)
!
! Subroutine to parse the command line arguments.
!
! Input:
!   cmd_words(:) -- Character(*), optional: If present then this is used
!                    in place of the command line.
! Output:
!   error -- Logical: Set True if there is an error. False otherwise.
!-

subroutine tao_parse_command_args (error, cmd_words)

implicit none

character(*), optional :: cmd_words(:)
character(80) arg0, base, switch
character(24) :: r_name = 'tao_parse_command_args'

integer n_arg, i_arg, ix
logical error

! Get command line input

error = .false.

call tao_hook_parse_command_args()
if (.not. s%com%parse_cmd_args) return

if (present(cmd_words)) then
  n_arg = size(cmd_words)
  if (cmd_words(1) == '') return
else
  n_arg = cesr_iargc()
  if (n_arg == 0) return
endif

! loop over all arguments

i_arg = 0

do 

  if (i_arg == n_arg) exit
  call get_next_arg (arg0)

  call match_word (arg0,            ['-?                       ', &
        '-init                    ', '-noinit                  ', '-beam_all                ', '-beam0                   ', &
        '-noplot                  ', '-lat                     ', '-log_startup             ', '-beam                    ', &
        '-var                     ', '-data                    ', '-building_wall           ', '-plot                    ', &
        '-startup                 ', 'help                     ', '-help                    ', '?                        ', &
        '-geometry                ', '-rf_on                   ', '-debug                   ', '-disable_smooth_line_calc', &
        '-color_prompt            '], &
              ix, .true., matched_name=switch)

  select case (switch)

  case ('-init')
    call get_next_arg (s%com%init_tao_file)
    s%com%init_tao_file_arg_set = .true.
    if (s%com%init_tao_file == '') then
      call out_io (s_fatal$, r_name, 'NO TAO INIT FILE NAME ON COMMAND LINE.')
      call err_exit
    endif
    ix = SplitFileName(s%com%init_tao_file, s%com%init_tao_file_path, base)

  case ('-beam')
    call get_next_arg (s%com%beam_file)

  case ('-beam_all')
    call get_next_arg (s%com%beam_all_file)

  case ('-beam0')
    call get_next_arg (s%com%beam0_file)

  case ('-building_wall')
    call get_next_arg (s%com%building_wall_file)

  case ('-color_prompt')
    s%global%prompt_color = 'BLUE'

  case ('-data')
    call get_next_arg (s%com%data_file)

  case ('-disable_smooth_line_calc')
    s%global%disable_smooth_line_calc = .true.

  case ('-debug')
    s%global%debug_on = .true.

  case ('-geometry')
    call get_next_arg (s%com%plot_geometry)

  case ('help', '-help', '?', '-?')
    call tao_print_command_line_info
    stop

  case ('-lat')
    call get_next_arg (s%com%lat_file)

  case ('-log_startup')
    s%com%log_startup = .true.

  case ('-noinit')
    s%com%init_tao_file = ''

  case ('-noplot')
    s%com%noplot_arg_set = .true.

  case ('-rf_on')
    s%global%rf_on = .true.

  case ('-plot')
    call get_next_arg (s%com%plot_file)

  case ('-startup')
    call get_next_arg (s%com%startup_file)

  case ('-var')
    call get_next_arg (s%com%var_file)

  case default
    call out_io (s_error$, r_name, 'BAD COMMAND LINE ARGUMENT: ' // arg0)
    call tao_print_command_line_info
    error = .true.
    return
  end select

enddo

!-----------------------------
contains

subroutine get_next_arg(arg)

character(*) arg

!

if (i_arg == n_arg) then
  call out_io (s_error$, r_name, 'MISSING COMMAND LINE ARGUMENT FOR: ' // arg0)
  error = .true.
  return
endif

i_arg = i_arg + 1

if (present(cmd_words)) then
  arg = cmd_words(i_arg)
else
  call cesr_getarg(i_arg, arg)
endif

end subroutine get_next_arg

end subroutine tao_parse_command_args

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine tao_print_command_line_info
!
! Routine to print a list of the command line options.
!-

subroutine tao_print_command_line_info

implicit none

character(40), parameter :: r_name = 'tao_print_command_line_info'

!

call out_io (s_blank$, r_name, [ &
        'Syntax:                                                                           ', &
        '  <path-to-tao-exe-directory>/tao {OPTIONS}                                                                   ', &
        'Options are:                                                                      ', &
        '  -beam <beam_file>               # Beam init particle positions                  ', &
        '  -beam0 <beam0_file>             # Beam init params (beam size, etc.)            ', &
        '  -beam_all <all_beam_file>       # Beam info from previous tracking              ', &
        '  -building_wall <wall_file>      # Define the building tunnel wall               ', &
        '  -color_prompt                   # Set color of prompt string to blue            ', &
        '  -data <data_file>               # Define data for plotting and optimization     ', &
        '  -disable_smooth_line_calc       # Disable the smooth line calc used in plotting ', &
        '  -geometry <width>x<height>      # Plot window geometry                          ', &
        '  -init <tao_init_file>           # Tao init file                                 ', &
        '  -lat <bmad_lattice_file>        # Bmad lattice file                             ', &
        '  -lat xsif::<xsif_lattice_file>  # XSIF lattice file                             ', &
        '  -log_startup                    # Write startup debugging info                  ', &
        '  -noinit                         # Do not useTao init file                       ', &
        '  -noplot                         # Do not open a plotting window                 ', &
        '  -plot <plot_file>               # Define plot setup info                        ', &
        '  -rf_on                          # Keep RF on (Default is to turn off)           ', &
        '  -startup <starup_command_file>  # Commands to run after parsing Tao init file   ', &
        '  -var <var_file>                 # Define variables for plotting and optimization'])


end subroutine tao_print_command_line_info

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Function tao_pointer_to_universe (ix_uni) result (u)
!
! Routine to set a pointer to a universe.
! If ix_uni is -1 then u(s%com%default_universe) will be used.
!
! Also see:
!   tao_pick_universe
!
! Input:
!   ix_uni -- Integer: Index to the s%u(:) array
!
! Output:
!   u      -- Tao_universe_struct, pointer: Universe pointer.
!               u will be nullified if ix_uni is out of range.
!               If not present then print a message and exit program.
!-

function tao_pointer_to_universe (ix_uni) result(u)

implicit none

type (tao_universe_struct), pointer :: u
integer ix_uni, ix_u
character(28) :: r_name = 'tao_pointer_to_universe'

!

ix_u = tao_universe_number(ix_uni)

if (ix_u < lbound(s%u, 1) .or. ix_u > ubound(s%u, 1)) then
  call out_io (s_fatal$, r_name, 'UNIVERSE INDEX OUT OF RANGE: \I0\ ', ix_u)
  nullify (u)
  return
endif

u => s%u(ix_u)

end function tao_pointer_to_universe

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine tao_ele_to_ele_track (ix_universe, ix_branch, ix_ele, ix_ele_track)
!
! Subroutine to compute ix_ele_track:
!   = ix_ele                 if ix_ele <= lat%branch(ix)branch)%n_ele_track
!   = ix_slave_at_exit_end   if ix_ele is a super_lord  
!   = -1                     otherwise
!
! Input:
!   ix_universe -- Integer: Universe index.
!   ix_branch   -- Integer: Branch index.
!   ix_ele      -- Integer: Element index
!
! Output:
!   ix_ele_track -- Integer: Corresponding element in the tracking 
!                         part of the lattice.
!-

subroutine tao_ele_to_ele_track (ix_universe, ix_branch, ix_ele, ix_ele_track)

implicit none

type (lat_struct), pointer :: lat

integer ix_universe, ix_branch, ix_ele, ix_ele_track
integer i_uni, ix_c

!

i_uni = tao_universe_number(ix_universe)
lat => s%u(i_uni)%model%lat

if (ix_ele < 0) then
  ix_ele_track = -1

elseif (ix_ele <= lat%branch(ix_branch)%n_ele_track) then
  ix_ele_track = ix_ele

elseif (lat%ele(ix_ele)%lord_status == super_lord$) then
  ix_c = lat%ele(ix_ele)%ix2_slave
  ix_ele_track = lat%control(ix_c)%slave%ix_ele ! element at exit end.

else  ! overlays, multipass_lords, etc.
  ix_ele_track = -1
endif

end subroutine tao_ele_to_ele_track

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine tao_string_to_element_id (str, ix_class, ele_name, err, print_err)
!
! Routine to split a string in the form str = "xxx::yyy" into an element class
! and an element name. Example: 
!   str = "quad::q*".
! gives
!   ix_class = quadrupole$
!   ele_name = "Q*"
!
! If "::" is not found then ix_class is set to 0 (all classes).
! If str is of the form: "*::yyy" then ix_class is set to 0 (all classes).
! Class abbreviations and lower case names accepted 
! ele_name will be converted to upper case
!
! Input:
!   str       -- Character(*): Character string to parse.
!   print_err -- Logical, optional: If True then print an error message if 
!                 there is a problem. Default is True.
!
! Output:
!   ix_class  -- Integer: Element class. 0 => all classes. -1 => Not an element [EG: str = "var::..."].
!   ele_name  -- Character(*): Element name.
!   err       -- Set true if there is a problem translating the element class.
!-

subroutine tao_string_to_element_id (str, ix_class, ele_name, err, print_err)

implicit none

integer ix, ix_class

character(*) str, ele_name
character(40) :: r_name = 'tao_string_to_element_id'
character(20) class

logical, optional :: print_err
logical err

! 

if (index(str, 'dat::') /= 0) then
  call out_io (s_error$, r_name, 'NAME USES OLD "dat::" SYNTAX. PLEASE CHANGE TO "data::": ' // str)
  call err_exit
endif

!

err = .false.
ix_class = -1

if (str(1:6) == 'data::') return
if (str(1:5) == 'var::') return
if (str(1:5) == 'lat::') return
if (str(1:6) == 'wall::') return

ix = index(str, '::')

if (ix == 0) then
  ix_class = 0
  ele_name = str
  call str_upcase (ele_name, ele_name)
  return
endif

class = str(:ix-1)
ele_name = str(ix+2:)
call str_upcase (ele_name, ele_name)

if (class == '*') then
  ix_class = 0
  return
endif

ix_class = key_name_to_key_index (class, .true.)
if (ix_class < 1) then
  if (logic_option (.true., print_err)) call out_io (s_error$, r_name, 'BAD CLASS NAME: ' // class)
  err = .true.
endif

end subroutine tao_string_to_element_id

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine floor_to_screen_coords (graph, floor, screen)
!
! Routine to project a 3D floor coordinate onto a 2D projection plane
! along with projecting the orientation vector.
!
! Input:
!   graph  -- tao_graph_struct: Graph defining the projection plane.
!   floor  -- floor_position_struct: 3D coordinate.
!
! Output:
!   screen  -- floor_position_struct: Projected point
!     %r(3)   -- projected (x, y) = (%r(1), %r(2))
!     %theta  -- angle in (x, y) plane of projection of the orientation vector.
!-

subroutine floor_to_screen_coords (graph, floor, screen)

implicit none

type (tao_graph_struct) graph
type (floor_position_struct) floor, screen
real(rp) orient(3), theta, phi, x, y

! Get projection position

call floor_to_screen (graph, floor%r, screen%r(1), screen%r(2))

! screen%theta does not depend upon floor%psi

theta = floor%theta
phi = floor%phi
orient = [sin(theta) * cos(phi), sin(phi), cos(theta) * cos(phi)]  ! orientation vector

select case (graph%floor_plan_view(1:1))
case ('x')
  x = orient(1)
case ('y')
  x = orient(2)
case ('z')
  x = orient(3)
end select

select case (graph%floor_plan_view(2:2))
case ('x')
  y = orient(1)
case ('y')
  y = orient(2)
case ('z')
  y = orient(3)
end select

screen%theta = atan2(y, x) + twopi * graph%floor_plan_rotation

end subroutine floor_to_screen_coords

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine floor_to_screen (graph, floor, x_screen, y_screen)
!
! Routine to project a 3D floor coordinate onto a 2D projection plane.
!
! Input:
!   graph  -- tao_graph_struct: Graph defining the projection plane.
!   floor  -- floor_position_struct: 3D coordinate.
!
! Output:
!   x_screen -- real(rp): x-coordinate of projected point.
!   y_screen -- real(rp): y-coordinate of projected point.
!-

subroutine floor_to_screen (graph, r_floor, x_screen, y_screen)

implicit none

type (tao_graph_struct) graph

real(rp) r_floor(3), x_screen, y_screen
real(rp) x, y
real(rp), save :: t, old_t = 0
real(rp), save :: cc, ss

! 

select case (graph%floor_plan_view(1:1))
case ('x')
  x = r_floor(1)
case ('y')
  x = r_floor(2)
case ('z')
  x = r_floor(3)
end select

select case (graph%floor_plan_view(2:2))
case ('x')
  y = r_floor(1)
case ('y')
  y = r_floor(2)
case ('z')
  y = r_floor(3)
end select

t = graph%floor_plan_rotation
if (t == 0) then
  x_screen = x
  y_screen = y
else
  if (t /= old_t) then
    cc = cos(twopi * t)
    ss = sin(twopi * t)
    old_t = t
  endif
  x_screen =  x * cc - y * ss
  y_screen =  x * ss + y * cc 
endif

end subroutine floor_to_screen

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_split_component (comp_str, comp, err)
!
! Routine to split a component string.
!
! Input:
!   comp_str -- Character(*): Components. EG: 'meas - design'
!
! Output:
!   comp(:) -- Tao_data_var_component_struct, allocatable: Array of individual components.
!   err     -- Logical: Set True if there is an error, False otherwise.
!-

subroutine tao_split_component (comp_str, comp, err)

type (tao_data_var_component_struct), allocatable :: comp(:)

integer i, n, ix, ix1, ix2

character(*) comp_str
character(60) str
character(40) :: r_name = 'tao_split_component'

logical err

! Count number of components.

err = .true.
call string_trim (comp_str, str, ix)
if (.not. (str(1:1) == '+' .or. str(1:1) == '-')) str = '+' // trim(str)

n = 0
do i = 1, len_trim(str)
  if (str(i:i) == '+' .or. str(i:i) == '-') n = n + 1
enddo

! Allocate space and transfer info.

if (allocated(comp)) then
  if (size(comp) /= n) deallocate (comp)
endif
if (.not. allocated(comp)) allocate (comp(n))

n = 0
do n = 1, size(comp)
  if (str(1:1) == '+') then
    comp(n)%sign = 1
  elseif (str(1:1) == '-') then
    comp(n)%sign = -1
  else
    call out_io (s_error$, r_name, 'BAD COMPONENT LIST: ' // comp_str)
    return
  endif

  call string_trim(str(2:), str, ix)

  if (n == size(comp)) then
    comp(n)%name = str
    exit
  endif

  ix1 = index(str, "+")
  ix2 = index(str, "-")
  if (ix1 /= 0 .and. ix2 /= 0) then
    ix = min(ix1, ix2) - 1
  else
    ix = max(ix1, ix2) - 1
  endif

  comp(n)%name = str(1:ix)
  call string_trim(str(ix+1:), str, ix)

enddo

err = .false.

end subroutine tao_split_component

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_turn_on_chrom_or_rad_int_calcs_if_needed_for_plotting ()
! 
! Routine to set u%dynch_rad_int_clac = T if needed for a plot.
!-

subroutine tao_turn_on_chrom_or_rad_int_calcs_if_needed_for_plotting ()

implicit none

type (tao_universe_struct), pointer :: u
type (tao_graph_struct), pointer :: graph
type (tao_curve_struct), pointer :: curve

integer i, j, k

character(*), parameter :: r_name = 'tao_turn_on_chrom_or_rad_int_calcs_if_needed_for_plotting'

! Go through all the plots and find which universes need chrom_or_rad_int_calcs.
! u%picked_uni = True => need calc.

s%u(:)%picked_uni  = s%u(:)%calc%rad_int_for_plotting
s%u(:)%picked2_uni = s%u(:)%calc%chrom_for_plotting

s%u(:)%calc%rad_int_for_plotting = .false.
s%u(:)%calc%chrom_for_plotting   = .false.

do i = 1, size(s%plot_page%region)
  if (.not. s%plot_page%region(i)%visible) cycle
  if (.not. allocated(s%plot_page%region(i)%plot%graph)) cycle

  do j = 1, size(s%plot_page%region(i)%plot%graph)
    graph => s%plot_page%region(i)%plot%graph(j)
    if (.not. allocated(graph%curve)) cycle

    do k = 1, size(graph%curve)
      curve => graph%curve(k)
      u => tao_pointer_to_universe(curve%ix_universe)

      if (.not. u%picked_uni .and. tao_rad_int_calc_needed(curve%data_type, curve%data_source)) then
        if (curve%ix_branch /= 0) then
          call out_io (s_fatal$, r_name, 'PLOTTING THIS: ' // curve%data_type, 'ON A BRANCH NOT YET IMPLEMENTED!')
          call err_exit
        endif
        u%calc%rad_int_for_plotting = .true.
        u%calc%lattice = .true.
        u%picked_uni = .true.
      endif

      if (.not. u%picked2_uni .and. tao_chrom_calc_needed(curve%data_type, curve%data_source)) then
        if (curve%ix_branch /= 0) then
          call out_io (s_fatal$, r_name, 'PLOTTING THIS: ' // curve%data_type, 'ON A BRANCH NOT YET IMPLEMENTED!')
          call err_exit
        endif
        u%calc%chrom_for_plotting = .true.
        u%calc%lattice = .true.
        u%picked2_uni = .true.
      endif

    enddo     ! graph%curve
  enddo     ! plot%graph
enddo     ! plotting%region

end subroutine tao_turn_on_chrom_or_rad_int_calcs_if_needed_for_plotting

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function tao_rad_int_calc_needed (data_type, data_source) result (do_rad_int)
! 
! Routine decide if a datum or plot curve needs the radiation integrals 
! to be evaluated.
!-

function tao_rad_int_calc_needed (data_type, data_source) result (do_rad_int)

implicit none

character(*) data_type, data_source
logical do_rad_int

!

do_rad_int = .false.

if (data_source /= 'lat') return

if (data_type  == 'sigma.pz') do_rad_int = .true. 
if (data_type(1:5)  == 'emit.') do_rad_int = .true. 
if (data_type(1:10) == 'norm_emit.') do_rad_int = .true. 
if (data_type(1:7)  == 'rad_int') do_rad_int = .true.
if (data_type(1:16)  == 'apparent_rad_int') do_rad_int = .true.

end function tao_rad_int_calc_needed

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function tao_chrom_calc_needed (data_type, data_source) result (do_chrom)
! 
! Routine decide if a datum or plot curve needs the chromaticity calculation.
!-

function tao_chrom_calc_needed (data_type, data_source) result (do_chrom)

implicit none

character(*) data_type, data_source
logical do_chrom

!

do_chrom = .false.

if (data_source /= 'lat') return
if (data_type(1:6)  == 'chrom.') do_chrom = .true. 

end function tao_chrom_calc_needed

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function tao_subin_uni_number (name_in, ix_uni, name_out) result (ok)
!
! Routine to mangle a file name based upon the universe number ix_uni.
! If there is a "#" in the name and the number of universes is more than one,
! the universe number will be substituted.
! If the number of universes is 1, the "#" will be stripped from the name.
!
! Input:
!   name_in -- Character(*): Input name with "#" character
!   ix_uni  -- Integer: Universe index.
!
! Output:
!   name_out -- Character(*): Output name.
!   ok       -- Logical: False if multiple universes and no "#" in name_in.
!                 True otherwise.
!-

function tao_subin_uni_number (name_in, ix_uni, name_out) result (ok)

character(*) name_in, name_out
character(28) :: r_name = 'tao_subin_uni_number'
integer ix, ix_uni
logical ok

!

ok = .true.
name_out = name_in

ix = index(name_out, '#')
if (size(s%u) > 1 .and. ix == 0) then
  call out_io (s_info$, r_name, 'FILE NAME DOES NOT HAVE A "#" CHARACTER!', &
    ' YOU NEED THIS TO GENERATE A UNIQUE FILE NAME FOR EACH UNIVERSE!')
  ok = .false.
endif

if (ix /= 0) write (name_out, '(a, i0, a)') &
                        name_out(1:ix-1), ix_uni, trim(name_out(ix+1:))

end function tao_subin_uni_number

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function tao_lat_emit_calc (plane, emit_type, ele, modes) result (emit)
!
! Routine to calculate the emittance.
! The "projected" emittance is:
!   emit = sqrt(sigma_xx * sigma_pp - sigma_xp^2)
! Where the sigmas are calculated including coupling effects but assuming that sigma_pz = 0.
! The "apparent" emittance is:
!   emit = sqrt(sigma_xx - sigma_xp^2 / sigma_pp) / beta
!
! Input:
!   plane     -- Integer: x_plane$ or y_plane$.
!   emit_type -- Integer: Either projected_emit$ or apparent_emit$
!   ele       -- ele_struct: Element holding the Twiss and coupling parameters.
!   modes     -- normal_modes_struct: Structure holding the emittances
!
! Output:
!   emit -- Real(rp): emittance.
!-

function tao_lat_emit_calc (plane, emit_type, ele, modes) result (emit)

implicit none

type (ele_struct) ele
type (normal_modes_struct) modes

real(rp) s_mat(4,4), v_mat(4,4), emit
real(rp), save :: a_mat(4,4) = 0
integer plane, emit_type

!

a_mat(1,1) =  modes%a%emittance * ele%a%beta
a_mat(1,2) = -modes%a%emittance * ele%a%alpha
a_mat(2,2) =  modes%a%emittance * ele%a%gamma
a_mat(2,1) = a_mat(1,2)

a_mat(3,3) =  modes%b%emittance * ele%b%beta
a_mat(3,4) = -modes%b%emittance * ele%b%alpha
a_mat(4,4) =  modes%b%emittance * ele%b%gamma
a_mat(4,3) = a_mat(3,4)

call make_v_mats (ele, v_mat)
s_mat = matmul(matmul(v_mat, a_mat), transpose(v_mat))

if (plane == x_plane$) then
  if (emit_type == projected_emit$) then
    emit = sqrt(s_mat(1,1) * s_mat(2,2) - s_mat(1,2)**2)
  elseif (emit_type == apparent_emit$) then
    emit = s_mat(1,1) / ele%a%beta
  else
    call err_exit    
  endif

elseif (plane == y_plane$) then
  if (emit_type == projected_emit$) then
    emit = sqrt(s_mat(3,3) * s_mat(4,4) - s_mat(3,4)**2)
  elseif (emit_type == apparent_emit$) then
    emit = s_mat(3,3) / ele%b%beta
  else
    call err_exit    
  endif

else
  call err_exit
endif

end function tao_lat_emit_calc

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function tao_beam_emit_calc (plane, emit_type, ele, bunch_params) result (emit)
!
! Routine to calculate the emittance from beam parameters.
! The "projected" emittance is:
!   emit = sqrt(sigma_xx * sigma_pp - sigma_xp^2)
! Where the sigmas are calculated including coupling effects but assuming that sigma_pz = 0.
! The "apparent" emittance is:
!   emit = sqrt(sigma_xx - sigma_xp^2 / sigma_pp) / beta
!
! Input:
!   plane        -- Integer: x_plane$ or y_plane$.
!   emit_type    -- Integer: Either projected_emit$ or apparent_emit$
!   ele          -- ele_struct: Element.
!   bunch_params -- bunch_params_struct: Bunch sigma matrix
!
! Output:
!   emit -- Real(rp): emittance.
!-

function tao_beam_emit_calc (plane, emit_type, ele, bunch_params) result (emit)

use beam_def_struct

implicit none

type (ele_struct) ele
type (bunch_params_struct) bunch_params
real(rp) emit
integer plane, emit_type

!

if (plane == x_plane$) then
  if (emit_type == projected_emit$) then
    emit = bunch_params%x%emit
  elseif (emit_type == apparent_emit$) then
    if (bunch_params%sigma(6,6) == 0) then
      emit = bunch_params%sigma(1,1) / ele%a%beta
    else
      emit = (bunch_params%sigma(1,1) - bunch_params%sigma(1,6)**2 / bunch_params%sigma(6,6)) / ele%a%beta
    endif
  else
    call err_exit    
  endif

elseif (plane == y_plane$) then
  if (emit_type == projected_emit$) then
    emit = bunch_params%y%emit
  elseif (emit_type == apparent_emit$) then
    if (bunch_params%sigma(6,6) == 0) then
      emit = bunch_params%sigma(3,3) / ele%b%beta
    else
      emit = (bunch_params%sigma(3,3) - bunch_params%sigma(3,6)**2 / bunch_params%sigma(6,6)) / ele%b%beta
    endif
  else
    call err_exit    
  endif

else
  call err_exit
endif

end function tao_beam_emit_calc

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function tao_pointer_to_ele_shape (ix_uni, ele, ele_shape, dat_var_name, dat_var_value, ix_shape) result (e_shape)
!
! Routine to return the shape associated with a lattice element
!
! Input:
!   ix_uni        -- integer: Universe index.
!   ele           -- ele_struct: Lattice element.
!   ele_shape(:)  -- tao_ele_shape_struct: Array of shapes to search.
!
! Output:
!   e_shape       -- tao_ele_shape_struct, pointer: Associated shape. 
!                       Nullified if there is no associated shape.
!   dat_var_name  -- character(*), optional: Name of datum or variable associated with e_shape. 
!                       Will be set to "" if there is no associated datum or variable.
!   dat_var_value -- real(rp), optional: Value of datum or variable associated with e_shape.
!                       Will be set to zero if there is no associated datum or variable.
!   ix_shape      -- integer, optional: Index of associated shape in ele_shape(:).
!-

function tao_pointer_to_ele_shape (ix_uni, ele, ele_shape, dat_var_name, dat_var_value, ix_shape) result (e_shape)

implicit none

type (ele_struct) ele

type (tao_ele_shape_struct), target :: ele_shape(:)
type (tao_ele_shape_struct), pointer :: e_shape, es
type (tao_data_array_struct), allocatable, target :: d_array(:)
type (tao_logical_array_struct), allocatable :: logic_array(:)
type (tao_data_struct), pointer :: datum
type (tao_var_array_struct), allocatable, target :: v_array(:)
type (tao_var_struct), pointer :: var
type (tao_real_pointer_struct), allocatable :: re_array(:)

real(rp), optional :: dat_var_value

integer, optional :: ix_shape
integer ix_uni
integer j, j2, k, n_ele_track

character(*), optional :: dat_var_name
character(*), parameter :: r_name = 'tao_pointer_to_ele_shape'

logical err

!

nullify(e_shape)
if (present(dat_var_name)) dat_var_name = ''
if (present(dat_var_value)) dat_var_value = 0

if (ele%lord_status == group_lord$) return
if (ele%lord_status == overlay_lord$) return
if (ele%slave_status == super_slave$) return

do k = 1, size(ele_shape)
  es => ele_shape(k)
  if (present(ix_shape)) ix_shape = k
  if (.not. es%draw) cycle

  ! Data

  if (index(es%ele_id, 'dat::') /= 0) then
    call out_io (s_error$, r_name, 'SHAPE USES OLD "dat::" SYNTAX. PLEASE CHANGE TO "data::": ' // es%ele_id)
    call err_exit
  endif

  if (es%ele_id(1:6) == 'data::') then
    call tao_find_data (err, es%ele_id, d_array = d_array, log_array = logic_array, re_array = re_array)
    if (err) cycle
    do j = 1, size(d_array)
      datum => d_array(j)%d
      if (datum%d1%d2%ix_uni /= ix_uni) cycle
      if (size(logic_array) /= 0) then
        if (.not. logic_array(j)%l) cycle
      endif
      if (ele%ix_branch /= datum%ix_branch .or. ele%ix_ele /=  datum%ix_ele) cycle
      e_shape => es
      if (present(dat_var_name)) dat_var_name = tao_datum_name(datum)
      if (present(dat_var_value)) then
        if (size(re_array) > 0) then
          dat_var_value = re_array(j)%r
        else
          dat_var_value = datum%model_value
        endif
      endif
      return
    enddo
    cycle
  endif

  ! Variables

  if (es%ele_id(1:5) == 'var::') then
    call tao_find_var (err, es%ele_id, v_array = v_array, log_array = logic_array, re_array = re_array)
    if (err) cycle

    do j = 1, size(v_array)
      var => v_array(j)%v
      if (size(logic_array) /= 0) then
        if (.not. logic_array(j)%l) cycle
      endif
      do j2 = 1, size(var%this)
        if (var%this(j2)%ix_uni /= ix_uni) cycle
        if (var%ele_name == 'BEAM_START') then
          if (ele%ix_branch /= 0 .or. ele%ix_ele /= 0) cycle
        else
          if (ele%ix_branch /= var%this(j2)%ix_branch .or. ele%ix_ele /=  var%this(j2)%ix_ele) cycle
        endif
        e_shape => es
        if (present(dat_var_name)) dat_var_name = tao_var1_name(var)
        if (present(dat_var_value)) then
          if (size(re_array) > 0) then
            dat_var_value = re_array(j)%r
          else
            dat_var_value = var%this(j2)%model_value
          endif
        endif
        return
      enddo
    enddo
    cycle
  endif

  ! All else

  if (es%ele_id == '') cycle
  if (es%ix_ele_key == -1) cycle
  if (es%ix_ele_key /= 0 .and. es%ix_ele_key /= ele%key) cycle
  if (.not. match_wild(ele%name, es%name_ele)) cycle

  e_shape => es
  return
enddo

end function tao_pointer_to_ele_shape 

!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
! This searches the lattice for the specified element and flags eles(:)

subroutine tao_init_find_elements (u, search_string, eles, attribute)

implicit none

type (tao_universe_struct), target :: u
type (ele_struct), pointer :: ele
type (ele_pointer_struct), allocatable :: eles(:)

character(*) search_string
character(*), optional :: attribute
character(80) string
character(40) ele_name, key_name_in
character(20) :: r_name = 'tao_init_find_elements'

integer key, found_key, ix_attrib, t
integer i, k, ix, ii, j, ix0, ix1, ix2, n_ele

logical no_slaves, no_lords, err, warn_given

! Sort switches

no_slaves = .false.
no_lords = .false.

call string_trim(search_string, string, ix)

do

  if (string(1:1) /= '-') exit

  select case (string(1:ix))
  case ('-no_lords') 
    no_lords = .true.
  case ('-no_slaves') 
    no_slaves = .true.
  case default
    call out_io (s_abort$, r_name, "BAD SEARCH SWITCH: " // search_string)
    call err_exit
  end select

  call string_trim (string(ix+1:), string, ix)

enddo

! Find elements

call tao_locate_elements (string, u%ix_uni, eles, err, print_err = .false.)

warn_given = .false.
n_ele = 0
do j = 1, size(eles)
  ele => eles(j)%ele
  t = ele%slave_status
  if ((t == multipass_slave$ .or. t == super_slave$) .and. no_slaves) cycle
  t = ele%lord_status 
  if ((t == girder_lord$ .or. t == overlay_lord$ .or. t == group_lord$ .or. &
       t == super_lord$ .or. t == multipass_lord$) .and. no_lords) cycle
  ! If attribute is not free then don't count it
  if (present(attribute)) then
    ix_attrib = attribute_index(ele, attribute)
    if (ix_attrib < 1) then
      call out_io (s_error$, r_name, &
                      'BAD ATTRIBUTE: ' // attribute, &
                      'FOR VARIABLE SEARCH: ' // search_string, &
                      'FOR ELEMENT: ' // ele%name)
      return
    endif
    if (.not. attribute_free (eles(j)%ele, attribute, .false.)) then
      if (.not. warn_given) then
        call out_io (s_info$, r_name, &
                  'Non-free attribute ' // attribute, &
                  'For variable search: ' // search_string, &
                  'For element: ' // ele%name)
        warn_given = .true.
      endif
      cycle
    endif
  endif
  ! Passes test so add it to the list
  n_ele = n_ele + 1
  eles(n_ele)%ele => eles(j)%ele
enddo

call re_allocate_eles (eles, n_ele, .true., .true.)

end subroutine tao_init_find_elements

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! Put the variables marked by key_bound in the key table.

subroutine tao_setup_key_table ()

implicit none

integer i, j

! Key table setup

call re_allocate (s%key, count(s%var%key_bound))
s%key = -1

j = 0
do i = 1, s%n_var_used
  if (s%var(i)%key_bound .and. s%var(i)%exists) then
    j = j + 1
    s%key(j) = i
    s%var(i)%key_val0 = s%var(i)%model_value
    s%var(i)%ix_key_table = j
  else
    s%var(i)%key_bound = .false.
    s%var(i)%ix_key_table = -1
  endif
enddo

end subroutine tao_setup_key_table

end module tao_utils
