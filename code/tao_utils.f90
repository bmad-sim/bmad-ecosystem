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
if (ios /= 0) then
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
! Subroutine tao_pointer_to_var_in_lattice (var, lattice, orbit, ptr, err)
! 
! Routine to set a pointer to the appropriate variable in a lattice
!
! Input:
!   var       -- Tao_var_struct: Structure has the info of where to point.
!   lattice   -- Ring_struct: lattice to set.
!   orbit(0:) -- Coord_struct: Orbit.
!
! Output:
!   ptr   -- Real(rp), pointer: Pointer set to appropriate variable.
!   err   -- Logical: Set True if there is an error. False otherwise.
!-

subroutine tao_pointer_to_var_in_lattice (var, lattice, orbit, ptr, err)

implicit none

type (tao_var_struct) var
type (ring_struct) lattice
type (coord_struct) orbit(0:)

real(rp), pointer :: ptr
integer ix
logical err
character(30) :: r_name = 'tao_pointer_to_var_in_lattice'

! locate element

call element_locator (var%ele_name, lattice, ix)
if (ix < 0) then
  call out_io (s_error$, r_name, 'ELEMENT NAME NOT FOUND: ' // var%ele_name)
  err = .true.
  return
endif

! locate attribute

call pointer_to_attribute (lattice%ele_(ix), var%attrib_name, .true., ptr, ix, err)

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_top10_print (s)
!
! Routine to print out the top10 contributors to the merit function.
!
! Input:
!   s - Tao_super_universe_struct
!-

subroutine tao_top10_print (s)

implicit none

type (tao_super_universe_struct) s
type (tao_top10_struct) top_merit(10)
type (tao_top10_struct) top_dmerit(10)
type (tao_top10_struct) top_delta(10)

real(rp) delta, a_max
integer i, j, n, nl, nu

character(16) name
character(80) fmt, lines(20)
character(20) :: r_name = 'tao_top10_print'

! Top merit

top_merit(:)%valid  = .false.
top_dmerit(:)%valid = .false.
top_delta(:)%valid  = .false.

nu = size(s%u)
do i = 1, nu
  do j = 1, size(s%u(i)%data)
    if (.not. s%u(i)%data(j)%useit_opt) cycle
    name = trim(s%u(i)%data(j)%d1%d2%class) // ':' // s%u(i)%data(j)%d1%sub_class
    if (nu > 1) write (name, '(2a, i1)') trim(name), '/', j
    call tao_to_top10 (top_merit, s%u(i)%data(j)%merit, name, &
                                                s%u(i)%data(j)%ix_d, 'max')
  enddo
enddo


if (s%global%parallel_vars) then
  do j = 1, size(s%u(1)%var)
    if (.not. s%u(1)%var(j)%useit_opt) cycle
    name = s%u(1)%var(j)%v1%class
    call tao_to_top10 (top_merit, s%u(1)%var(j)%merit, name, s%u(1)%var(j)%ix_v, 'max')
    call tao_to_top10 (top_dmerit, s%u(1)%var(j)%dmerit_dvar, name, &
                                                  s%u(1)%var(j)%ix_v, 'max')
    delta = s%u(1)%var(j)%model_value - s%u(1)%var(j)%design_value
    call tao_to_top10 (top_dmerit, delta, name, s%u(1)%var(j)%ix_v, 'max')
  enddo
else
  call out_io (s_abort$, r_name, 'SERIES VARIABLES NOT YET IMPLEMENTED')
  call err_exit
endif

! write results


a_max = max(1.1, maxval(abs(top_delta(:)%value)))
n = max(0, 6 - int(log10(a_max)))

write (fmt, '(a, i1, a)') &
    '((1x, a10, i3, f10.1, 2x), (a8, i4, 1pe12.3, 2x), a8, i4, 0pf11.', n, ')'


nl = 0
lines(nl+1) = ' '
lines(nl+2) = '       Top10 merit      |    Top10 derivative     |     Top10 delta'
lines(nl+3) = '  Name     ix     Value | Name     ix  Derivative |  Name    ix      delta'
nl = nl + 3

do i = 1, 10
  nl = nl + 1
  write (lines(nl), fmt) &
      top_merit(i)%name,  top_merit(i)%index,  top_merit(i)%value, &
      top_dmerit(i)%name, top_dmerit(i)%index, top_dmerit(i)%value,  &
      top_delta(i)%name,  top_delta(i)%index,  top_delta(i)%value
enddo

nl = nl + 1
write (lines(nl), *) 'Merit:  ', tao_merit(s)

call out_io (s_blank$, r_name, lines(1:nl))

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_to_top10 (top10, value, name, c_index, order)
!
! Routine to order the largest contributors to the merit function in
! a list. Call this routine for each contributor.
!
! Note: Before first calling this routine set:
!   top10(:)%valid = .false.
!
! Input:
!   value   -- Real(rp): value of the contributor.
!   name    -- Character(16): Name of the contributor..
!   c_index -- Integer: Index of the contributor.
!   order   -- Character(16): Ordering of the list. Possibilities are:
!                 'max' -- #1 has the maximum value.
!                 'min' -- #1 has the minimum value.
!                 'abs' -- #1 has the maximum aplitude.
!
! Output:
!   top10(:) -- Tao_top10_struct: List of top contributors.
!                 Note that the list is not limited to 10 entries.
!-

subroutine tao_to_top10 (top10, value, name, c_index, order)

implicit none

type (tao_top10_struct) top10(:)

integer c_index, ix, n
real(rp) value

character(*) name, order
character(20) :: r_name = 'tao_to_top10'

! Find where in list the current contributor is.

n = size(top10)
do ix = n, 1, -1
  if (.not. top10(ix)%valid) cycle
  select case (order)
  case ('max')
    if (value < top10(ix)%value) exit
  case ('min')
    if (value < top10(ix)%value) exit
  case ('abs')  
    if (abs(value) > abs(top10(ix)%value)) exit
  case default
    call out_io (s_abort$, r_name, 'BAD "ORDER" ARGUMENT: ' // order)
  end select
enddo

ix = ix + 1          ! place to put current contributor.
if (ix > n) return   ! not big enough to be in list.

! Move the people below the current contributor down to make room and
! then put the contributor in.

top10(ix+1:n) = top10(ix:n-1) 
top10(ix) = tao_top10_struct(name, value, c_index, .true.)

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
!     %class        -- Match name if by_who = "BY_TYPE"
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

type (tao_plot_struct) :: plots(:)
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
    if (plot_name == plots(i)%class) exit
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
! Subroutine tao_dModel_dVar_calc (s)
!
! Subroutine to calculate the dModel_dVar derivative matrix.
!
! Input:
!   s       -- Super_universe_struct.
!
! Output:
!   s       -- Super_universe_struct.
!    %u(:)%dModel_dVar(:,:)  -- Derivative matrix
!-

subroutine tao_dModel_dvar_calc (s)

implicit none

type (tao_super_universe_struct), target :: s
type (tao_universe_struct), pointer :: u

integer i, j, m, n, k
integer n_data, n_var
character(20) :: r_name = 'tao_dmodel_dvar_calc'
logical reinit

! make sure size of matrix is correct.

reinit = .false.

do i = 1, size(s%u)

  u => s%u(i)
  n_data = count (u%data%useit_opt)
  n_var = count (u%var%useit_opt)

  if (.not. associated(u%dModel_dVar)) then
    allocate (u%dModel_dVar(n_data, n_var))
    reinit = .true.
    cycle
  endif

  if (size(u%dModel_dVar, 1) /= n_data .or. size(u%dModel_dVar, 2) /= n_var) then
    deallocate (u%dModel_dVar)
    allocate (u%dModel_dVar(n_data, n_var))
    reinit = .true.
    cycle
  endif

enddo

if (.not. reinit) return
call out_io (s_info$, r_name, 'Remaking dModel_dVar derivative matrix') 

! Calculate matrix

do i = 1, size(s%u)

  u => s%u(i)
  n_data = count (u%data%useit_opt)
  n_var = count (u%var%useit_opt)

  call tao_lattice_calc (s)
  m = 0
  do j = 1, size(u%data)
    if (u%data(j)%useit_opt) cycle
    m = m + 1
    u%data(j)%ix_dModel = m
    u%data(j)%old_value = u%data(j)%delta
  enddo



  m = 0
  do j = 1, size(u%var)
    if (.not. u%var(j)%useit_opt) cycle
    m = m + 1
    u%var(j)%ix_dVar = m
    if (u%var(j)%step == 0) then
      call out_io (s_error$, r_name, 'VARIABLE STEP SIZE IS ZERO FOR: ' // u%var(j)%name)
      call err_exit
    endif
    call tao_lattice_calc (s)
    u%var(j)%model_value = u%var(j)%model_value + u%var(j)%step
    call tao_lattice_calc (s)
    n = 0
    do k = 1, size(u%data)
      if (u%data(k)%useit_opt) cycle
      n = n + 1
      u%dModel_dVar(m,n) = (u%data(k)%delta - u%data(k)%old_value) / u%var(j)%step
    enddo
  enddo

enddo

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
!   plot     -- Tao_super_universe_struct:
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

data%useit_plot = data%exists .and. data%good_user
if (any(plot%who%name == 'data')) data%useit_plot = data%useit_plot .and. data%good_data
if (any(plot%who%name == 'ref'))  data%useit_plot = data%useit_plot .and. data%good_ref

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_find_data (err, u, data_class, d2_ptr, d1_ptr, &
!                                                     		data_number, d_ptr)
!
! Routine to set data pointers to the correct data.
!
! Input:
!   u		          -- tao_universe_struct
!   data_class    -- character(*): the data class type. Eg: "orbit:x"
!   data_number   -- character(*), optional: the data point index.
!                     If data_number = 'null' then d_ptr will be nullified.
!
! Output:
!   err 		-- logical: err condition
!   d2_ptr	-- tao_d2_data_struct, optional, pointer: to the d2 data
!   d1_ptr 	-- tao_d1_data_struct, optional: pointer to the d1 data
!   d_ptr   -- tao_data_struct, optional: pointer to the data point
!-

subroutine tao_find_data (err, u, data_class, d2_ptr, d1_ptr, &
                                                     			data_number, d_ptr)

implicit none

type (tao_universe_struct), target           :: u
type (tao_d2_data_struct), pointer, optional :: d2_ptr
type (tao_d1_data_struct), pointer, optional :: d1_ptr
type (tao_data_struct), pointer, optional    :: d_ptr
type (tao_d2_data_struct), pointer :: d2_pointer
type (tao_d1_data_struct), pointer :: d1_pointer

character(*)                                :: data_class
character(*), optional                      :: data_number
character(20) :: r_name = 'TAO_FIND_DATA'
character(16) name, class, sub_class

integer :: data_num, ios
integer i, ix, ix_plane

logical err

! init

err = .false.

if (present(d2_ptr)) nullify(d2_ptr)
if (present(d1_ptr)) nullify(d1_ptr)
if (present(d_ptr)) nullify(d_ptr)

ix = index(data_class, ':')
if (ix == 0) then
  class = data_class
  sub_class = ' '
else
  class = data_class(1:ix-1)
  sub_class = data_class(ix+1:)
endif

! Point to the correct d2 data type 

  call string_trim (class, name, ix) ! Strip off all whitespace

  do i = 1, size(u%d2_data)
    if (name == u%d2_data(i)%class) then
      d2_pointer => u%d2_data(i) 
      if (present(d2_ptr)) d2_ptr => d2_pointer
      exit
    endif
    if (i == size(u%d2_data)) then
      call out_io (s_error$, r_name, "Couldn't find data class")
      err = .true.
      return
    endif
  enddo

  if (sub_class == ' ') return

! strip off all whitespace

  call string_trim (sub_class, name, ix)
  
! point to the correct d1 data type

  do i = 1, size(d2_pointer%d1)
    if (sub_class == d2_pointer%d1(i)%sub_class) then
      d1_pointer => d2_pointer%d1(i)
      if (present(d1_ptr)) d1_ptr => d1_pointer
      exit
    endif
    if (i .eq. size(d2_pointer%d1)) then
      call out_io (s_error$, r_name, &
                      "COULDN'T FIND DATA SUB_CLASS: " // sub_class)
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
    call out_io (s_error$, r_name, "DATA_NUMBER OUT OF RANGE: " // data_number)
    err = .true.
    return	
  endif

  if (present(d_ptr)) d_ptr => d1_pointer%d(data_num)

end subroutine tao_find_data

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine: tao_find_var
!
! find a v1 variable type, and variable data then point to it
!
!Input:
! u		-- tao_universe_struct
! var_class     -- character(*): the variable class type
! var_number    -- integer: (optional) the variable data point.
!                     If var_number = 'null' then d_ptr will be nullified.
!
!Output:
! err		-- logical: err condition
! v1_ptr:	-- tao_v1_var_struct: pointer to the v1 variable
! v_ptr:        -- tao_var_struct: (optional) pointer to the variable data point
!-

subroutine tao_find_var (err, u, var_class, v1_ptr, var_number, v_ptr)

implicit none

logical err
type (tao_universe_struct), target           :: u
type (tao_v1_var_struct), pointer, optional  :: v1_ptr
type (tao_v1_var_struct), pointer :: v1
type (tao_var_struct), pointer, optional     :: v_ptr

character(*)                                 :: var_class
character(*), optional                       :: var_number
character(20) :: r_name = 'tao_find_var'
character(32) v_class

integer i, ix, n_var, ios

err = .false.

! Strip off all whitespace

if (present(v1_ptr)) nullify (v1_ptr)
if (present(v_ptr)) nullify (v_ptr)

call string_trim(var_class, v_class, ix)
if (ix == 0) then
  err = .true.
  call out_io (s_error$, r_name, 'VARIABLE CLASS NAME IS BLANK')
  return
endif

! Point to the correct v1 data type 

  do i = 1, size(u%v1_var)
    if (v_class == u%v1_var(i)%class) then
      v1 => u%v1_var(i) 
      if (present(v1_ptr)) v1_ptr => u%v1_var(i) 
      exit
    endif
    if (i .eq. size(u%v1_var)) then
      call out_io (s_error$, r_name, &
                        "COULD NOT FIND VARIABLE CLASS: " // v_class)
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
! Subroutine tao_var_target_calc (s)
! 
! Subroutine to calculate the variable target values (the values that they need
! to be set to to do a correction of the orbit, phase, etc.
!
! Input:
!   s       -- Tao_super_universe_struct:
!
! Output:
!   s       -- Tao_super_universe_struct:
!-

subroutine tao_var_target_calc (s)

implicit none

type (tao_super_universe_struct), target :: s
type (tao_var_struct), pointer :: var

integer i, j

!

do i = 1, size(s%u)
  do j = 1, size(s%u(i)%var)

    var => s%u(i)%var(j)
    var%target_value = var%data_value + (var%design_value - var%model_value)

  enddo
enddo

end subroutine

end module tao_utils
