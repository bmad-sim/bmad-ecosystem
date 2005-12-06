module tao_set_mod

use tao_mod
use quick_plot
use tao_lattice_calc_mod

contains

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! Subroutine tao_set_lattice_cmd (set_lattice, to_lattice)
!
! Sets a lattice equal to another. This will also update the data structs
! If the 
!
! Input:
!   set_lattice -- Character(*): Maybe: 'model', 'design', or 'base' with 
!                     optional ';n' at end to indicate the universe
!   to_lattice  -- Character(*): Maybe: 'model', 'design', or 'base' 
!
!  Output:
!    s%u(n) -- ring_struct: changes specified lattice in specified universe 
!-

subroutine tao_set_lattice_cmd (set_lattice, to_lattice)

implicit none

character(*) set_lattice, to_lattice
character(16) set_lat_name
character(20) :: r_name = 'tao_set_lattice_cmd'

integer i

logical, automatic :: this_u(size(s%u))
logical err

call tao_pick_universe (set_lattice, set_lat_name, this_u, err)
if (err) return

do i = 1, size(s%u)
  if (.not. this_u(i)) cycle
  call set_lat (s%u(i))
  if (err) return
enddo

!-------------------------------------------
contains

subroutine set_lat (u)

implicit none

type (tao_universe_struct), target :: u
type (ring_struct), pointer :: set_this_lat
type (ring_struct), pointer :: to_this_lat
real(rp), pointer :: set_this_data(:)
real(rp), pointer :: to_this_data(:)

err = .false.

select case (set_lat_name)
  case ('model')
    set_this_lat => u%model
    set_this_data => u%data%model_value
  case ('base')
    set_this_lat => u%base
    set_this_data => u%data%base_value
  case ('design')
    set_this_lat => u%design
    set_this_data => u%data%design_value
  case default
    call out_io (s_error$, r_name, 'BAD LATTICE: ' // set_lattice)
    err = .true.
    return
end select

select case (to_lattice)
  case ('model')
    ! make sure model data is up to date
    s%global%lattice_recalc = .true.
    call tao_lattice_calc ()
    to_this_lat => u%model
    to_this_data => u%data%model_value
  case ('base')
    to_this_lat => u%base
    to_this_data => u%data%base_value
  case ('design')
    to_this_lat => u%design
    to_this_data => u%data%design_value
  case default
    call out_io (s_error$, r_name, 'BAD LATTICE: ' // to_lattice)
    err = .true.
    return
end select
  
set_this_lat = to_this_lat

set_this_data = to_this_data

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
!    %global  -- Global variables structure.
!-

subroutine tao_set_global_cmd (who, set_value)

implicit none

type (tao_global_struct) global

character(*) who, set_value
character(20) :: r_name = 'tao_set_global_cmd'

integer iu, ios

namelist / params / global

! open a scratch file for a namelist read

iu = lunget()
open (iu, status = 'scratch')
write (iu, *) '&params'
write (iu, *) ' global%' // trim(who) // ' = ' // trim(set_value)
write (iu, *) '/'
rewind (iu)
global = s%global  ! set defaults
read (iu, nml = params, iostat = ios)
close (iu)

if (ios == 0) then
  s%global = global
  if (trim(who) .eq. 'track_type') s%global%lattice_recalc = .true.
else
  call out_io (s_error$, r_name, 'BAD COMPONENT OR NUMBER')
endif

end subroutine tao_set_global_cmd

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
character(20) :: r_name = 'tao_set_plot_page_cmd'

real(rp) x, y
integer ix

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

type (tao_curve_struct), pointer :: curve

integer i, ios, i_uni
integer, allocatable :: ix_ele(:)

character(*) curve_name, component, set_value
character(20) :: r_name = 'tao_set_curve_cmd'

logical err

!

call tao_find_plot_by_region (err, curve_name, curve = curve)
if (err) return

i_uni = curve%ix_universe
if (i_uni == 0) i_uni = s%global%u_view

select case (component)

  case ('ele2_name')
    curve%ele2_name = set_value
    call tao_locate_element (curve%ele2_name, i_uni, ix_ele, .true.)
    curve%ix_ele2 = ix_ele(1)

  case ('ix_ele2')
    read (set_value, '(i)', iostat = ios) i
    if (ios /= 0) then
      call out_io (s_error$, r_name, 'BAD IX_ELE2 VALUE')
      return
    endif
    curve%ix_ele2 = i      
    curve%ele2_name = s%u(i_uni)%model%ele_(curve%ix_ele2)%name

  case default
    
    call out_io (s_error$, r_name, "BAD CURVE COMPONENT")
    return
    
end select

s%global%lattice_recalc = .true.

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_set_var_cmd (name, component, set_value, list)
!
! Routine to set var values.
!
! Input:
!   name       -- Character(*): Which var name to set.
!   component  -- Character(*): Which component to set.
!   set_value  -- Character(*): What value to set it to.
!   list       -- Character(*): If not blank then gives which indexes to apply to.
!
!  Output:
!-

subroutine tao_set_var_cmd (name, component, set_value, list)

use tao_mod
use quick_plot

implicit none

type (tao_v1_var_struct), pointer :: v1_ptr

integer i, j

character(*) name, component, set_value, list
character(20) :: r_name = 'tao_set_var_cmd'

logical err

! Find the var to set

if (name == 'all') then
  call set_this_var (s%var, .false.)
else
  call tao_find_var(err, name, v1_ptr)  
  if (err) return
  call set_this_var (v1_ptr%v, .true.)
endif


!----------------------------------------------------------------------------
contains

subroutine set_this_var (var, can_list)

type (tao_var_struct), target :: var(:)

real(rp), pointer :: r_ptr
real(rp) value

integer n_set, n1, n2, iv, ix, ios

character(1) using

logical, pointer :: l_ptr
logical, allocatable :: set_it(:)
logical multiply, divide, can_list, change_lattice

change_lattice = .false.

! parse the list of vars to set.
! The result is the set_it(:) array.

n1 = var(1)%ix_v1
n2 = n1 + size(var) - 1
allocate (set_it(n1:n2))

if (list == ' ') then
  set_it = .true.
elseif (.not. can_list) then
  call out_io (s_error$, r_name, 'A LIST DOES NOT MAKE SENSE HERE.')
  err = .true.
  return
else
  call location_decode (list, set_it, n1, n_set)
endif

! loop over all vars and do the sets.

do iv = 1, size(var)

  ix = var(iv)%ix_v1  ! input index.
  if (.not. set_it(ix)) cycle
  if (.not. var(iv)%exists) cycle

! select component

  select case (component)
  case ('weight')
    r_ptr => var(iv)%weight
    using = 'r'
  case ('step')
    r_ptr => var(iv)%step
    using = 'r'
  case ('model')
    r_ptr => var(iv)%model_value
    using = 'r'
    change_lattice = .true.
  case ('base')
    r_ptr => var(iv)%base_value
    using = 'r'
  case ('good_user')
    l_ptr => var(iv)%good_user
    using = 'l'
  case ('good_var')
    l_ptr => var(iv)%good_var
    using = 'l'
  case default
    err = .true.
    call out_io (s_error$, r_name, 'UNKNOWN COMPONENT NAME: ' // component)
    return
  end select

! select value and set.

  select case (set_value)
  case ('meas')
    call check_using (using, 'r', err); if (err) return
    r_ptr = var(iv)%meas_value
  case ('ref')
    call check_using (using, 'r', err); if (err) return
    r_ptr = var(iv)%ref_value
  case ('model')
    call check_using (using, 'r', err); if (err) return
    r_ptr = var(iv)%model_value
  case ('base')
    call check_using (using, 'r', err); if (err) return
    r_ptr = var(iv)%base_value
  case ('design')
    call check_using (using, 'r', err); if (err) return
    r_ptr = var(iv)%design_value
  case ('f', 'F')
    call check_using (using, 'l', err); if (err) return
    l_ptr = .false.
  case ('t', 'T')
    call check_using (using, 'l', err); if (err) return
    l_ptr = .true.
  case default
    call check_using (using, 'r', err); if (err) return
    multiply = .false.
    divide = .false.
    ix = 1
    if (set_value(1:1) == '*') then
      multiply = .true.
      ix = 2
    elseif (set_value(1:1) == '/') then
      divide = .true.
      ix = 2
    endif
    read (set_value(ix:), *, iostat = ios) value
    if (ios /= 0 .or. set_value(ix:) == ' ') then
      err = .true.
      call out_io (s_error$, r_name, 'BAD VALUE: ' // set_value)
      return
    endif
    if (multiply) then
      r_ptr = r_ptr * value
    elseif (divide) then
      r_ptr = r_ptr / value
    else
      r_ptr = value
    endif
  end select

! the model lattice has changed
  if (change_lattice) then
    call check_using (using, 'r', err); if (err) return
    call tao_set_var_model_value (var(iv), r_ptr)
  endif
    
enddo

! cleanup

deallocate (set_it)

end subroutine set_this_var

!----------------------------------------------------------------------------
! contains

subroutine check_using (using, should_be_using, err)
character(1) using, should_be_using
logical err

err = .false.
if (using /= should_be_using) then
  err = .true.
  call out_io (s_error$, r_name, 'VARIABLE COMPONENT/SET_VALUE TYPE MISMATCH.')
endif

end subroutine check_using

end subroutine tao_set_var_cmd
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_set_data_cmd (name, component, set_value, list)
!
! Routine to set data values.
!
! Input:
!   name       -- Character(*): Which data name to set.
!   component  -- Character(*): Which component to set.
!   set_value  -- Character(*): What value to set it to.
!   list       -- Character(*): If not blank then gives which indexes to apply to.
!
!  Output:
!-

subroutine tao_set_data_cmd (name, component, set_value, list)

use tao_mod
use quick_plot

implicit none


integer i

character(*) name, component, set_value, list
character(20) :: r_name = 'tao_set_data_cmd'

logical, automatic :: this_u(size(s%u))
logical err

!

call tao_pick_universe (name, name, this_u, err)
if (err) return

do i = 1, size(s%u)
  if (.not. this_u(i)) cycle
  call set_data (s%u(i))
  if (err) return
enddo

!----------------------------------------------------------------------------
contains

subroutine set_data (u)

type (tao_universe_struct), target :: u
type (tao_d2_data_struct), pointer :: d2_ptr
type (tao_d1_data_struct), pointer :: d1_ptr

integer j

! Find the data to set

if (name == 'all') then
  call set_this_data (u%data, .false.)
else
  call tao_find_data(err, u, name, d2_ptr, d1_ptr)  
  if (err) return
  if (associated(d1_ptr)) then
    call set_this_data (d1_ptr%d, .true.)
  else
    do j = 1, size(d2_ptr%d1)
      call set_this_data (d2_ptr%d1(j)%d, .true.)
      if (err) return
    enddo
  endif
endif

call tao_set_data_useit_opt ()

end subroutine set_data

!----------------------------------------------------------------------------
! contains

subroutine set_this_data (data, can_list)

type (tao_data_struct), target :: data(:)

real(rp), pointer :: r_ptr(:)
real(rp) value

integer n_set, n1, n2, ix, ios

character(1) using

logical, pointer :: l_ptr(:)
logical, allocatable :: set_it(:), good(:)
logical multiply, divide, can_list

! parse the list of data to set

err = .true.

n1 = data(1)%ix_d1
n2 = n1 + size(data) - 1
allocate (set_it(n1:n2), good(n1:n2))

if (list == ' ') then
  set_it = .true.
elseif (.not. can_list) then
  call out_io (s_error$, r_name, 'A LIST DOES NOT MAKE SENSE HERE.')
  return
else
  call location_decode (list, set_it, n1, n_set)
  if (n_set < 0) return
endif

set_it = set_it .and. data(:)%exists

! select component

select case (component)
case ('weight')
  r_ptr => data(:)%weight
  using = 'r'
case ('meas')
  r_ptr => data(:)%meas_value
  using = 'r'
case ('ref')
  r_ptr => data(:)%ref_value
  using = 'r'
case ('good_user')
  l_ptr => data(:)%good_user
  using = 'l'
case ('good_meas')
  l_ptr => data(:)%good_meas
  using = 'l'
case default
  call out_io (s_error$, r_name, 'UNKNOWN COMPONENT NAME: ' // component)
  return
end select

! select value and set.

select case (set_value)
case ('meas')
  call check_using (using, 'r', err); if (err) return
  where (set_it) r_ptr = data(:)%meas_value
  good = data(:)%good_meas
case ('ref')
  call check_using (using, 'r', err); if (err) return
  where (set_it) r_ptr = data(:)%ref_value
  good = data(:)%good_ref
case ('model')
  call check_using (using, 'r', err); if (err) return
  where (set_it) r_ptr = data(:)%model_value
  good = data(:)%exists
case ('base')
  call check_using (using, 'r', err); if (err) return
  where (set_it) r_ptr = data(:)%base_value
  good = data(:)%exists
case ('design')
  call check_using (using, 'r', err); if (err) return
  where (set_it) r_ptr = data(:)%design_value
  good = data(:)%exists
case ('f', 'F')
  call check_using (using, 'l', err); if (err) return
  where (set_it) l_ptr = .false.
case ('t', 'T')
  call check_using (using, 'l', err); if (err) return
  where (set_it) l_ptr = .true.
case default
  call check_using (using, 'r', err); if (err) return
  multiply = .false.
  divide = .false.
  ix = 1
  if (set_value(1:1) == '*') then
    multiply = .true.
    ix = 2
  elseif (set_value(1:1) == '/') then
    divide = .true.
    ix = 2
  endif
  read (set_value(ix:), *, iostat = ios) value
  if (ios /= 0 .or. set_value(ix:) == ' ') then
    call out_io (s_error$, r_name, 'BAD VALUE: ' // set_value)
    deallocate (set_it)
    return
  endif
  if (multiply) then
    where (set_it) r_ptr = r_ptr * value
  elseif (divide) then
    where (set_it) r_ptr = r_ptr / value
  else
    where (set_it) r_ptr = value
  endif
  good = data(:)%exists
end select

! set good

select case (component)
case ('meas')
  where (set_it) data(:)%good_meas = good
case ('ref')
  where (set_it) data(:)%good_ref = good
end select

! cleanup

err = .false.
deallocate (set_it)

end subroutine set_this_data

!----------------------------------------------------------------------------
! contains

subroutine check_using (using, should_be_using, err)
character(1) using, should_be_using
logical err

err = .false.
if (using /= should_be_using) then
  err = .true.
  call out_io (s_error$, r_name, 'VARIABLE COMPONENT/SET_VALUE TYPE MISMATCH.')
endif

end subroutine check_using

end subroutine tao_set_data_cmd

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine tao_set_uni_cmd (uni, on_off, recalc)
!
! turns a universe of or off
!
! Input:
!  uni       -- Integer: which universe; 0 => all off (if you really want to)
!  on_off    -- Character(*): "on" or "off"
!  recalc    -- Logical: Recalculate lattices
!
! Output:
!  s%u(uni)%is_on
!
!-

subroutine tao_set_uni_cmd (uni, on_off, recalc)

implicit none

integer uni, i

character(*) on_off
character(20) :: r_name = "tao_set_universe_cmd"

logical is_on, recalc

  call str_upcase (on_off, on_off)

  if (on_off(1:2) .eq. 'ON') then
    is_on = .true.
  elseif (on_off(1:3) .eq. 'OFF') then
    is_on = .false.
  else
    call out_io (s_warn$, r_name, &
                 "Syntax Error: Can only turn universe 'on' or 'off'")
    return
  endif

  if (uni .lt. 0 .or. uni .gt. size(s%u)) then
    call out_io (s_warn$, r_name, &
                 "Invalid Universe specifier")
    return
  endif
  
  if (uni .eq. 0) then
    call out_io (s_blank$, r_name, &
        "Changing all universes!")
    s%u(:)%is_on = is_on
  else
    s%u(uni)%is_on = is_on
  endif

  ! make sure lattice calculation is up to date if turning lattice on
  if (recalc) s%global%lattice_recalc = .true.
  
end subroutine tao_set_uni_cmd


end module tao_set_mod
