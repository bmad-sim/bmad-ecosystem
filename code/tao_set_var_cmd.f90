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
!----------------------------------------------------------------------------
contains

subroutine set_this_var (var, can_list)

type (tao_var_struct), target :: var(:)

real(rp), pointer :: c_ptr
real(rp) value

integer n_set, n1, n2, iv, ix, ios

logical, allocatable :: set_it(:)
logical multiply, divide, can_list

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
    c_ptr => var(iv)%weight
  case ('step')
    c_ptr => var(iv)%step
  case ('model')
    c_ptr => var(iv)%model_value
  case ('base')
    c_ptr => var(iv)%base_value
  case default
    err = .true.
    call out_io (s_error$, r_name, 'UNKNOWN COMPONENT NAME: ' // component)
  end select

! select value and set.

  select case (set_value)
  case ('meas')
    c_ptr = var(iv)%meas_value
  case ('ref')
    c_ptr = var(iv)%ref_value
  case ('model')
    c_ptr = var(iv)%model_value
  case ('base')
    c_ptr = var(iv)%base_value
  case ('design')
    c_ptr = var(iv)%design_value
  case default
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
      c_ptr = c_ptr * value
    elseif (divide) then
      c_ptr = c_ptr / value
    else
      c_ptr = value
    endif
  end select

enddo

! cleanup

deallocate (set_it)

end subroutine 

end subroutine

