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

logical err

!

do i = 1, size(s%u)
  call set_data (s%u(i))
  if (err) return
enddo

!----------------------------------------------------------------------------
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

end subroutine

!----------------------------------------------------------------------------
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
case default
  call out_io (s_error$, r_name, 'UNKNOWN COMPONENT NAME: ' // component)
  return
end select

! select value and set.

select case (set_value)
case ('meas')
  call check_using (using, 'r', err); if (err) return
  where (set_it) r_ptr = data(:)%meas_value
  good = data(:)%good_data
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
  where (set_it) data(:)%good_data = good
case ('ref')
  where (set_it) data(:)%good_ref = good
end select

! cleanup

err = .false.
deallocate (set_it)

end subroutine 

!----------------------------------------------------------------------------
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

end subroutine

end subroutine
