!+
! Suburoutine tao_var_repoint()
!
! Routine to reset tao variable pointers.
!-

subroutine tao_var_repoint ()

use tao_struct

implicit none

type (tao_var_struct), pointer :: var_ptr
type (tao_var_slave_struct), pointer :: sl
type (tao_universe_struct), pointer :: u
type (all_pointer_struct) :: a_ptr
type (ele_struct), pointer :: ele

integer i, j
logical err

! Remember: If the "-no_grouping" swithch has been used, there can be multiple elements with the same name
! but a variable, in this case, will only have one of them as a slave.
! This makes finding the correct slave element difficult if the slave element has shifted position. 
! Therefore, it is forbidden for any Tao command (for example the "read lattice" command) to modify element positions.

do i = 1, size(s%var)
  var_ptr => s%var(i)
  if (size(var_ptr%slave) == 0) cycle

  do j = 1, size(var_ptr%slave)
    sl => var_ptr%slave(j)
    u => s%u(sl%ix_uni)

    ele => pointer_to_ele(u%model%lat, sl%ix_ele, sl%ix_branch)
    call pointer_to_attribute (ele, var_ptr%attrib_name, .true., a_ptr, err, .false.)
    sl%model_value => a_ptr%r

    ele => pointer_to_ele(u%base%lat, sl%ix_ele, sl%ix_branch)
    call pointer_to_attribute (ele, var_ptr%attrib_name, .true., a_ptr, err, .false.)
    sl%base_value => a_ptr%r
  enddo

  var_ptr%model_value => var_ptr%slave(1)%model_value
  var_ptr%base_value => var_ptr%slave(1)%base_value
enddo

end subroutine
