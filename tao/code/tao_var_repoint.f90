!+
! Suburoutine tao_var_repoint()
!
! Routine to reset tao variable pointers.
!-

subroutine tao_var_repoint ()

use tao_interface, dummy => tao_var_repoint

implicit none

type (tao_var_struct), pointer :: var
type (tao_var_slave_struct), pointer :: sl
type (tao_universe_struct), pointer :: u
type (all_pointer_struct) :: a_ptr
type (all_pointer_struct), allocatable :: all_ptr(:)
type (ele_struct), pointer :: ele

integer i, j
logical err, print_err

! Remember: If the "-no_grouping" swithch has been used, there can be multiple elements with the same name
! but a variable, in this case, will only have one of them as a slave.
! This makes finding the correct slave element difficult if the slave element has shifted position. 
! Therefore, it is forbidden for any Tao command (for example the "read lattice" command) to modify element positions.

do i = 1, size(s%var)
  var => s%var(i)
  if (size(var%slave) == 0) cycle

  do j = 1, size(var%slave)
    sl => var%slave(j)
    u => s%u(sl%ix_uni)
    print_err = associated(sl%model_value)

    ! Something like 'PARTICLE_START' will not have an associated lattice element.
    if (sl%ix_ele < 0) then
      call pointers_to_attribute (u%model%lat, var%ele_name, var%attrib_name, .true., all_ptr, err, print_err)
      sl%model_value => all_ptr(1)%r

      call pointers_to_attribute (u%base%lat, var%ele_name, var%attrib_name, .true., all_ptr, err, print_err)
      sl%base_value => all_ptr(1)%r

    else
      ele => pointer_to_ele(u%model%lat, sl%ix_ele, sl%ix_branch)
      call pointer_to_attribute (ele, var%attrib_name, .true., a_ptr, err, print_err)
      sl%model_value => a_ptr%r

      ele => pointer_to_ele(u%base%lat, sl%ix_ele, sl%ix_branch)
      call pointer_to_attribute (ele, var%attrib_name, .true., a_ptr, err, print_err)
      sl%base_value => a_ptr%r
    endif

  enddo

  var%model_value => var%slave(1)%model_value
  var%base_value => var%slave(1)%base_value
enddo

end subroutine
