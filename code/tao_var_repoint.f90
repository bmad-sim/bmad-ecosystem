!+
! Suburoutine tao_var_repoint()
!
! Routine to reset tao variable pointers.
! This is needed if the lattice is perturbed.
! For example, this routine needs to be called after a read command (which calls bmad_parser2) since the
! lattice file called may use superposition.
!-

subroutine tao_var_repoint ()

use tao_struct

implicit none

type (tao_var_struct), pointer :: var_ptr
type (tao_var_slave_struct), pointer :: sl
type (tao_universe_struct), pointer :: u
type (all_pointer_struct), allocatable :: a_ptr(:), b_ptr(:)

integer i, j, num
logical err

!

do i = 1, size(s%var)
  var_ptr => s%var(i)
  num = 0
  if (size(var_ptr%slave) == 0) cycle

  do 
    if (num == size(var_ptr%slave)) exit
    sl => var_ptr%slave(num+1)
    u => s%u(sl%ix_uni)
    call pointers_to_attribute (u%model%lat, var_ptr%ele_name, var_ptr%attrib_name, .true., a_ptr, err, .false.)
    call pointers_to_attribute (u%model%lat, var_ptr%ele_name, var_ptr%attrib_name, .true., b_ptr, err, .false.)
    do j = 1, size(a_ptr)
      var_ptr%slave(num+j)%model_value => a_ptr(j)%r
      var_ptr%slave(num+j)%base_value => b_ptr(j)%r
    enddo
    num = num + size(a_ptr)
  enddo

  var_ptr%model_value => var_ptr%slave(1)%model_value
  var_ptr%base_value => var_ptr%slave(1)%base_value
enddo

end subroutine
