!+
! Function tao_pointer_to_datum (d1, ele_name) result (datum_ptr)
!
! Routine to return a pointer to a datum using the lattice element name.
!
! Input:
!   d1        -- tao_d1_data_struct: D1 data struct to search.
!   ele_name  -- character(*): Name of lattice element to match to.
!
! Output:
!   datum_ptr   -- data_struct, pointer: Pointer to the matched datum. 
!                   Will be null if no match found.
!-

function tao_pointer_to_datum (d1, ele_name) result (datum_ptr)

use tao_interface, except_dummy => tao_pointer_to_datum

implicit none

type (tao_d1_data_struct), target :: d1
type (tao_data_struct), pointer :: datum_ptr

integer i
character(*) ele_name

!

do i = lbound(d1%d, 1), ubound(d1%d, 1)
  if (d1%d(i)%ele_name /= ele_name) cycle
  datum_ptr => d1%d(i)
  return
enddo

nullify(datum_ptr)

end function
