!+
! Function pointer_to_girder(ele) result (girder, ix_slave_back)
!
! Routine to determine if an element is being supported by a girder element.
!
! Input:
!   ele            -- ele_struct: Element to check.
!
! Output:
!   ix_slave_back -- integer, optional: Index back to ele. That is, pointer_to_slave(girder, ix_slave_back) 
!                        will point back to ele. Set to -1 if no girder present
!   girder        -- ele_struct, pointer :: Pointer to the girder. Null if ele is not girder supported.
!- 

function pointer_to_girder (ele, ix_slave_back) result (girder)

use equal_mod, dummy => pointer_to_girder

implicit none

type (ele_struct), target :: ele
type (ele_struct), pointer :: girder
integer, optional :: ix_slave_back
integer i

!

do i = 1, ele%n_lord
  girder => pointer_to_lord(ele, i, ix_slave_back = ix_slave_back)
  if (girder%key /= girder$) cycle
  return
enddo

if (present(ix_slave_back)) ix_slave_back = -1
nullify(girder)

end function pointer_to_girder
