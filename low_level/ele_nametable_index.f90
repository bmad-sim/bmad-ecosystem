!+
! Function ele_nametable_index(ele) result(ix_nt)
!
! Routine to return the index in the nametable corresponding to ele.
! The reverse routine is: pointer_to_ele.
!
! Input:
!   ele       -- ele_struct: Element in a lattice.
!
! Output:
!   ix_nt     -- integer: Nametable index. lat%nametable%name(ix_nt) and lat%nametable%index(ix_nt)
!                 correspond with ele. 
!                 Set to -1 if ele is not a lattice element. For example, a slice_slave is not a lattice element.
!-

function ele_nametable_index(ele) result(ix_nt)

use bmad_routine_interface, dummy => ele_nametable_index

implicit none

type (ele_struct), target :: ele
type (lat_struct), pointer :: lat
integer ix_nt, ib

!

ix_nt = -1
if (ele%ix_ele < 0) return
if (.not. associated(ele%branch)) return

lat => ele%branch%lat
ix_nt = ele%ix_ele
do ib = 0, ele%ix_branch-1
  ix_nt = ix_nt + lat%branch(ib)%n_ele_max + 1
enddo

end function ele_nametable_index

