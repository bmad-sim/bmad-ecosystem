!+
! Subroutine set_ele_name(ele, name)
!
! Routine to set a lattice element name and update the lat%nametable as appropriate.
!
! Input:
!   ele     -- ele_struct: Element whose name is to be set.
!   name    -- character(*): Name to set.
!
! Output:
!   ele     -- ele_struct: Element with name set.
!-

subroutine set_ele_name(ele, name)

use bmad_routine_interface, dummy => set_ele_name

implicit none

type (ele_struct) ele
character(*) name
integer ix

!

ele%name = name

if (.not.associated(ele%branch)) return
ix = ele_nametable_index(ele)
call nametable_change1(ele%branch%lat%nametable, ele%name, ix)

end subroutine set_ele_name
