!+
! Subroutine ELEMENT_LOCATOR (ELE_NAME, RING, IX_ELE)
!
! Subroutine to locate an element in a ring. 
!
! Modules Needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!   ELE_NAME -- Character*16: Name of the element to find.
!   RING     -- Ring_struct: Ring to search through
!
! Output:
!   IX_ELE -- Integer: Index of element. IX_ELE set to -1 if not found.
!-

subroutine element_locator (ele_name, ring, ix_ele)

  use bmad_struct

  implicit none

  type (ring_struct) ring
  integer ix_ele
  character*(*) ele_name

!

  do ix_ele = 1, ring%n_ele_max
    if (ring%ele_(ix_ele)%name == ele_name) return
  enddo

  ix_ele = -1

end subroutine
