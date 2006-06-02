!+
! Subroutine element_locator (ele_name, ring, ix_ele)
!
! Subroutine to locate an element in a lattice. 
!
! Modules Needed:
!   use bmad
!
! Input:
!   ele_name -- Character(40): Name of the element to find.
!   ring     -- Ring_struct: Lattice to search through.
!
! Output:
!   ix_ele -- Integer: Index of element in ring%ele_(:) array. 
!               ix_ele set to -1 if not found.
!-

#include "CESR_platform.inc"

subroutine element_locator (ele_name, ring, ix_ele)

  use bmad_struct
  use bmad_interface, except => element_locator
  
  implicit none

  type (ring_struct) ring
  integer ix_ele
  character(*) ele_name

! Very simple.

  do ix_ele = 0, ring%n_ele_max
    if (ring%ele_(ix_ele)%name == ele_name) return
  enddo

  ix_ele = -1

end subroutine
