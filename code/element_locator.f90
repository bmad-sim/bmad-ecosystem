!+
! Subroutine element_locator (ele_name, lat, ix_ele)
!
! Subroutine to locate an element in a lattice. 
!
! Modules Needed:
!   use bmad
!
! Input:
!   ele_name -- Character(40): Name of the element to find.
!   lat     -- lat_struct: Lattice to search through.
!
! Output:
!   ix_ele -- Integer: Index of element in lat%ele(:) array. 
!               ix_ele set to -1 if not found.
!-

#include "CESR_platform.inc"

subroutine element_locator (ele_name, lat, ix_ele)

  use bmad_struct
  use bmad_interface, except_dummy => element_locator
  
  implicit none

  type (lat_struct) lat
  integer ix_ele
  character(*) ele_name

! Very simple.

  do ix_ele = 0, lat%n_ele_max
    if (lat%ele(ix_ele)%name == ele_name) return
  enddo

  ix_ele = -1

end subroutine
