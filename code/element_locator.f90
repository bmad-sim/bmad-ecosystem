!+
! Subroutine ELEMENT_LOCATOR (ELE_NAME, RING, IX_ELE)
!
! Subroutine to locate an element in a ring. 
!
! Modules Needed:
!   use bmad
!
! Input:
!   ELE_NAME -- Character*16: Name of the element to find.
!   RING     -- Ring_struct: Ring to search through
!
! Output:
!   IX_ELE -- Integer: Index of element. IX_ELE set to -1 if not found.
!-

!$Id$
!$Log$
!Revision 1.4  2002/02/23 20:32:15  dcs
!Double/Single Real toggle added
!
!Revision 1.3  2002/01/08 21:44:38  dcs
!Aligned with VMS version  -- DCS
!
!Revision 1.2  2001/09/27 18:31:51  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine element_locator (ele_name, ring, ix_ele)

  use bmad
  
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
