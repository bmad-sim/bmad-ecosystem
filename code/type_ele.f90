!+
! Subroutine type_ele (ele, type_zero_attrib, type_mat6, type_twiss, 
!                                                           type_control, ring)
!
! Subroutine to type out information on an element. 
! See also the subroutine type2_ele.
!
! Modules needed:
!   use bmad
!
! Input:
!   ELE              -- Ele_struct: Element
!   TYPE_ZERO_ATTRIB -- Logical: If true then type all attributes even if the
!                          attribute value is 0.
!   TYPE_MAT6        -- Integer:
!                          TYPE_MAT6 = 0   => Do not type ELE%MAT6
!                          TYPE_MAT6 = 4   => Type 4X4 XY submatrix
!                          TYPE_MAT6 = 6   => Type full 6x6 matrix
!   TYPE_TWISS       -- Logical: If true then type the twiss parameters
!                          at the end of the element.
!   TYPE_CONTROL     -- Logical: If true then type control status.
!   RING             -- Ring_struct, Optional: Needed for control typeout.
!
!-
!$Id$
!$Log$
!Revision 1.5  2002/02/23 20:32:30  dcs
!Double/Single Real toggle added
!
!Revision 1.4  2001/11/29 19:39:54  helms
!Updates from DCS including (*) -> (:)
!
!Revision 1.3  2001/10/12 20:53:35  rwh24
!DCS changes and two files added
!
!Revision 1.2  2001/09/27 18:32:01  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine type_ele (ele, type_zero_attrib, type_mat6, type_twiss,  &
                                                       type_control, ring)

  use bmad_struct
  use bmad_interface, only: type2_ele

  implicit none

  type (ele_struct)  ele
  type (ring_struct), optional :: ring

  integer type_mat6, n_lines, i
  logical type_twiss, type_control, type_zero_attrib

  character*80 lines(50)

!

  if (type_twiss) then
    call type2_ele (ele, type_zero_attrib, type_mat6, radians$, &
                                          type_control, lines, n_lines, ring)
  else
    call type2_ele (ele, type_zero_attrib, type_mat6, 0, &
                                          type_control, lines, n_lines, ring)
  endif

  do i = 1, n_lines
    print '(1x, a)', trim(lines(i))
  enddo

end subroutine
