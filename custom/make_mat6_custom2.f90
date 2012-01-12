!+
! Subroutine make_mat6_custom2 (ele, param, c0, c1)
!
! Dummy subroutine for custom calculation of the transfer matrix. 
! If called, this routine will generate an error message and quit.
! This routine needs to be replaced for a custom calculation.
!
! Note: This routine is not to be confused with make_mat6_custom.
! See the Bmad manual for more details.
!
! General rule: Your code may NOT modify any argument that is not listed as
! an output agument below."
!
! Modules needed:
!   use bmad
!
! Input:
!   ele    -- Ele_struct: Element with transfer matrix
!   param  -- lat_param_struct: Parameters are needed for some elements.
!   c0     -- Coord_struct: Coordinates at the beginning of element. 
!
! Output:
!   ele    -- Ele_struct: Element with transfer matrix.
!     %mat6  -- 6x6 transfer matrix.
!   c1     -- Coord_struct: Coordinates at the end of element.
!+

#include "CESR_platform.inc"

subroutine make_mat6_custom2 (ele, param, c0, c1)

  use bmad_struct
  use bmad_interface, except_dummy => make_mat6_custom2

  implicit none

  type (ele_struct), target :: ele
  type (coord_struct) :: c0, c1
  type (lat_param_struct)  param

!

  print *, 'ERROR: DUMMY MAKE_MAT6_CUSTOM2 CALLED FOR: ', ele%name
  print *, '       EITHER CUSTOM MAT6_CALC_METHOD WAS CALLED BY MISTAKE,'
  print *, '       OR THE CORRECT ROUTINE WAS NOT LINKED IN!'
  call err_exit

end subroutine
