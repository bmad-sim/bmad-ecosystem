!+
! Subroutine make_mat6_symp_lie_ptc (ele, param, c0, c1)
!
! Subroutine to make the 6x6 transfer matrix for an element. 
!
! Modules needed:
!   use bmad
!
! Input:
!   ele    -- Ele_struct: Element with transfer matrix
!   param  -- Param_struct: Parameters are needed for some elements.
!   c0     -- Coord_struct: Coordinates at the beginning of element. 
!
! Output:
!   ele    -- Ele_struct: Element with transfer matrix.
!     %mat6  -- 6x6 transfer matrix.
!   c1     -- Coord_struct: Coordinates at the end of element.
!-

#include "CESR_platform.inc"

subroutine make_mat6_symp_lie_ptc (ele, param, c0, c1)

  use ptc_interface_mod

  implicit none

  type (ele_struct), target :: ele
  type (coord_struct) :: c0, c1
  type (param_struct)  param

!

  call ele_to_taylor(ele, c0, param)
  call taylor_to_mat6 (ele%taylor, c0%vec, ele%mat6, c1%vec)


end subroutine

