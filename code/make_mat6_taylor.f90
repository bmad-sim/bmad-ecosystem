!+
! Subroutine make_mat6_taylor (ele, param, c0)
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
!-

#include "CESR_platform.inc"

subroutine make_mat6_taylor (ele, param, c0)

  use ptc_interface_mod, except => make_mat6_taylor

  implicit none

  type (ele_struct), target :: ele
  type (coord_struct) :: c0
  type (param_struct)  param

!

  if (.not. associated(ele%taylor(1)%term)) call ele_to_taylor(ele, param, c0)
  call taylor_to_mat6 (ele%taylor, c0%vec, ele%mat6)

end subroutine

