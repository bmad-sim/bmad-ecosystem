!+
! Subroutine make_mat6_tracking (ele, param, c0, c1)
!
! Subroutine to make the 6x6 transfer matrix for an element using the
! Present tracking method.
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

subroutine make_mat6_tracking (ele, param, c0, c1)

  use bmad

  implicit none

  type (ele_struct), target :: ele
  type (coord_struct) :: c0, c1
  type (param_struct)  param

  real (rdef) error
  integer temp_method
  logical temp_symplectify

!

  temp_symplectify = ele%symplectify
  ele%symplectify = .false.   ! don't do this twice.

  call transfer_mat_from_tracking (ele, param, c0, bmad_com%d_orb, c1, error)

  ele%symplectify = temp_symplectify

end subroutine

