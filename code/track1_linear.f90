!+
! Subroutine track1_linear (start, ele, param, end)
!
! Particle tracking through a single element assuming linearity.
! That is, just using ele%mat6.
!
! Modules Needed:
!   use bmad
!
! Input:
!   start  -- Coord_struct: Starting position
!   ele    -- Ele_struct: Element
!   param  -- Param_struct:
!
! Output:
!   end   -- Coord_struct: End position
!   param -- Param_struct:
!-

#include "CESR_platform.inc"

subroutine track1_linear (start, ele, param, end)

  use bmad_struct
  use bmad_interface

  implicit none

  type (coord_struct), intent(in)  :: start
  type (coord_struct), intent(out) :: end
  type (ele_struct),   intent(inout)  :: ele
  type (param_struct), intent(inout) :: param

! 

  end%vec = matmul (ele%mat6, start%vec) + ele%vec0

end subroutine
