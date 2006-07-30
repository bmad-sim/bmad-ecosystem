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
  use bmad_interface, except => track1_linear

  implicit none

  type (coord_struct) :: start
  type (coord_struct) :: end
  type (ele_struct) :: ele
  type (param_struct) :: param

! 

  end%vec = matmul (ele%mat6, start%vec) + ele%vec0

end subroutine
