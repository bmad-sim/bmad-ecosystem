!+
! Subroutine track1_custom (start, ele, param, end)
!
! Default routine for custom_tracking. 
! This routine will do Runge Kutta tracking.
! You need to supply field_rk_custom.
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
!   param
!-

#include "CESR_platform.inc"

subroutine track1_custom (start, ele, param, end)

  use runge_kutta_mod
  use bmad_interface

  implicit none

  type (coord_struct), intent(in)  :: start
  type (coord_struct), intent(out) :: end
  type (ele_struct),   intent(inout)  :: ele
  type (param_struct), intent(inout) :: param

!

  call track1_runge_kutta (start, ele, param, end)

end subroutine
