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

!$Id$
!$Log$
!Revision 1.2  2002/08/20 20:34:54  dcs
!symp_lie_bmad / symp_lie_ptc added
!
!Revision 1.1  2002/06/13 15:07:22  dcs
!Merged with FPP/PTC
!
!Revision 1.2  2001/09/27 18:31:50  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"

subroutine track1_custom (start, ele, param, end)

  use bmad

  implicit none

  type (coord_struct), intent(in)  :: start
  type (coord_struct), intent(out) :: end
  type (ele_struct),   intent(inout)  :: ele
  type (param_struct), intent(inout) :: param

!

  call track1_runge_kutta (start, ele, param, end)

end subroutine
