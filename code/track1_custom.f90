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
!     %aperture_limit_on -- If .true. then %lost will be set if the
!                 particle is outsile the aperture.
!
! Output:
!   end   -- Coord_struct: End position
!   param
!     %lost -- Set .true. If the particle is outside the aperture and
!                %aperture_limit_on is set. Also: %LOST is set .true. if
!                the particle does not make it through a bend irregardless
!                of the the setting of %aperture_limit_on.
!-

!$Id$
!$Log$
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
