!+
! Subroutine ring_geometry (ring)
!
! Subroutine to calculate the physical placement of all the elements in a ring.
! That is, the layout on the floor. Note: At present this routine assumes no 
! vertical bends.
!
! Modules Needed:
!   use bmad
!
! Input:
!   ring -- ring_struct: The ring
!
! Output:
!     ring
!       %ele_(i)
!         %x_position        -- X position at end of element
!         %y_position        -- Y position at end of element
!         %z_position        -- Z position at end of element
!         %theta_position    -- Orientation angle at end of element in X-Y
!
! Note: The starting point is taken to be the position in RING%ELE_(0)
!-

!$Id$
!$Log$
!Revision 1.5  2002/06/13 14:54:28  dcs
!Interfaced with FPP/PTC
!
!Revision 1.4  2002/02/23 20:32:23  dcs
!Double/Single Real toggle added
!
!Revision 1.3  2002/01/08 21:44:42  dcs
!Aligned with VMS version  -- DCS
!
!Revision 1.2  2001/09/27 18:31:57  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"

subroutine ring_geometry (ring)

  use bmad
  implicit none

  type (ring_struct)  ring
  type (ele_struct)  ele

  integer i

  real*8 theta_tot, x_pos, y_pos, chord_len, angle, leng
  real*8 z_pos

!

  theta_tot = 0D0
  x_pos = 0.0D0
  y_pos = 0.0D0
  z_pos = 0.0D0

  do i = 1, ring%n_ele_ring

    ele = ring%ele_(i)
    leng = dble(ele%value(l$))

    if (ele%key == sbend$ .or. ele%key == rbend$) then
      angle = dble(ele%value(l$)) * dble(ele%value(g_design$))

      if (angle == 0) then
        chord_len = 0
      else
        chord_len = 2 * leng * sin(angle/2) / angle
      endif

      if (abs(modulo(ele%value(tilt$), twopi) - pi) < 1e-6) then
        angle = -angle
      elseif (ele%value(tilt$) /= 0) then
        print *, 'ERROR IN RING_GEOMETRY: NON-PLAINER GEOMETRY NOT YET IMPLEMENTED!'
        print *, '      FOR BEND: ', ele%name
        print *, '      WITH TILT: ', ele%value(tilt$)
      endif

    else
      angle = 0.0D0
      chord_len = leng
    endif

    theta_tot = theta_tot + angle / 2
    x_pos = x_pos + chord_len * cos(theta_tot)
    y_pos = y_pos + chord_len * sin(theta_tot)
    theta_tot = theta_tot + angle / 2

    ring%ele_(i)%theta_position = theta_tot
    ring%ele_(i)%x_position = x_pos
    ring%ele_(i)%y_position = y_pos

  enddo

  return
  end
