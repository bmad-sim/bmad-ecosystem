!+
! Subroutine RING_GEOMETRY (RING)
!
! Subroutine to calculate the physical placement of all the elements in a ring.
! That is, the layout on the floor.
!
! Modules Needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!     RING -- ring_struct: The ring
!
! Output:
!     RING
!       .ELE_(I)
!         .X_POSITION        -- X position at end of element
!         .Y_POSITION        -- Y position at end of element
!         .Z_POSITION        -- Z position at end of element
!         .THETA_POSITION    -- Orientation angle at end of element in X-Y
!
! Note: The starting point is taken to be the position in RING.ELE_(0)
!-

!$Id$
!$Log$
!Revision 1.2  2001/09/27 18:31:57  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"



subroutine ring_geometry (ring)

  use bmad_struct
  implicit none

  type (ring_struct)  ring
  type (ele_struct)  ele

  integer i

  real*8 theta_tot, x_pos, y_pos, chord_len, angle, leng
  real*8 z_pos

!

  theta_tot = 0.0D0
  x_pos = 0.0D0
  y_pos = 0.0D0
  z_pos = 0.0D0

  do i = 1, ring%n_ele_ring

    ele = ring%ele_(i)
    leng = dble(ele%value(l$))

    if (ele%key == sbend$ .or. ele%key == rbend$) then
      angle = dble(ele%value(angle$))
      chord_len = 2 * leng * sin(angle/2) / angle
    else
      chord_len = leng
      angle = 0.0D0
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
