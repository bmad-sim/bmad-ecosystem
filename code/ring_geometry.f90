!+
! Subroutine ring_geometry (ring)
!
! Subroutine to calculate the physical placement of all the elements in a ring.
! That is, the layout on the floor. 
! This is the same as the MAD convention. See the MAD manual for more details.
!
! Note: At present this routine assumes no  vertical bends. That is, 
! y_position is always 0 and the ring in in the X-Z plane.
!
! Modules Needed:
!   use bmad
!
! Input:
!   ring -- ring_struct: The ring
!
! Output:
!   ring
!     %ele_(i)
!       %x_position        -- X position at end of element
!       %y_position        -- Y position at end of element
!       %z_position        -- Z position at end of element
!       %theta_position    -- Orientation angle at end of element in X-Z plane
!
! Note: The starting point is taken to be the position in RING%ELE_(0)
!-

#include "CESR_platform.inc"

subroutine ring_geometry (ring)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ring_struct)  ring
  type (ele_struct)  ele

  integer i

  real(8) theta_tot, x_pos, y_pos, chord_len, angle, leng
  real(8) z_pos

!

  theta_tot = ring%ele_(0)%theta_position
  x_pos = ring%ele_(0)%x_position
  y_pos = ring%ele_(0)%y_position
  z_pos = ring%ele_(0)%z_position

  do i = 1, ring%n_ele_ring

    if (ele%key == patch$) then
      print *, 'ERROR IN RING_GEOMETRY: PATCH NOT YET IMPLEMENTED!'
      call err_exit      
      cycle
    endif

    ele = ring%ele_(i)
    leng = dble(ele%value(l$))

    if (ele%key == sbend$ .or. ele%key == rbend$) then
      angle = dble(ele%value(l$)) * dble(ele%value(g$))

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

    theta_tot = theta_tot - angle / 2
    x_pos = x_pos + chord_len * sin(theta_tot)
    z_pos = z_pos + chord_len * cos(theta_tot)
    theta_tot = theta_tot - angle / 2

    ring%ele_(i)%theta_position = theta_tot
    ring%ele_(i)%x_position = x_pos
    ring%ele_(i)%y_position = 0
    ring%ele_(i)%z_position = z_pos

  enddo

end subroutine
