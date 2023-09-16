!+
! Function coords_relative_to_floor (floor0, dr, theta, phi, psi) result (floor1)
!
! Starting from a given reference frame specified by its orientation and
! position in the global (floor) coordinates, and given a shift in position
! and angular orientation with respect to this reference frame, return the 
! resulting reference frame orientation and position.
!
! Also see: coords_floor_to_relative
!
! Input:
!   floor0   -- floor_position_struct: Initial reference frame.
!   dr(3)    -- real(rp): (x, y, z) positional shift of the reference frame.
!   theta, phi, psi
!            -- real(rp), optional: Angular shift of the reference frame. See the 
!                 Bmad manual on the Global Coordinate system for more details.
!                 All angles must either be absent or present.
!
! Output:
!   floor1   -- floor_position_struct: Shifted reference frame.
!-

function coords_relative_to_floor (floor0, dr, theta, phi, psi) result (floor1)

use bmad_interface, dummy => coords_relative_to_floor

implicit none

type (floor_position_struct) floor0, floor1
real(rp) dr(3)
real(rp), optional :: theta, phi, psi
real(rp) w_mat(3,3), w0_mat(3,3)

!
floor1%r = matmul(floor0%w, dr) + floor0%r

if (present(theta)) then
  call floor_angles_to_w_mat (theta, phi, psi, w_mat)
  floor1%W = matmul(floor0%W, w_mat)
  call update_floor_angles (floor1, floor0)
endif

end function coords_relative_to_floor
