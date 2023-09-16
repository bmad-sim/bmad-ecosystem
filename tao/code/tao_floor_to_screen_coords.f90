!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_floor_to_screen_coords (graph, floor, screen)
!
! Routine to project a 3D floor coordinate onto a 2D projection plane
! along with projecting the orientation vector.
!
! Input:
!   graph  -- tao_graph_struct: Graph defining the projection plane.
!   floor  -- floor_position_struct: 3D coordinate.
!
! Output:
!   screen  -- floor_position_struct: Projected point
!     %r(3)   -- projected (x, y) = (%r(1), %r(2))
!     %theta  -- angle in (x, y) plane of projection of the orientation vector.
!-

subroutine tao_floor_to_screen_coords (graph, floor, screen)

use tao_struct

implicit none

type (tao_graph_struct) graph
type (floor_position_struct) floor, screen
real(rp) orient(3), theta, phi, x, y

! Get projection position

call tao_floor_to_screen (graph, floor%r, screen%r(1), screen%r(2))

! screen%theta does not depend upon floor%psi

theta = floor%theta
phi = floor%phi
orient = [sin(theta) * cos(phi), sin(phi), cos(theta) * cos(phi)]  ! orientation vector

select case (graph%floor_plan%view(1:1))
case ('x')
  x = orient(1)
case ('y')
  x = orient(2)
case ('z')
  x = orient(3)
end select

select case (graph%floor_plan%view(2:2))
case ('x')
  y = orient(1)
case ('y')
  y = orient(2)
case ('z')
  y = orient(3)
end select

screen%theta = atan2(y, x) + twopi * graph%floor_plan%rotation

end subroutine tao_floor_to_screen_coords
