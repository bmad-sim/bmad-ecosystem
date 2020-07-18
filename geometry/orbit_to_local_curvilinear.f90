!+
! Function orbit_to_local_curvilinear (orbit, ele_s_start, in_ele_frame) result (local_position)
!
! Routine to return the local curvilinear position and orientation of a particle.
!
! Input:
!   orbit         -- coord_struct: Particle orbit.
!   ele_s_start   -- real(rp): S_start of the element the particle is in.
!   in_ele_frame  -- logical, optional: True => orbit is in ele body frame and includes misalignments.
!                      Default is False.
!
! Output:
!   local_position  -- floor_position_struct: Position in local coordinates.
!-

function orbit_to_local_curvilinear (orbit, ele_s_start, in_ele_frame) result (local_position)

use bmad_routine_interface, dummy => orbit_to_local_curvilinear

implicit none

type (coord_struct) orbit
type (floor_position_struct) local_position

real(rp) ele_s_start
real(rp) px, py, phi, theta
logical, optional :: in_ele_frame

!

local_position%r = [orbit%vec(1), orbit%vec(3), orbit%s - ele_s_start]
px = orbit%vec(2) / (1 + orbit%vec(6))
py = orbit%vec(4) / (1 + orbit%vec(6))

local_position%phi = asin(py)
local_position%theta = asin(px / cos(phi))
call floor_angles_to_w_mat (local_position%theta, local_position%phi, 0.0_rp, local_position%w)

end function orbit_to_local_curvilinear

