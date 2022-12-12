!+
! Function orbit_to_local_curvilinear (orbit, ele, z_direction) result (local_position)
!
! Routine to return the local curvilinear position and orientation of a particle.
! The orbit is in laboratory coords so element misalignments do not affect the calculation.
!
! Input:
!   orbit         -- coord_struct: Particle orbit in laboratory (not body) coordinates.
!   ele           -- ele_struct: Lattice element particle is in.
!   z_direction   -- integer, optional: Set to +1 or -1.  Z-direction of particle velocity
!                     relative to element z-axis. Default is ele%orientation * orbit%direction.
!
! Output:
!   local_position  -- floor_position_struct: Position in local coordinates.
!-

function orbit_to_local_curvilinear (orbit, ele, z_direction) result (local_position)

use bmad_routine_interface, dummy => orbit_to_local_curvilinear

implicit none

type (coord_struct) orbit
type (ele_struct) ele
type (floor_position_struct) local_position

real(rp) px, py, pz, r2
integer, optional :: z_direction

!

select case (ele%orientation)
case (+1);  local_position%r = [orbit%vec(1), orbit%vec(3), orbit%s - ele%s_start]
case (-1);  local_position%r = [orbit%vec(1), orbit%vec(3), ele%s - orbit%s]
end select

px = orbit%vec(2) / (1 + orbit%vec(6))
py = orbit%vec(4) / (1 + orbit%vec(6))
r2 = 1.0_rp - px**2 - py**2
pz = integer_option(ele%orientation * orbit%direction, z_direction) * sqrt(max(0.0_rp, r2))

local_position%phi = atan2(py, sqrt(px**2  + py**2))
local_position%theta = atan2(px, pz)
call floor_angles_to_w_mat (local_position%theta, local_position%phi, 0.0_rp, local_position%w)

end function orbit_to_local_curvilinear

