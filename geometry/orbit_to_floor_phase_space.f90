!+
! Function orbit_to_floor_phase_space (orbit, ele) result (floor_phase_space)
!
! Routine to return the global floor phase space coordinates of a particle:
!       (x, y, z, vx, vy, vz)
! where v_vec = (vx, vy, vz) is the normalized velocity with |v_vec| = 1.
!
! The orbit is in laboratory coords so element misalignments do not affect the calculation.
!
! Input:
!   orbit         -- coord_struct: Particle orbit in local (not element) coordinates.
!   ele           -- ele_struct: Lattice element particle is in.
!
! Output:
!   floor_phase_space(6)  -- real(rp): Floor phase space
!-

function orbit_to_floor_phase_space (orbit, ele) result (floor_phase_space)

use bmad_routine_interface, dummy => orbit_to_floor_phase_space

implicit none

type (coord_struct) orbit
type (ele_struct) ele
type (floor_position_struct) local_position, floor

real(rp) floor_phase_space(6), v(3), w_mat(3,3)

!

local_position%r = [orbit%vec(1), orbit%vec(3), orbit%s - ele%s_start]

v(1) = orbit%vec(2) / (1 + orbit%vec(6))
v(2) = orbit%vec(4) / (1 + orbit%vec(6))
v(3) = ele%orientation * orbit%direction * sqrt(max(0.0_rp, 1.0_rp - v(1)**2 - v(2)**2))

floor = coords_local_curvilinear_to_floor(local_position, ele, .false., w_mat)
floor_phase_space(1:3) = floor%r
floor_phase_space(4:6) = matmul(w_mat, v)

end function orbit_to_floor_phase_space

