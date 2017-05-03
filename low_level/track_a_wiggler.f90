!+
! Subroutine track_a_wiggler (orbit, ele, param, mat6, make_matrix)
!
! Bmad_standard tracking through a wiggler or undulator element.
!
! Input:
!   orbit       -- Coord_struct: Starting position.
!   ele         -- ele_struct: Wiggler element.
!   param       -- lat_param_struct: Lattice parameters.
!   make_matrix -- logical, optional: Propagate the transfer matrix? Default is false.
!
! Output:
!   orbit      -- coord_struct: End position.
!   mat6(6,6)  -- real(rp), optional: Transfer matrix through the element.
!-

subroutine track_a_wiggler (orbit, ele, param, mat6, make_matrix)

use track1_mod, except_dummy => track_a_wiggler

implicit none

type (coord_struct) :: orbit
type (ele_struct), target :: ele
type (lat_param_struct) :: param

real(rp), optional :: mat6(6,6)
real(rp) mat2(2,2), z_start, beta_ref, p_factor, k1_factor, k1, k1l, length, k_z, charge_dir, rel_p, t_start

logical, optional :: make_matrix

! Only periodic type wigglers are handled here.
! In the horizontal plane the tracking looks like a drift.
! The tracking in the vertical plane is:
!   1) 1/2 the octupole kick at the entrance face.
!   2) Track as a quadrupole through the body
!   3) 1/2 the octupole kick at the exit face.

z_start = orbit%vec(5)
t_start = orbit%t
beta_ref = ele%value(p0c$) / ele%value(e_tot$)
length = ele%value(l$)
charge_dir = rel_tracking_charge_to_mass(orbit, param) * ele%orientation
rel_p = 1 + orbit%vec(6)

call offset_particle (ele, param, set$, orbit)

if (ele%value(l_pole$) == 0) then
  k_z = 1d100    ! Something large
else
  k_z = pi / ele%value(l_pole$)
endif

k1_factor = -charge_dir * 0.5 * (c_light * ele%value(b_max$) / ele%value(p0c$))**2
k1 = k1_factor / rel_p**2
k1l = length * k1

p_factor = 1 - (orbit%vec(2) / rel_p**2)**2 - (orbit%vec(4) / rel_p**2)**2
if (p_factor < 0) then
  orbit%state = lost_z_aperture$
  return
endif

p_factor = sqrt(p_factor)
orbit%vec(5) = orbit%vec(5) + 0.5 * length * (orbit%beta / beta_ref - 1/p_factor) - &
                                    k1_factor * (1/(c_light * beta_ref) - 1 / (c_light * orbit%beta * rel_p**2)) / (4 * k_z**2)
orbit%t = orbit%t + length / (2 * c_light * orbit%beta * p_factor) + k1l / (4 * k_z**2 * c_light * orbit%beta)

! 1/2 of the octupole kick at the entrance face.

orbit%vec(4) = orbit%vec(4) + k1l * rel_p * k_z**2 * orbit%vec(3)**3 / 3

! Quadrupole body

call quad_mat2_calc (k1, length, rel_p, mat2)
orbit%vec(1) = orbit%vec(1) + length * orbit%vec(2) / rel_p
orbit%vec(3:4) = matmul (mat2, orbit%vec(3:4))

! 1/2 of the octupole kick at the exit face.

orbit%vec(4) = orbit%vec(4) + k1l * rel_p * k_z**2 * orbit%vec(3)**3 / 3

p_factor = 1 - (orbit%vec(2) / rel_p**2)**2 - (orbit%vec(4) / rel_p**2)**2
if (p_factor < 0) then
  orbit%state = lost_z_aperture$
  return
endif

p_factor = sqrt(p_factor)
orbit%vec(5) = orbit%vec(5) + 0.5 * length * (orbit%beta / beta_ref - 1/p_factor) - &
                                    k1_factor * (1/(c_light * beta_ref) - 1 / (c_light * orbit%beta * rel_p**2)) / (4 * k_z**2)
orbit%t = orbit%t + length / (2 * c_light * orbit%beta * p_factor) + k1l / (4 * k_z**2 * c_light * orbit%beta)

call offset_particle (ele, param, unset$, orbit)

end subroutine
