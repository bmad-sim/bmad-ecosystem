!+
! Subroutine track_a_elseparator (orbit, ele, param, mat6, make_matrix)
!
! Particle tracking through a elseparator element. 
!
! Modules Needed:
!   use bmad
!
! Input:
!   orbit       -- Coord_struct: Starting position.
!   ele         -- Ele_struct: Elseparator element.
!   param       -- lat_param_struct: Lattice parameters.
!   make_matrix -- logical, optional: Propagate the transfer matrix? Default is false.
!
! Output:
!   orbit      -- Coord_struct: End position.
!   mat6(6,6)  -- Real(rp), optional: Transfer matrix through the element.
!-

subroutine track_a_elseparator (orbit, ele, param, mat6, make_matrix)

use track1_mod, except_dummy => track_a_elseparator

implicit none

type (coord_struct) :: orbit, start_orb
type (ele_struct), target :: ele
type (lat_param_struct) :: param
type (fringe_edge_info_struct) fringe_info

real(rp), optional :: mat6(6,6)
real(rp) an(0:n_pole_maxx), bn(0:n_pole_maxx), an_elec(0:n_pole_maxx), bn_elec(0:n_pole_maxx)
real(rp) hk, vk, rtc, angle_E, pc, E_tot, E_rel, mc2, px, p_factor, k_E, ps, alpha, coef, cos_E, sin_E, dx
real(rp) sinh_k, cos1_k, dt, beta_ref, rel_tracking_charge, charge_dir, length, cosh1_k, r_step, step_len, kick

integer orientation, n_step, ix_pole_max, ix_elec_max

logical, optional :: make_matrix

!

start_orb = orbit

beta_ref = ele%value(p0c$) / ele%value(e_tot$)
mc2 = mass_of(orbit%species)
orientation = ele%orientation * start_orb%direction
rel_tracking_charge = rel_tracking_charge_to_mass(start_orb, param)
charge_dir = rel_tracking_charge * orientation
length = ele%value(l$)
if (length == 0) length = 1d-50  ! To avoid divide by zero

rtc = abs(rel_tracking_charge) * sign(1, charge_of(orbit%species))

hk = ele%value(hkick$) * rtc
vk = ele%value(vkick$) * rtc
kick = sqrt(hk**2 + vk**2) 

angle_E = atan2(vk, hk)

if (kick == 0) then
  cos_E = 1
  sin_E = 0
else
  cos_E = hk / kick
  sin_E = vk / kick
endif

k_E = kick / length

call multipole_ele_to_ab (ele, .false., ix_pole_max, an,      bn,      magnetic$)
call multipole_ele_to_ab (ele, .false., ix_elec_max, an_elec, bn_elec, electric$)

!

if (logic_option(.false., make_matrix)) call mat_make_unit(mat6)

call offset_particle (ele, param, set$, orbit, set_hvkicks = .false., set_multipoles = .false., mat6 = mat6, make_matrix = make_matrix)

nullify(fringe_info%hard_ele)
fringe_info%particle_at = first_track_edge$
call apply_element_edge_kick(orbit, fringe_info, ele, param, .false., mat6, make_matrix)
if (orbit%state /= alive$) return

!

n_step = 1
if (ix_pole_max > -1 .or. ix_elec_max > -1) n_step = max(nint(ele%value(l$) / ele%value(ds_step$)), 1)

r_step = 1.0_rp / n_step
step_len = ele%value(l$) * r_step

!


pc = ele%value(p0c$) * (1 + orbit%vec(6))
call convert_pc_to (pc, orbit%species, E_tot = E_tot)
E_rel = E_tot / ele%value(p0c$)

call tilt_coords (angle_E, orbit%vec, mat6, make_matrix)
px = orbit%vec(2)
p_factor = (mc2 / ele%value(p0c$))**2 + orbit%vec(2)**2 + orbit%vec(4)**2

if (E_rel**2 < p_factor) then
  orbit = start_orb
  orbit%state = lost_z_aperture$
  return
endif

! Track

ps = sqrt(E_rel**2 - p_factor)
alpha = length / ps
coef = k_E * length / ps

if (abs(coef) > 10) then ! lost
  orbit%state = lost$
  return
endif

if (abs(coef) < 1d-3) then
  sinh_k = alpha * (1 + coef**2 / 6 + coef**4/120)
  cosh1_k = alpha * coef * (1.0_rp / 2 + coef**2 / 24 + coef**4 / 720)
else
  sinh_k = sinh(coef) / k_E
  cosh1_k = (cosh(coef) - 1) / k_E
endif

dx = E_rel * cosh1_k + px * sinh_k
orbit%vec(1) = orbit%vec(1) + dx
orbit%vec(2) = E_rel * sinh(coef) + px * cosh(coef)
orbit%vec(3) = orbit%vec(3) + length * orbit%vec(4) / ps

dt = (E_rel * sinh_k + px * cosh1_k) / c_light
orbit%t = orbit%t + dt
orbit%vec(5) = orbit%vec(5) + orbit%beta * (length / beta_ref - c_light * dt)

call apply_energy_kick (dx*k_E*ele%value(p0c$), orbit, [k_E, 0.0_rp], mat6, make_matrix)

call tilt_coords (-angle_E, orbit%vec, mat6, make_matrix)

fringe_info%particle_at = second_track_edge$
call apply_element_edge_kick(orbit, fringe_info, ele, param, .false., mat6, make_matrix)
if (orbit%state /= alive$) return

call offset_particle (ele, param, unset$, orbit, set_hvkicks = .false., set_multipoles = .false., mat6 = mat6, make_matrix = make_matrix)

end subroutine
