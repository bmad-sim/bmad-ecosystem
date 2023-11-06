!+
! Subroutine track_a_thick_multipole (orbit, ele, param, mat6, make_matrix)
!
! Bmad_standard tracking through a thick multipole element (sextupole, kicker, etc). 
!
! Input:
!   orbit       -- Coord_struct: Starting position.
!   ele         -- ele_struct: Thick multipole element.
!   param       -- lat_param_struct: Lattice parameters.
!   mat6(6,6)   -- Real(rp), optional: Transfer matrix before the element.
!   make_matrix -- logical, optional: Propagate the transfer matrix? Default is false.
!
! Output:
!   orbit      -- coord_struct: End position.
!   mat6(6,6)  -- real(rp), optional: Transfer matrix through the element.
!-

subroutine track_a_thick_multipole (orbit, ele, param, mat6, make_matrix)

use bmad_interface, except_dummy => track_a_thick_multipole

implicit none

type (coord_struct) :: orbit, start_orb
type (ele_struct), target :: ele
type (lat_param_struct) :: param
type (fringe_field_info_struct) fringe_info

real(rp), optional :: mat6(6,6)
real(rp) rel_tracking_charge, charge_dir, r_step, step_len, s_off, mass, length
real(rp) an(0:n_pole_maxx), bn(0:n_pole_maxx), an_elec(0:n_pole_maxx), bn_elec(0:n_pole_maxx)
real(rp) rtc, hk, vk, kick, angle_E, k_E, beta_ref, mc2, ac_amp, s_pos

integer i, n_step, orientation, ix_mag_max, ix_elec_max

logical, optional :: make_matrix
logical drifting, include_kicks

!

start_orb = orbit
length = ele%value(l$) * start_orb%time_dir
orientation = ele%orientation * start_orb%direction
rel_tracking_charge = rel_tracking_charge_to_mass(start_orb, param%particle)
charge_dir = rel_tracking_charge * orientation
include_kicks = .true.
kick = 0

if (ele%key == elseparator$) then
  include_kicks = .false.
  rtc = abs(rel_tracking_charge) * sign_of(charge_of(orbit%species))
  beta_ref = ele%value(p0c$) / ele%value(e_tot$)
  mc2 = mass_of(orbit%species)

  hk = ele%value(hkick$) * rtc
  vk = ele%value(vkick$) * rtc
  kick = sqrt(hk**2 + vk**2) 

  angle_E = atan2(vk, hk)

  if (length == 0) then
    k_E = 1  ! Something non-zero
  else
    k_E = kick / ele%value(l$)
  endif
endif

!

call multipole_ele_to_ab (ele, .false., ix_mag_max, an,      bn,      magnetic$, include_kicks$)
call multipole_ele_to_ab (ele, .false., ix_elec_max, an_elec, bn_elec, electric$)

if (kick == 0 .and. ix_mag_max == -1 .and. ix_elec_max == -1) then
  n_step = 1
else
  n_step = max(nint(abs(length / ele%value(ds_step$))), 1)
endif
r_step = real(orbit%time_dir, rp) / n_step
step_len = length / n_step

! Entrance edge

call offset_particle (ele, set$, orbit, set_hvkicks = .false., mat6 = mat6, make_matrix = make_matrix)

call init_fringe_info (fringe_info, ele)
if (fringe_info%has_fringe) then
  fringe_info%particle_at = first_track_edge$
  call apply_element_edge_kick(orbit, fringe_info, ele, param, .false., mat6, make_matrix)
  if (orbit%state /= alive$) return
endif

! Multipole kicks. Notice that the magnetic multipoles have already been normalized by the length.

ac_amp = ac_kicker_amp(ele, orbit)

if (ix_mag_max > -1)  call ab_multipole_kicks (an,      bn,      ix_mag_max,  ele, orbit, magnetic$, ac_amp*r_step/2,   mat6, make_matrix)
if (ix_elec_max > -1) call ab_multipole_kicks (an_elec, bn_elec, ix_elec_max, ele, orbit, electric$, ac_amp*step_len/2, mat6, make_matrix)

! Body

do i = 1, n_step

  if (ele%key == elseparator$ .and. kick /= 0) then
    call track_this_elsep(step_len)
    if (orbit%state /= alive$) return
  else
    call track_a_drift (orbit, step_len, mat6, make_matrix)
  endif

  s_pos = i * step_len / n_step
  ac_amp = ac_kicker_amp(ele, orbit)

  if (i == n_step) then
    if (ix_mag_max > -1)  call ab_multipole_kicks (an,      bn,      ix_mag_max,  ele, orbit, magnetic$, ac_amp*r_step/2,   mat6, make_matrix)
    if (ix_elec_max > -1) call ab_multipole_kicks (an_elec, bn_elec, ix_elec_max, ele, orbit, electric$, ac_amp*step_len/2, mat6, make_matrix)
  else
    if (ix_mag_max > -1)  call ab_multipole_kicks (an,      bn,      ix_mag_max,  ele, orbit, magnetic$, ac_amp*r_step,   mat6, make_matrix)
    if (ix_elec_max > -1) call ab_multipole_kicks (an_elec, bn_elec, ix_elec_max, ele, orbit, electric$, ac_amp*step_len, mat6, make_matrix)
  endif

enddo

! Exit edge

if (fringe_info%has_fringe) then
  fringe_info%particle_at = second_track_edge$
  call apply_element_edge_kick(orbit, fringe_info, ele, param, .false., mat6, make_matrix)
  if (orbit%state /= alive$) return
endif

call offset_particle (ele, unset$, orbit, set_hvkicks = .false., mat6 = mat6, make_matrix = make_matrix)

orbit%t = start_orb%t + orbit%direction*orbit%time_dir*ele%value(delta_ref_time$) + (start_orb%vec(5) - orbit%vec(5)) / (orbit%beta * c_light)

!---------------------------------------------
contains

subroutine track_this_elsep (step_len)

real(rp) step_len, sinh_k, cosh_k, cosh1_k, dx, dct, pc, E_tot, E_rel, px, py, p_factor, ps, ps2, alpha, coef
real(rp) kmat(6,6), dErel_dpz, dcoef_dpx, dcoef_dpy, dcoef_dpz, dbeta_dpz
real(rp) beta, pz, z1, E_rel_new, E_new, dbeta_dE, f

!

call tilt_coords (angle_E, orbit%vec, mat6, make_matrix)

beta = orbit%beta
px = orbit%vec(2)
py = orbit%vec(4)
pz = orbit%vec(6)
pc = ele%value(p0c$) * (1 + pz)

E_tot = pc / beta
E_rel = E_tot / ele%value(p0c$)

p_factor = (mc2 / ele%value(p0c$))**2 + px**2 + py**2
if (E_rel**2 < p_factor) then
  orbit = start_orb
  orbit%state = lost_pz$
  return
endif

! Track

ps2 = E_rel**2 - p_factor
ps = sqrt(ps2)
coef = k_E * step_len / ps

if (abs(coef) > 10) then ! lost
  orbit%state = lost$
  return
endif

if (abs(coef) < 1d-3) then
  alpha = step_len / ps
  sinh_k = alpha * (1 + coef**2 / 6 + coef**4/120)
  cosh1_k = alpha * coef * (1.0_rp / 2 + coef**2 / 24 + coef**4 / 720)
  cosh_k = cosh(coef) / k_E
else
  cosh_k = cosh(coef) / k_E
  sinh_k = sinh(coef) / k_E
  cosh1_k = cosh_k - 1 / k_E
endif

dx = E_rel * cosh1_k + px * sinh_k
orbit%vec(1) = orbit%vec(1) + dx
orbit%vec(2) = k_E * (E_rel * sinh_k + px * cosh_k)
orbit%vec(3) = orbit%vec(3) + step_len * py / ps

dct = (E_rel * sinh_k + px * cosh1_k)
orbit%t = orbit%t + dct / c_light
orbit%vec(5) = orbit%vec(5) + orbit%beta * (step_len / beta_ref - dct)
z1 = orbit%vec(5)

call apply_energy_kick (dx*k_E*ele%value(p0c$), orbit, [0.0_rp, 0.0_rp])

if (logic_option(.false., make_matrix)) then
  dErel_dpz = (1 + orbit%vec(6)) / E_rel
  dcoef_dpx = px * coef / ps2
  dcoef_dpy = py * coef / ps2
  dcoef_dpz = -E_rel * dErel_dpz * coef / ps2
  dbeta_dpz = 1 / (E_rel * (E_tot / mc2)**2) 
  E_rel_new = (1 + orbit%vec(6)) / orbit%beta
  E_new = E_rel_new * orbit%p0c
  dbeta_dE = 1 / ((1 + orbit%vec(6)) * orbit%p0c * (E_new / mc2)**2)
  f = k_E * ele%value(p0c$)

  call mat_make_unit(kmat)
  kmat(1,2) = (E_rel * sinh_k + px * cosh_k) * dcoef_dpx + sinh_k 
  kmat(1,4) = (E_rel * sinh_k + px * cosh_k) * dcoef_dpy
  kmat(1,6) = (E_rel * sinh_k + px * cosh_k) * dcoef_dpz + dErel_dpz * cosh1_k

  kmat(2,2) = k_E * (E_rel * cosh_k + px * sinh_k) * dcoef_dpx + k_E * cosh_k
  kmat(2,4) = k_E * (E_rel * cosh_k + px * sinh_k) * dcoef_dpy
  kmat(2,6) = k_E * (E_rel * cosh_k + px * sinh_k) * dcoef_dpz + k_E * dErel_dpz * sinh_k

  kmat(3,2) = step_len * py * px / (ps * ps2)
  kmat(3,4) = step_len * (py * py / (ps * ps2) + 1 / ps)
  kmat(3,6) = -step_len * py * E_rel * dErel_dpz / (ps * ps2)

  kmat(5,2) = -orbit%beta * ((E_rel * cosh_k + px * sinh_k) * dcoef_dpx + cosh1_k) + &
                                                            z1 * dbeta_dE * f * kmat(1,2) / beta
  kmat(5,4) = -orbit%beta *  (E_rel * cosh_k + px * sinh_k) * dcoef_dpy + &
                                                            z1 * dbeta_dE * f * kmat(1,4) / beta
  kmat(5,5) = orbit%beta / beta
  kmat(5,6) = -orbit%beta * ((E_rel * cosh_k + px * sinh_k) * dcoef_dpz + dErel_dpz * sinh_k) + &
                dbeta_dpz * (step_len / beta_ref - dct) + &
                z1 * dbeta_dE * (f * kmat(1,6) + beta * orbit%p0c) / beta - z1 * orbit%beta * dbeta_dpz / beta**2 

  kmat(6,2) = f * kmat(1,2) / (orbit%p0c * orbit%beta)
  kmat(6,4) = f * kmat(1,4) / (orbit%p0c * orbit%beta)
  kmat(6,6) = (f * kmat(1,6) + orbit%p0c * beta) / (orbit%p0c * orbit%beta)

  mat6 = matmul(kmat, mat6)
endif

call tilt_coords (-angle_E, orbit%vec, mat6, make_matrix)

end subroutine track_this_elsep

end subroutine
