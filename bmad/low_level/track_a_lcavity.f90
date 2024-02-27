!+
! Subroutine track_a_lcavity (orbit, ele, param, mat6, make_matrix)
!
! Bmad_standard tracking through a lcavity element.
!
! Modified version of:
!       J. Rosenzweig and L. Serafini
!       Phys Rev E, Vol. 49, p. 1599, (1994)
! with b_0 = b_-1 = 1. See the Bmad manual for more details.
!
! One must keep in mind that we are NOT using good canonical coordinates since
!   the energy of the reference particle is changing.
! This means that the resulting matrix will NOT be symplectic.
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

subroutine track_a_lcavity (orbit, ele, param, mat6, make_matrix)

use bmad_interface, except_dummy => track_a_lcavity
use super_recipes_mod, only: super_zbrent

implicit none

type (coord_struct) :: orbit
type (ele_struct), target :: ele
type (lat_param_struct) :: param
type (em_field_struct) field

! Used to get around ifort bug preventing, when debugging, seeing vars in contained code.
type internal_state_struct
  real(rp) gradient_max, step_len, E_start, E_end, phase0, dphase, pc_start, pc_end
  real(rp) beta_start, beta_end, cdt_ref
end type

type (internal_state_struct), target :: iss
type (internal_state_struct), pointer :: issp

real(rp), optional :: mat6(6,6)
real(rp) length, coef, k2
real(rp) alpha, sin_a, cos_a, r_mat(2,2), dph, E_ref_end
real(rp) phase, voltage_max, sqrt_8, f, k1
real(rp) sqrt_beta12, dsqrt_beta12(6), f_ave, E1, E2
real(rp) pxy2, xp1, xp2, yp1, yp2, mc2, om, om_g, m2(2,2), kmat(6,6), ds, r_step
real(rp) c_min, c_plu, dc_min, dc_plu, cos_term, drp1_dr0, drp1_drp0, drp2_dr0, drp2_drp0
real(rp) an(0:n_pole_maxx), bn(0:n_pole_maxx), an_elec(0:n_pole_maxx), bn_elec(0:n_pole_maxx)
real(rp) E_ref_start, p0c_start, p0c_end, phase1, phase2, ph_err(2), rf_phase
real(rp), parameter :: phase_abs_tol = 1e-4_rp

integer i, ix_mag_max, ix_elec_max, n_step, status

logical, optional :: make_matrix

character(*), parameter :: r_name = 'track_a_lcavity'

! 

issp => iss

if (ele%value(rf_frequency$) == 0  .and. (ele%value(voltage$) /= 0 .or. ele%value(voltage_err$) /= 0)) then
  call out_io (s_error$, r_name, 'LCAVITY ELEMENT HAS ZERO RF_FREQUENCY: ' // ele%name)
  orbit%state = lost$
  return
endif

length = orbit%time_dir * ele%value(l$)
if (length == 0) return
n_step = max(1, nint(ele%value(num_steps$)))
r_step = rp8(orbit%time_dir) / n_step
iss%step_len = length / n_step

call multipole_ele_to_ab (ele, .false., ix_mag_max, an,      bn,      magnetic$, include_kicks$)
call multipole_ele_to_ab (ele, .false., ix_elec_max, an_elec, bn_elec, electric$)

call offset_particle (ele, set$, orbit, mat6 = mat6, make_matrix = make_matrix)

if (ix_mag_max > -1)  call ab_multipole_kicks (an,      bn,      ix_mag_max,  ele, orbit, magnetic$, 0.5_rp*r_step,   mat6, make_matrix)
if (ix_elec_max > -1) call ab_multipole_kicks (an_elec, bn_elec, ix_elec_max, ele, orbit, electric$, 0.5_rp*iss%step_len, mat6, make_matrix)

if (orbit%time_dir * orbit%direction == 1) then
  E_ref_start = ele%value(E_tot_start$)
  E_ref_end   = ele%value(E_tot$)
  p0c_start   = ele%value(p0c_start$)
  p0c_end     = ele%value(p0c$)
else
  E_ref_start = ele%value(E_tot$)
  E_ref_end   = ele%value(E_tot_start$)
  p0c_start   = ele%value(p0c$)
  p0c_end     = ele%value(p0c_start$)
endif

mc2 = mass_of(orbit%species)

call reference_energy_correction(ele, orbit, first_track_edge$, mat6, make_matrix)

iss%pc_start = p0c_end * (1 + orbit%vec(6))
iss%beta_start = orbit%beta
iss%E_start = iss%pc_start / iss%beta_start 

! Traveling_wave fringe (standing_wave fringe is built-in to the body formulas)

if (nint(ele%value(cavity_type$)) == traveling_wave$ .and. fringe_here(ele, orbit, first_track_edge$)) then
  ds = bmad_com%significant_length / 10  ! Make sure inside field region
  call em_field_calc (ele, param, ds, orbit, .true., field, logic_option(.false., make_matrix))
  f = charge_of(orbit%species) / (2 * p0c_end)

  if (logic_option(.false., make_matrix)) then
    call mat_make_unit(kmat)
    kmat(2,1) = -f * (field%dE(3,1) * orbit%vec(1) + field%E(3))
    kmat(2,3) = -f * field%dE(3,2) * orbit%vec(1) 
    kmat(2,5) = -f * field%dE(3,3) * orbit%vec(1) * iss%beta_start
    kmat(2,6) =  f * field%E(3) * orbit%vec(1) * f / iss%pc_start
    kmat(4,1) = -f * field%dE(3,1) * orbit%vec(3)
    kmat(4,3) = -f * (field%dE(3,2) * orbit%vec(3) + field%E(3))
    kmat(4,5) = -f * field%dE(3,3) * orbit%vec(3) * iss%beta_start
    kmat(4,6) =  f * field%E(3) * orbit%vec(1) * f / iss%pc_start
    mat6 = matmul(kmat, mat6)
  endif

  orbit%vec(2) = orbit%vec(2) - f * field%E(3) * orbit%vec(1)
  orbit%vec(4) = orbit%vec(4) - f * field%E(3) * orbit%vec(3)
endif

do i = 1, n_step

  E1 = E_ref_start + (E_ref_end - E_ref_start) * real(i - 1, rp) / n_step
  E2 = E_ref_start + (E_ref_end - E_ref_start) * real(i, rp) / n_step
  iss%cdt_ref = iss%step_len * (E1 + E2) / (sqrt(E1**2 - mc2**2) + sqrt(E2**2 - mc2**2))

  ! The RF phase is defined with respect to the time at the beginning of the element.
  ! So if dealing with a slave element and absolute time tracking then need to correct.
  ! Note: phi0_autoscale is not used here since bmad_standard tracking by design gives the correct tracking.
  ! In fact, using phi0_autoscale would be a mistake if, say, tracking_method = runge_kutta, mat6_calc_method = bmad_standard.

  phase = twopi * (ele%value(phi0_err$) + ele%value(phi0$) + ele%value(phi0_multipass$) + &
             (particle_rf_time (orbit, ele, .false.) - rf_ref_time_offset(ele)) * ele%value(rf_frequency$))
  if (bmad_com%absolute_time_tracking .and. ele%orientation*orbit%time_dir*orbit%direction == -1) then
    phase = phase - twopi * ele%value(rf_frequency$) * ele%value(delta_ref_time$)
  endif
  phase = modulo2(phase, pi)

  if (i == 1) call rf_coupler_kick (ele, param, first_track_edge$, phase, orbit, mat6, make_matrix)

  !

  iss%gradient_max = e_accel_field(ele, gradient$, .true.)
  iss%phase0 = phase
  phase1 = phase

  ! The idea is to find the "average" phase so that forward tracking followed by backwards time 
  ! tracking comes back to the original position.
  ! If the phase change from start to finish is more than pi then this whole calculation is garbage.
  ! In this case, the particle is considered lost

  ph_err(1) = phase_func(phase1, status)
  if (abs(iss%dphase) > pi) then
    orbit%state = lost_pz$
    return

  elseif (abs(iss%dphase) < phase_abs_tol) then
    rf_phase = phase1 + 0.5_rp * iss%dphase
    call track_this_lcavity(rf_phase, iss, orbit, make_matrix)

  else
    phase2 = iss%phase0 + iss%dphase
    ph_err(2) = phase_func(phase2, status)
    ! At very low energies can happen that iss%phase0 and iss%phase0+iss%dphase do not bracket the solution.
    ! In this case try extrapolating. Factor of 2 to try to make sure root is bracketed.
    do
      if (ph_err(1)*ph_err(2) <= 0) exit
      dph = phase2 - phase1
      phase2 = phase1 + 2 * sign_of(dph) * max(abs(dph), 0.1)
      if (abs(phase2 - phase1) > pi) then
        orbit%state = lost_pz$
        return
      endif
      ph_err(2) = phase_func(phase2, status)
    enddo

    rf_phase = super_zbrent (phase_func, phase1, phase2, 0.0_rp, phase_abs_tol, status, ph_err)
    call track_this_lcavity(rf_phase, iss, orbit, make_matrix)
  endif

  if (i /= n_step) then
    if (ix_mag_max > -1)  call ab_multipole_kicks (an,      bn,      ix_mag_max,  ele, orbit, magnetic$, r_step,   mat6, make_matrix)
    if (ix_elec_max > -1) call ab_multipole_kicks (an_elec, bn_elec, ix_elec_max, ele, orbit, electric$, iss%step_len, mat6, make_matrix)
  endif

  iss%pc_start = iss%pc_end
  iss%beta_start = iss%beta_end
  iss%E_start = iss%E_end
enddo

!

if (nint(ele%value(cavity_type$)) == traveling_wave$ .and. fringe_here(ele, orbit, second_track_edge$)) then
  ds = bmad_com%significant_length / 10  ! Make sure inside field region
  call em_field_calc (ele, param, length - ds, orbit, .true., field, logic_option(.false., make_matrix))
  f = -charge_of(orbit%species) / (2 * p0c_end)

  if (logic_option(.false., make_matrix)) then
    call mat_make_unit(kmat)
    kmat(2,1) = -f * (field%dE(3,1) * orbit%vec(1) + field%E(3))
    kmat(2,3) = -f * field%dE(3,2) * orbit%vec(1) 
    kmat(2,5) = -f * field%dE(3,3) * orbit%vec(1) * iss%beta_end
    kmat(2,6) =  f * field%E(3) * orbit%vec(1) * f / iss%pc_end
    kmat(4,1) = -f * field%dE(3,1) * orbit%vec(3)
    kmat(4,3) = -f * (field%dE(3,2) * orbit%vec(3) + field%E(3))
    kmat(4,5) = -f * field%dE(3,3) * orbit%vec(3) * iss%beta_end
    kmat(4,6) =  f * field%E(3) * orbit%vec(1) * f / iss%pc_end
    mat6 = matmul(kmat, mat6)
  endif

  orbit%vec(2) = orbit%vec(2) - f * field%E(3) * orbit%vec(1)
  orbit%vec(4) = orbit%vec(4) - f * field%E(3) * orbit%vec(3)
endif

! Coupler kick

call rf_coupler_kick (ele, param, second_track_edge$, phase, orbit, mat6, make_matrix)

if (ix_elec_max > -1) call ab_multipole_kicks (an_elec, bn_elec, ix_elec_max, ele, orbit, electric$, 0.5_rp*iss%step_len, mat6, make_matrix)
if (ix_mag_max > -1)  call ab_multipole_kicks (an,      bn,      ix_mag_max,  ele, orbit, magnetic$, 0.5_rp*r_step,   mat6, make_matrix)

call offset_particle (ele, unset$, orbit, mat6 = mat6, make_matrix = make_matrix)

!---------------------------------------------------------------------------------------------
contains

! Returns rf_phase_at_end - rf_phase_at_start
! Note: rf_phase_at_start = iss%phase0 global variable
!       rf_phase is the phase used to calculate the accelerating gradient.

function dphase_end_minus_start(rf_phase, iss) result (dphase)

type (coord_struct) this_orb
type (internal_state_struct) iss
real(rp) rf_phase, dphase

!

this_orb = orbit
call track_this_lcavity (rf_phase, iss, this_orb, .false.)
dphase = twopi * (ele%value(rf_frequency$) / c_light) * (orbit%vec(5) / orbit%beta - this_orb%vec(5) / this_orb%beta)

end function dphase_end_minus_start

!---------------------------------------------------------------------------------------------
! contains
! Used with zbrent.

function phase_func (rf_phase, status) result (phase_err)

real(rp), intent(in) :: rf_phase
real(rp) phase_err
integer status

!

iss%dphase = dphase_end_minus_start(rf_phase, iss)
phase_err = rf_phase - (iss%phase0 + 0.5_rp * iss%dphase)
if (abs(phase_err) < phase_abs_tol) phase_err = 0

end function phase_func

!---------------------------------------------------------------------------------------------
! contains

! rf_phase is the phase used to calculate the change in energy.

subroutine track_this_lcavity (rf_phase, iss, orbit, make_matrix)

type (coord_struct) orbit
type (internal_state_struct) iss
real(rp) rf_phase, cdt, rel_p, dE, gradient_net, E_ratio, cos_phi, sin_phi, dcos_phi
real(rp) dbeta1_dE1, dbeta2_dE2, dalpha_dt1, dalpha_dE1, dcoef_dt1, dcoef_dE1, z21, z22
real(rp) dz_factor
logical, optional :: make_matrix

!

cos_phi = cos(rf_phase)
sin_phi = sin(rf_phase)
gradient_net = iss%gradient_max * cos_phi + gradient_shift_sr_wake(ele, param)

dE = gradient_net * iss%step_len
iss%E_end = iss%E_start + dE
if (iss%E_end <= mass_of(orbit%species)) then
  orbit%state = lost_pz$
  orbit%vec(6) = -1.01  ! Something less than -1
  if (present(mat6)) mat6 = 0
  return
endif

call convert_total_energy_to (iss%E_end, orbit%species, pc = iss%pc_end, beta = iss%beta_end)
E_ratio = iss%E_end / iss%E_start
sqrt_beta12 = sqrt(iss%beta_start/iss%beta_end)
mc2 = mass_of(orbit%species)
rel_p = 1 + orbit%vec(6)

! Body tracking longitudinal

cdt = iss%step_len * (iss%E_start + iss%E_end) / (iss%pc_end + iss%pc_start) 

if (logic_option(.false., make_matrix)) then
  om = twopi * ele%value(rf_frequency$) / c_light
  om_g = om * iss%gradient_max * iss%step_len

  dbeta1_dE1 = mc2**2 / (iss%pc_start * iss%E_start**2)
  dbeta2_dE2 = mc2**2 / (iss%pc_end * iss%E_end**2)

  ! Convert from (x, px, y, py, z, pz) to (x, x', y, y', c(t_ref-t), E) coords 
  mat6(2,:) = mat6(2,:) / rel_p - orbit%vec(2) * mat6(6,:) / rel_p**2
  mat6(4,:) = mat6(4,:) / rel_p - orbit%vec(4) * mat6(6,:) / rel_p**2

  m2(1,:) = [1/orbit%beta, -orbit%vec(5) * mc2**2 * orbit%p0c / (iss%pc_start**2 * iss%E_start)]
  m2(2,:) = [0.0_rp, orbit%p0c * orbit%beta]
  mat6(5:6,:) = matmul(m2, mat6(5:6,:))

  ! longitudinal track
  call mat_make_unit (kmat)
  kmat(6,5) = om_g * sin_phi
endif

! Convert to (x', y', c(t_ref-t), E) coords

orbit%vec(2) = orbit%vec(2) / rel_p    ! Convert to x'
orbit%vec(4) = orbit%vec(4) / rel_p    ! Convert to y'
orbit%vec(5) = orbit%vec(5) / iss%beta_start
orbit%vec(6) = rel_p * orbit%p0c / orbit%beta
orbit%t = orbit%t + cdt / c_light

! Body tracking. Note: Transverse kick only happens with standing wave cavities.

!--------------------------------------------------------------------
! Traveling wave
if (nint(ele%value(cavity_type$)) == traveling_wave$) then
  f_ave = (iss%pc_start + iss%pc_end) / (2 * iss%pc_end)
  dz_factor = (orbit%vec(2)**2 + orbit%vec(4)**2) * f_ave**2 * cdt / 2

  if (logic_option(.false., make_matrix)) then
    if (abs(dE) <  1d-4*(iss%pc_end+iss%pc_start)) then
      kmat(5,5) = 1 - iss%step_len * (-mc2**2 * kmat(6,5) / (2 * iss%pc_start**3) + mc2**2 * dE * kmat(6,5) * iss%E_start / iss%pc_start**5)
      kmat(5,6) = -iss%step_len * (-dbeta1_dE1 / iss%beta_start**2 + 2 * mc2**2 * dE / iss%pc_start**4 + &
                      (mc2 * dE)**2 / (2 * iss%pc_start**5) - 5 * (mc2 * dE)**2 / (2 * iss%pc_start**5))
    else
      kmat(5,5) = 1 - kmat(6,5) / (iss%beta_end * gradient_net) + kmat(6,5) * (iss%pc_end - iss%pc_start) / (gradient_net**2 * iss%step_len)
      kmat(5,6) = -1 / (iss%beta_end * gradient_net) + 1 / (iss%beta_start * gradient_net)
    endif

    kmat(1,2) = iss%step_len * f_ave
    kmat(1,5) = -orbit%vec(2) * kmat(5,6) * iss%step_len * iss%pc_start / (2 * iss%pc_end**2)
    kmat(1,6) =  orbit%vec(2) * (1 - kmat(6,6) * iss%pc_start / iss%pc_end) / (2 * iss%pc_end)

    kmat(2,2) = iss%pc_start / iss%pc_end
    kmat(2,5) = -orbit%vec(2) * kmat(5,6) * iss%pc_start / iss%pc_end**2
    kmat(2,6) =  orbit%vec(2) * (1 - kmat(6,6) * iss%pc_start / iss%pc_end) / iss%pc_end 

    kmat(3,4) = iss%step_len * f_ave
    kmat(3,5) = -orbit%vec(4) * kmat(5,6) * iss%step_len * iss%pc_start / (2 * iss%pc_end**2)
    kmat(3,6) =  orbit%vec(4) * (1 - kmat(6,6) * iss%pc_start / iss%pc_end) / (2 * iss%pc_end)

    kmat(4,4) = iss%pc_start / iss%pc_end
    kmat(4,5) = -orbit%vec(4) * kmat(5,6) * iss%pc_start / iss%pc_end**2
    kmat(4,6) =  orbit%vec(4) * (1 - kmat(6,6) * iss%pc_start / iss%pc_end) / iss%pc_end 

    kmat(5,2) = -orbit%vec(2) * f_ave**2 * cdt 
    kmat(5,4) = -orbit%vec(4) * f_ave**2 * cdt 

    mat6 = matmul(kmat, mat6)
  endif

  orbit%vec(1) = orbit%vec(1) + orbit%vec(2) * iss%step_len * f_ave
  orbit%vec(2) = orbit%vec(2) * iss%pc_start / iss%pc_end 
  orbit%vec(3) = orbit%vec(3) + orbit%vec(4) * iss%step_len * f_ave
  orbit%vec(4) = orbit%vec(4) * iss%pc_start / iss%pc_end 

  orbit%vec(5) = orbit%vec(5) - dz_factor
  orbit%t = orbit%t + dz_factor / c_light

!--------------------------------------------------------------------
! Standing wave
else
  sqrt_8 = 2 * sqrt_2
  voltage_max = iss%gradient_max * iss%step_len

  if (abs(voltage_max * cos_phi) < 1d-5 * iss%E_start) then
    f = voltage_max / iss%E_start
    alpha = f * (1 + f * cos_phi / 2)  / sqrt_8
    coef = iss%step_len * (1 - voltage_max * cos_phi / (2 * iss%E_start))
  else
    alpha = log(E_ratio) / (sqrt_8 * cos_phi)
    coef = sqrt_8 * iss%E_start * sin(alpha) / iss%gradient_max
  endif

  cos_a = cos(alpha)
  sin_a = sin(alpha)

  if (logic_option(.false., make_matrix)) then
    if (abs(voltage_max * cos_phi) < 1d-5 * iss%E_start) then
      dalpha_dt1 = f * f * om * sin_phi / (2 * sqrt_8) 
      dalpha_dE1 = -(voltage_max / iss%E_start**2 + voltage_max**2 * cos_phi / iss%E_start**3) / sqrt_8
      dcoef_dt1 = -iss%step_len * sin_phi * om_g / (2 * iss%E_start)
      dcoef_dE1 = iss%step_len * voltage_max * cos_phi / (2 * iss%E_start**2)
    else
      dalpha_dt1 = kmat(6,5) / (iss%E_end * sqrt_8 * cos_phi) - log(E_ratio) * om * sin_phi / (sqrt_8 * cos_phi**2)
      dalpha_dE1 = 1 / (iss%E_end * sqrt_8 * cos_phi) - 1 / (iss%E_start * sqrt_8 * cos_phi)
      dcoef_dt1 = sqrt_8 * iss%E_start * cos(alpha) * dalpha_dt1 / iss%gradient_max
      dcoef_dE1 = coef / iss%E_start + sqrt_8 * iss%E_start * cos(alpha) * dalpha_dE1 / iss%gradient_max
    endif

    z21 = -iss%gradient_max / (sqrt_8 * iss%E_end)
    z22 = iss%E_start / iss%E_end  

    c_min = cos_a - sqrt_2 * sin_a * cos_phi
    c_plu = cos_a + sqrt_2 * sin_a * cos_phi
    dc_min = -sin_a - sqrt_2 * cos_a * cos_phi 
    dc_plu = -sin_a + sqrt_2 * cos_a * cos_phi 

    cos_term = 1 + 2 * cos_phi**2
    dcos_phi = om * sin_phi

    kmat(1,1) =  c_min
    kmat(1,2) =  coef 
    kmat(2,1) =  z21 * cos_term * sin_a
    kmat(2,2) =  c_plu * z22

    kmat(1,5) = orbit%vec(1) * (dc_min * dalpha_dt1 - sqrt_2 * sin_a * dcos_phi) + & 
                 orbit%vec(2) * (dcoef_dt1)

    kmat(1,6) = orbit%vec(1) * dc_min * dalpha_dE1 + orbit%vec(2) * dcoef_dE1

    kmat(2,5) = orbit%vec(1) * z21 * (4 * cos_phi * dcos_phi * sin_a + cos_term * cos_a * dalpha_dt1) + &
                 orbit%vec(1) * (-kmat(2,1) * kmat(6,5) / (iss%E_end)) + &
                 orbit%vec(2) * z22 * (dc_plu * dalpha_dt1 + sqrt_2 * sin_a * dcos_phi - c_plu * kmat(6,5) / iss%E_end)

    kmat(2,6) = orbit%vec(1) * z21 * (cos_term * cos_a * dalpha_dE1) + &
                 orbit%vec(1) * (-kmat(2,1) / (iss%E_end)) + &
                 orbit%vec(2) * z22 * dc_plu * dalpha_dE1 + &
                 orbit%vec(2) * c_plu * (iss%E_end - iss%E_start) / iss%E_end**2

    kmat(3:4,3:4) = kmat(1:2,1:2)

    kmat(3,5) = orbit%vec(3) * (dc_min * dalpha_dt1 - sqrt_2 * iss%beta_start * sin_a * dcos_phi) + & 
                 orbit%vec(4) * (dcoef_dt1)

    kmat(3,6) = orbit%vec(3) * (dc_min * dalpha_dE1 - sqrt_2 * dbeta1_dE1 * sin_a * cos_phi) + &
                 orbit%vec(4) * (dcoef_dE1)

    kmat(4,5) = orbit%vec(3) * z21 * (4 * cos_phi * dcos_phi * sin_a + cos_term * cos_a * dalpha_dt1) + &
                 orbit%vec(3) * (-kmat(2,1) * kmat(6,5) / (iss%E_end)) + &
                 orbit%vec(4) * z22 * (dc_plu * dalpha_dt1 + sqrt_2 * sin_a * dcos_phi - c_plu * kmat(6,5) / iss%E_end)

    kmat(4,6) = orbit%vec(3) * z21 * (cos_term * cos_a * dalpha_dE1) + &
                 orbit%vec(3) * (-kmat(2,1) / iss%E_end) + &
                 orbit%vec(4) * z22 * dc_plu * dalpha_dE1 + &
                 orbit%vec(4) * c_plu * (iss%E_end - iss%E_start) / iss%E_end**2

    ! Correction to z for finite x', y'
    ! Note: Corrections to kmat(5,5) and kmat(5,6) are ignored since these are small (quadratic in the transvers coords).

    c_plu = sqrt_2 * cos_phi * cos_a + sin_a

    drp1_dr0  = -gradient_net / (2 * iss%E_start)
    drp1_drp0 = 1

    drp2_dr0  = (c_plu * z21)
    drp2_drp0 = (cos_a * z22)

    xp2 = drp2_dr0 * orbit%vec(1) + drp2_drp0 * orbit%vec(2)
    yp2 = drp2_dr0 * orbit%vec(3) + drp2_drp0 * orbit%vec(4)

    kmat(5,1) = -(orbit%vec(1) * (drp1_dr0**2 + drp1_dr0*drp2_dr0 + drp2_dr0**2) + &
                   orbit%vec(2) * (drp1_dr0 * drp1_drp0 + drp2_dr0 * drp2_drp0 + &
                                (drp1_dr0 * drp2_drp0 + drp1_drp0 * drp2_dr0) / 2)) * cdt / 3

    kmat(5,2) = -(orbit%vec(2) * (drp1_drp0**2 + drp1_drp0*drp2_drp0 + drp2_drp0**2) + &
                   orbit%vec(1) * (drp1_dr0 * drp1_drp0 + drp2_dr0 * drp2_drp0 + &
                                (drp1_dr0 * drp2_drp0 + drp1_drp0 * drp2_dr0) / 2)) * cdt / 3

    kmat(5,3) = -(orbit%vec(3) * (drp1_dr0**2 + drp1_dr0*drp2_dr0 + drp2_dr0**2) + &
                   orbit%vec(4) * (drp1_dr0 * drp1_drp0 + drp2_dr0 * drp2_drp0 + &
                                (drp1_dr0 * drp2_drp0 + drp1_drp0 * drp2_dr0) / 2)) * cdt / 3

    kmat(5,4) = -(orbit%vec(4) * (drp1_drp0**2 + drp1_drp0*drp2_drp0 + drp2_drp0**2) + &
                   orbit%vec(3) * (drp1_dr0 * drp1_drp0 + drp2_dr0 * drp2_drp0 + &
                                (drp1_dr0 * drp2_drp0 + drp1_drp0 * drp2_dr0) / 2)) * cdt / 3

    ! beta /= 1 corrections

    dsqrt_beta12 = -0.5_rp * sqrt_beta12 * dbeta2_dE2 * kmat(6,:) / iss%beta_end
    dsqrt_beta12(6) = dsqrt_beta12(6) + 0.5_rp * sqrt_beta12 * dbeta1_dE1 / iss%beta_start 

    kmat(1:4,1:4) = sqrt_beta12 * kmat(1:4,1:4)

    !

    mat6 = matmul(kmat, mat6)
  endif

  k1 = -gradient_net / (2 * iss%E_start)
  orbit%vec(2) = orbit%vec(2) + k1 * orbit%vec(1)    ! Entrance kick
  orbit%vec(4) = orbit%vec(4) + k1 * orbit%vec(3)    ! Entrance kick

  xp1 = orbit%vec(2)
  yp1 = orbit%vec(4)

  r_mat(1,1) =  cos_a
  r_mat(1,2) =  coef 
  r_mat(2,1) = -sin_a * iss%gradient_max / (sqrt_8 * iss%E_end)
  r_mat(2,2) =  cos_a * iss%E_start / iss%E_end

  orbit%vec(1:2) = sqrt_beta12 * matmul(r_mat, orbit%vec(1:2))   ! Modified R&S Eq 9.
  orbit%vec(3:4) = sqrt_beta12 * matmul(r_mat, orbit%vec(3:4))

  xp2 = orbit%vec(2)
  yp2 = orbit%vec(4)

  ! Correction of z for finite transverse velocity assumes a uniform change in slope.

  dz_factor = (xp1**2 + xp2**2 + xp1*xp2 + yp1**2 + yp2**2 + yp1*yp2) * cdt / 6
  orbit%vec(5) = orbit%vec(5) - dz_factor
  orbit%t = orbit%t + dz_factor / c_light

  !

  k2 = gradient_net / (2 * iss%E_end) 
  orbit%vec(2) = orbit%vec(2) + k2 * orbit%vec(1)         ! Exit kick
  orbit%vec(4) = orbit%vec(4) + k2 * orbit%vec(3)         ! Exit kick
endif

! 

orbit%vec(5) = orbit%vec(5) - (cdt - iss%cdt_ref)
orbit%vec(6) = (iss%pc_end - p0c_end) / p0c_end

! Convert back from (x', y', c(t_ref-t), E) coords

if (logic_option(.false., make_matrix)) then
  rel_p = iss%pc_end / p0c_end
  mat6(2,:) = rel_p * mat6(2,:) + orbit%vec(2) * mat6(6,:) / (p0c_end * iss%beta_end)
  mat6(4,:) = rel_p * mat6(4,:) + orbit%vec(4) * mat6(6,:) / (p0c_end * iss%beta_end)

  m2(1,:) = [iss%beta_end, orbit%vec(5) * mc2**2 / (iss%pc_end * iss%E_end**2)]
  m2(2,:) = [0.0_rp, 1 / (p0c_end * iss%beta_end)]

  mat6(5:6,:) = matmul(m2, mat6(5:6,:))
endif

orbit%vec(2) = orbit%vec(2) * (1 + orbit%vec(6))  ! Convert back to px
orbit%vec(4) = orbit%vec(4) * (1 + orbit%vec(6))  ! Convert back to py
orbit%vec(5) = orbit%vec(5) * iss%beta_end

orbit%beta = iss%beta_end

end subroutine track_this_lcavity

end subroutine track_a_lcavity
