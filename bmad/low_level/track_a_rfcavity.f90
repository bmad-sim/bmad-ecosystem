!+
! Subroutine track_a_rfcavity (orbit, ele, param, mat6, make_matrix)
!
! Bmad_standard tracking through an rfcavity element.
!
! Input:
!   orbit       -- Coord_struct: Starting position.
!   ele         -- ele_struct: rfcavity element.
!   param       -- lat_param_struct: Lattice parameters.
!   make_matrix -- logical, optional: Propagate the transfer matrix? Default is false.
!
! Output:
!   orbit      -- coord_struct: End position.
!   mat6(6,6)  -- real(rp), optional: Transfer matrix through the element.
!-

subroutine track_a_rfcavity (orbit, ele, param, mat6, make_matrix)

use bmad_interface, except_dummy => track_a_rfcavity

implicit none

type (coord_struct) :: orbit
type (ele_struct), target :: ele
type (lat_param_struct) :: param

real(rp), optional :: mat6(6,6)
real(rp) beta_ref, dt_ref, voltage, phase0, phase, dE, z, charge_dir, length
real(rp) mc2, p0c, m2(2,2), t0, factor, pc, E, ff, new_beta, r_step, step_len
real(rp) an(0:n_pole_maxx), bn(0:n_pole_maxx), an_elec(0:n_pole_maxx), bn_elec(0:n_pole_maxx)

integer i, n_slice, n_step
integer ix_mag_max, ix_elec_max

logical, optional :: make_matrix
logical drifting

character(*), parameter :: r_name = 'track_a_rfcavity'

!

if (ele%value(rf_frequency$) == 0 .and. ele%value(voltage$) /= 0) then
  call out_io (s_error$, r_name, 'RFCAVITY ELEMENT HAS ZERO RF_FREQUENCY: ' // ele%name)
  orbit%state = lost$
  return
endif

length = ele%value(l$) * orbit%time_dir
n_step = 1
r_step = rp8(orbit%time_dir) / n_step
step_len = length / n_step

beta_ref = ele%value(p0c$) / ele%value(e_tot$)
n_slice = max(1, nint(abs(length) / ele%value(ds_step$))) 
dt_ref = length / (c_light * beta_ref)
charge_dir = rel_tracking_charge_to_mass(orbit, param%particle) * ele%orientation
mc2 = mass_of(orbit%species)
p0c = orbit%p0c

call multipole_ele_to_ab (ele, .false., ix_mag_max, an,      bn,      magnetic$, include_kicks$)
call multipole_ele_to_ab (ele, .false., ix_elec_max, an_elec, bn_elec, electric$)

call offset_particle (ele, set$, orbit, mat6 = mat6, make_matrix = make_matrix)

if (ix_mag_max > -1)  call ab_multipole_kicks (an,      bn,      ix_mag_max,  ele, orbit, magnetic$, 0.5_rp*r_step,   mat6, make_matrix)
if (ix_elec_max > -1) call ab_multipole_kicks (an_elec, bn_elec, ix_elec_max, ele, orbit, electric$, 0.5_rp*step_len, mat6, make_matrix)


! The cavity field is modeled as a standing wave antisymmetric wrt the center.
! Thus if the cavity is flipped (orientation = -1), the wave of interest, which is 
! always the accelerating wave, is the "backward" wave. And the phase of the backward 
! wave is different from the phase of the forward wave by a constant dt_ref * freq.
! Note: phi0_autoscale is not used here since bmad_standard tracking by design gives the correct tracking.
! In fact, using phi0_autoscale would be a mistake if, say, tracking_method = runge_kutta, mat6_calc_method = bmad_standard.

voltage = e_accel_field(ele, voltage$, .true.) * charge_dir

phase0 = twopi * (ele%value(phi0$) - &
                          (particle_rf_time (orbit, ele, .false., rf_freq = ele%value(rf_frequency$)) - &
                          rf_ref_time_offset(ele)) * ele%value(rf_frequency$))
if (.not. bmad_com%absolute_time_tracking) phase0 = phase0 + twopi * ele%value(phi0_multipass$)

if (bmad_com%absolute_time_tracking .and. ele%orientation*orbit%time_dir*orbit%direction == -1) then
  phase0 = phase0 - twopi * ele%value(rf_frequency$) * dt_ref
endif
phase = modulo2(phase0, pi)

call rf_coupler_kick (ele, param, first_track_edge$, phase, orbit, mat6, make_matrix)

! Track through slices.
! The phase of the accelerating wave traveling in the same direction as the particle is
! assumed to be traveling with a phase velocity the same speed as the reference velocity.

t0 = orbit%t

do i = 0, n_slice

  if (logic_option(.false., make_matrix)) then
    factor = voltage / n_slice
    if (i == 0 .or. i == n_slice) factor = factor / 2

    dE = factor * sin(phase)
    pc = (1 + orbit%vec(6)) * p0c 
    E = pc / orbit%beta
    call convert_total_energy_to (E + dE, orbit%species, beta = new_beta)
    ff = twopi * factor * ele%value(rf_frequency$) * cos(phase) / (p0c * new_beta * c_light)

    m2(2,1) = ff / orbit%beta
    m2(2,2) = orbit%beta / new_beta - ff * orbit%vec(5) *mc2**2 * p0c / (E * pc**2) 
    m2(1,1) = new_beta / orbit%beta + orbit%vec(5) * (mc2**2 * p0c * m2(2,1) / (E+dE)**3) / orbit%beta
    m2(1,2) = orbit%vec(5) * mc2**2 * p0c * (m2(2,2) / ((E+dE)**3 * orbit%beta) - new_beta / (pc**2 * E))
    if (orbit%time_dir == -1) call mat_inverse(m2, m2)

    mat6(5:6, :) = matmul(m2, mat6(5:6, :))
  endif

  dE = orbit%time_dir * voltage * sin(phase) / n_slice
  if (i == 0 .or. i == n_slice) dE = dE / 2

  call apply_energy_kick (dE, orbit, [0.0_rp, 0.0_rp])
  
  if (orbit%vec(6) == -1) then
    orbit%state = lost_pz$
    return
  endif

  if (i /= n_slice) then
    z = orbit%vec(5)
    call track_a_drift (orbit, length/n_slice, mat6, make_matrix)
    phase = phase + twopi * ele%value(rf_frequency$) * (orbit%vec(5) - z) / (c_light * orbit%beta)
  endif

enddo

! coupler kick, multipoles, back to lab coords.

call rf_coupler_kick (ele, param, second_track_edge$, phase, orbit, mat6, make_matrix)
orbit%phase(1) = modulo2(phase, 0.5_rp)


if (ix_elec_max > -1) call ab_multipole_kicks (an_elec, bn_elec, ix_elec_max, ele, orbit, electric$, 0.5_rp*step_len, mat6, make_matrix)
if (ix_mag_max > -1)  call ab_multipole_kicks (an,      bn,      ix_mag_max,  ele, orbit, magnetic$, 0.5_rp*r_step,   mat6, make_matrix)

call offset_particle (ele, unset$, orbit, mat6 = mat6, make_matrix = make_matrix)

end subroutine
