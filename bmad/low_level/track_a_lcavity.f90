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

implicit none

type (coord_struct) :: orbit
type (ele_struct), target :: ele
type (ele_struct), pointer :: lord
type (lat_param_struct) :: param
type (rf_stair_step_struct), pointer :: step

real(rp), optional :: mat6(6,6)
real(rp) s_now, s_end, kmat(6,6), phase, ds, mc2
real(rp) an(0:n_pole_maxx), bn(0:n_pole_maxx), an_elec(0:n_pole_maxx), bn_elec(0:n_pole_maxx), l_active

integer ix_mag_max, ix_elec_max, ix_step_start, ix_step_end, n_steps, s_dir
integer ix_step, body_dir

logical, optional :: make_matrix
logical make_mat, track_spin

character(*), parameter :: r_name = 'track_a_lcavity'

! 

if (ele%value(rf_frequency$) == 0  .and. (ele%value(voltage$) /= 0 .or. ele%value(voltage_err$) /= 0)) then
  call out_io (s_error$, r_name, 'LCAVITY ELEMENT HAS ZERO RF_FREQUENCY: ' // ele%name)
  orbit%state = lost$
  return
endif

if (ele%value(l$) == 0) return

! The lord will have the step information.
! See the documentation for the rf_ele_struct for some details.

l_active = orbit%time_dir * ele%value(l_active$)
track_spin = (bmad_com%spin_tracking_on .and. ele%spin_tracking_method == tracking$)

make_mat = logic_option(.false., make_matrix)
lord => pointer_to_super_lord(ele)
mc2 = mass_of(orbit%species)
n_steps = ubound(lord%rf%steps, 1) - 1

s_dir = orbit%time_dir * orbit%direction        ! Longitudinal propagation direction 
body_dir = orbit%direction * ele%orientation    ! Forward time direction of travel with respect to body coordinates

if (s_dir == 1) then
  s_now = ele%s_start - lord%s_start
  s_end = ele%s       - lord%s_start
  ix_step_start = ele_rf_step_index(ele%value(E_tot_start$), s_now, lord)
  ix_step_end   = ele_rf_step_index(ele%value(E_tot$), s_end, lord)
else
  s_now = ele%s       - lord%s_start
  s_end = ele%s_start - lord%s_start
  ix_step_start = ele_rf_step_index(ele%value(E_tot$), s_now, lord)
  ix_step_end   = ele_rf_step_index(ele%value(E_tot_start$), s_end, lord)
endif

call multipole_ele_to_ab (ele, .false., ix_mag_max,  an,      bn,      magnetic$, include_kicks$)
call multipole_ele_to_ab (ele, .false., ix_elec_max, an_elec, bn_elec, electric$)

!

call offset_particle (ele, set$, orbit, mat6 = mat6, make_matrix = make_mat)

! Forward tracking

if (s_dir == 1) then
  do ix_step = ix_step_start, ix_step_end
    step => lord%rf%steps(ix_step)
    ! Drift to edge of step
    if (ix_step == ix_step_end) then
      ds = s_end - s_now
      s_now = s_end
    else
      ds = step%s - s_now
      s_now = step%s
    endif
    call step_drift(orbit, ds, step, lord, param, mat6, make_mat)
    ! Entrence fringe?
    if (ix_step == 0 .and. s_now == step%s .and. fringe_here(lord, orbit, first_track_edge$)) then
      phase = this_rf_phase(orbit, ele, lord, lord%rf%steps(0))
      call rf_coupler_kick (ele, param, first_track_edge$, phase, orbit, mat6, make_mat)
      call fringe_kick(orbit, lord, +1, phase, body_dir, mc2, mat6, make_mat)
    endif
    ! Stair step kick
    if (ix_step /= ix_step_end) call this_energy_kick(orbit, lord, step, body_dir, mat6, make_mat)
    ! Exit fringe?
    if (ix_step == n_steps .and. s_now == step%s .and. fringe_here(lord, orbit, second_track_edge$)) then
      phase = this_rf_phase(orbit, ele, lord, lord%rf%steps(n_steps))
      call fringe_kick(orbit, lord, -1, phase, body_dir, mc2, mat6, make_mat)
      call rf_coupler_kick (ele, param, second_track_edge$, phase, orbit, mat6, make_mat)
    endif
  enddo
endif

! Reverse tracking

if (s_dir == -1) then
  do ix_step = ix_step_start, ix_step_end, -1
    step => lord%rf%steps(ix_step)
    ! Drift to edge of step
    if (ix_step == ix_step_end) then
      ds = s_end - s_now
      s_now = s_end
    else
      ds = step%s0 - s_now
      s_now = step%s0
    endif
    call step_drift(orbit, ds, step, lord, param, mat6, make_mat)
    ! Entrence fringe?
    if (ix_step == n_steps+1 .and. s_now == step%s0 .and. fringe_here(ele, orbit, first_track_edge$)) then
      phase = this_rf_phase(orbit, ele, lord, lord%rf%steps(n_steps))
      call rf_coupler_kick (ele, param, first_track_edge$, phase, orbit, mat6, make_mat)
      call fringe_kick(orbit, lord, +1, phase, body_dir, mc2, mat6, make_mat)
    endif
    ! Stair step kick
    if (ix_step /= ix_step_end) call this_energy_kick(orbit, lord, lord%rf%steps(ix_step-1), body_dir, mat6, make_mat)
    ! Exit fringe?
    if (ix_step == 1 .and. s_now == step%s0 .and. fringe_here(ele, orbit, second_track_edge$)) then
      phase = this_rf_phase(orbit, ele, lord, lord%rf%steps(0))
      call fringe_kick(orbit, lord, -1, phase, body_dir, mc2, mat6, make_mat)
      call rf_coupler_kick (ele, param, second_track_edge$, phase, orbit, mat6, make_mat)
    endif
  enddo
endif

!

call offset_particle (ele, unset$, orbit, mat6 = mat6, make_matrix = make_mat)

!---------------------------------------------------------------------------------------
contains

subroutine step_drift(orbit, ds, step, lord, param, mat6, make_mat)

type (coord_struct) orbit
type (ele_struct) lord
type (lat_param_struct) param
type (rf_stair_step_struct) :: step
type (em_field_struct) field

real(rp) ds, ks_rel, s_omega(3)
real(rp), optional :: mat6(6,6)
logical make_mat

!

if (lord%value(ks$) == 0) then
  call track_a_drift(orbit, ds, mat6, make_mat, lord%orientation)
else
  if (track_spin) then
    field = em_field_struct()
    field%b(3) = 0.5_rp * ds * lord%value(bs_field$)
    s_omega = spin_omega(field, orbit, orbit%direction*lord%orientation)
    call rotate_spin(s_omega, orbit%spin)
    call solenoid_track_and_mat (lord, ds, param, orbit, orbit, mat6, make_mat)
    s_omega = spin_omega(field, orbit, orbit%direction*lord%orientation)
    call rotate_spin(s_omega, orbit%spin)
  else
    call solenoid_track_and_mat (lord, ds, param, orbit, orbit, mat6, make_mat)
  endif
endif

end subroutine step_drift

!---------------------------------------------------------------------------------------
! contains

subroutine fringe_kick(orbit, lord, edge, phase, body_dir, mc2, mat6, make_mat)

type (coord_struct) orbit
type (ele_struct) lord
type (em_field_struct) field

real(rp) phase, mc2
real(rp), optional :: mat6(6,6)
real(rp) f, ff, dE, pc, ez_field, dez_dz_field, rf_omega, pz_end, s_omega(3), gradient_tot

integer edge  ! +1 -> entrance end, -1 -> exit end.
integer body_dir
logical make_mat

! The fringe kick only exists if there is a wave traveling in the same direction as the particle.

if (nint(lord%value(cavity_type$)) == traveling_wave$ .and. body_dir == -1) return

! Init

gradient_tot = body_dir * orbit%time_dir * lord%value(gradient_tot$) * lord%value(field_autoscale$) * lord%value(l$) / lord%value(l_active$)
ff = edge * orbit%time_dir * charge_of(orbit%species) / (2.0_rp * charge_of(lord%ref_species))
f = ff / orbit%p0c
pc = orbit%p0c * (1 + orbit%vec(6))
ez_field = gradient_tot * cos(phase)
rf_omega = twopi * ele%value(rf_frequency$) / c_light
dez_dz_field = gradient_tot * sin(phase) * rf_omega
dE = -ff * 0.5_rp * dez_dz_field * (orbit%vec(1)**2 + orbit%vec(3)**2)
pz_end = orbit%vec(6) + dpc_given_dE(pc, mc2, dE) / orbit%p0c

! Spin

if (track_spin) then
  field = em_field_struct()
  field%E(1:2) = (-0.5_rp * ff * ez_field / charge_of(orbit%species)) * [orbit%vec(1), orbit%vec(3)]
  s_omega = spin_omega(field, orbit, orbit%direction*lord%orientation)
  call rotate_spin(s_omega, orbit%spin)
endif

! Kick

call to_energy_coords(orbit, mc2, mat6, make_mat)

if (make_mat) then
  call mat_make_unit(kmat)
  kmat(2,1) = -f * ez_field
  kmat(2,5) = -f * dez_dz_field * orbit%vec(1)
  kmat(4,3) = -f * ez_field
  kmat(4,5) = -f * dez_dz_field * orbit%vec(3)
  kmat(6,1) = -ff * dez_dz_field * orbit%vec(1) / orbit%p0c
  kmat(6,3) = -ff * dez_dz_field * orbit%vec(3) / orbit%p0c
  kmat(6,5) = -ff * 0.5_rp * ez_field * (orbit%vec(1)**2 + orbit%vec(3)**2) * rf_omega**2 / orbit%p0c
  mat6 = matmul(kmat, mat6)
endif

orbit%vec(2) = orbit%vec(2) - f * ez_field * orbit%vec(1)
orbit%vec(4) = orbit%vec(4) - f * ez_field * orbit%vec(3)
orbit%vec(6) = orbit%vec(6) + dE / orbit%p0c

call to_momentum_coords(orbit, pz_end, mc2, mat6, make_mat)

! Spin

if (track_spin) then
  s_omega = spin_omega(field, orbit, orbit%direction*lord%orientation)
  call rotate_spin(s_omega, orbit%spin)
endif

end subroutine fringe_kick

!---------------------------------------------------------------------------------------
! contains

subroutine this_energy_kick(orbit, lord, step, body_dir, mat6, make_mat)

type (coord_struct) orbit
type (ele_struct) lord
type (rf_stair_step_struct) :: step
type (em_field_struct) field

real(rp), optional :: mat6(6,6)
real(rp) scale, t_ref, phase, rel_p, mc2, m2(2,2), dE, dE_amp, pz_end, p1c
real(rp) pc_start, pc_end, om, r_pc, r2_pc, dp_dE, s_omega(3)

integer body_dir
logical make_mat

! Multipole half kicks

scale = 0.5_rp * step%scale

if (ix_mag_max > -1)  call ab_multipole_kicks (an,      bn,      ix_mag_max,  lord, orbit, magnetic$, rp8(orbit%time_dir)*scale,   mat6, make_mat)
if (ix_elec_max > -1) call ab_multipole_kicks (an_elec, bn_elec, ix_elec_max, lord, orbit, electric$, lord%value(l$)*scale, mat6, make_mat)

!-------------------------------------------------
! Standing wave transverse half kick

call pondermotive_transverse_kick(orbit, lord, scale, body_dir, mat6, make_mat)

!-------------------------------------------------
! Calc some stuff

if (s_dir == 1) then
  p1c = step%p1c
else
  p1c = step%p0c
endif

phase = this_rf_phase(orbit, ele, lord, step)

dE_amp = s_dir * (lord%value(voltage$) + lord%value(voltage_err$)) * step%scale * lord%value(field_autoscale$) ! Amplitude of RF field
dE = dE_amp * cos(phase)

rel_p = 1 + orbit%vec(6)
pc_start = rel_p * orbit%p0c
mc2 = mass_of(orbit%species)
pz_end = orbit%vec(6) + dpc_given_dE(orbit%p0c*rel_p, mc2, dE) / orbit%p0c
pc_end = (1 + pz_end) * orbit%p0c

!-------------------------------------------------
! Spin

if (track_spin) then
  field = em_field_struct()
  field%e(3) = dE / charge_of(orbit%species)
  s_omega = spin_omega(field, orbit, orbit%direction*lord%orientation)
  call rotate_spin(s_omega, orbit%spin)
endif

!-------------------------------------------------
! Convert to energy coords

call to_energy_coords(orbit, mc2, mat6, make_mat)

!-------------------------------------------------
! Energy kick

orbit%vec(6) = orbit%vec(6) + dE / orbit%p0c

if (make_mat) then
  om = twopi * ele%value(rf_frequency$) / c_light
  mat6(6,:) = mat6(6,:) + (om * dE_amp * sin(phase) / orbit%p0c) * mat6(5,:)
endif

!-------------------------------------------------
! Convert to momentum coords

call to_momentum_coords(orbit, pz_end, mc2, mat6, make_mat)
call orbit_reference_energy_correction(orbit, p1c, mat6, make_mat)

!-------------------------------------------------
! Spin

if (track_spin) then
  s_omega = spin_omega(field, orbit, orbit%direction*lord%orientation)
  call rotate_spin(s_omega, orbit%spin)
endif

!-------------------------------------------------
! Standing wave transverse half kick

call pondermotive_transverse_kick(orbit, lord, scale, body_dir, mat6, make_mat)

!-------------------------------------------------
! Multipole half kicks

if (ix_mag_max > -1)  call ab_multipole_kicks (an,      bn,      ix_mag_max,  lord, orbit, magnetic$, rp8(orbit%time_dir)*scale,   mat6, make_mat)
if (ix_elec_max > -1) call ab_multipole_kicks (an_elec, bn_elec, ix_elec_max, lord, orbit, electric$, lord%value(l$)*scale, mat6, make_mat)

end subroutine this_energy_kick

!---------------------------------------------------------------------------------------
! contains

subroutine pondermotive_transverse_kick(orbit, lord, scale, body_dir, mat6, make_mat)

type (coord_struct) orbit
type (ele_struct) lord

real(rp) scale, coef, kmat(6,6), rel_p
real(rp), optional :: mat6(6,6)

integer body_dir
logical make_mat

! The pondermotive force only occurs if there is a EM wave in the opposite direction from the direction of travel.

if (nint(lord%value(cavity_type$)) == traveling_wave$ .and. body_dir == 1) return

!

rel_p = 1.0_rp + orbit%vec(6)
coef = (lord%value(field_autoscale$)*lord%value(gradient_tot$)*lord%value(l$)/lord%value(l_active$))**2 * &
                                      scale * orbit%time_dir * lord%value(l$) / (8.0_rp * orbit%p0c**2 * rel_p)

if (make_mat) then
  call mat_make_unit(kmat)
  kmat(2,1) = -coef
  kmat(2,6) =  coef * orbit%vec(1) / rel_p
  kmat(4,3) = -coef
  kmat(4,6) =  coef * orbit%vec(3) / rel_p
  kmat(5,1) = -coef * orbit%vec(1) / rel_p
  kmat(5,3) = -coef * orbit%vec(3) / rel_p
  kmat(5,6) =  coef * (orbit%vec(1)**2 + orbit%vec(3)**2) / rel_p**2
  mat6 = matmul(kmat, mat6)
endif

orbit%vec(2) = orbit%vec(2) - coef * orbit%vec(1)
orbit%vec(4) = orbit%vec(4) - coef * orbit%vec(3)
orbit%vec(5) = orbit%vec(5) - 0.5_rp * coef * (orbit%vec(1)**2 + orbit%vec(3)**2) / rel_p

end subroutine pondermotive_transverse_kick

!---------------------------------------------------------------------------------------
! contains

function this_rf_phase(orbit, ele, lord, step) result (phase)

type (coord_struct) orbit
type (ele_struct), target :: ele, lord
type (ele_struct), pointer :: mlord
type (rf_stair_step_struct) step
type (ele_pointer_struct), allocatable :: chain(:)
real(rp) phase, particle_time
integer ix_pass, n_links

! Lord will be something like the super_lord to a super_slave. 
! The lord can be a multipass_slave but will not be a multipass_lord.

mlord => lord
if (lord%slave_status == multipass_slave$) mlord => pointer_to_lord(lord, 1)

if (absolute_time_tracking(lord)) then
  particle_time = modulo2(orbit%t, 0.5_qp / lord%value(rf_frequency$))
  if (bmad_com%absolute_time_ref_shift) then
    if (lord%slave_status == multipass_slave$) then
      call multipass_chain(lord, ix_pass, n_links, chain)
      particle_time = particle_time - chain(1)%ele%value(ref_time_start$)
    else
      particle_time = particle_time - lord%value(ref_time_start$)
    endif
  endif

  particle_time = particle_time - mlord%rf%steps(step%ix_step)%time

else  ! Relative time tracking
  particle_time = particle_rf_time (orbit, lord, .false.) + step%time - mlord%rf%steps(step%ix_step)%time
endif

!

phase = twopi * (lord%value(phi0_err$) + lord%value(phi0$) + particle_time * lord%value(rf_frequency$))
if (.not. bmad_com%absolute_time_tracking) phase = phase + twopi * lord%value(phi0_multipass$)

if (bmad_com%absolute_time_tracking .and. lord%orientation*orbit%time_dir*orbit%direction == -1) then
  phase = phase - twopi * lord%value(rf_frequency$) * lord%value(delta_ref_time$)
endif

phase = modulo2(phase, pi)
orbit%phase(1) = phase/twopi

end function this_rf_phase

!---------------------------------------------------------------------------------------
! contains

subroutine to_energy_coords(orbit, mc2, mat6, make_mat)

type (coord_struct) orbit
real(rp), optional :: mat6(6,6)
real(rp) mc2, m2(2,2), pc
logical make_mat

! Convert from (z, pz) to (c(t0-t), E) coords 

if (make_mat) then
  pc = (1 + orbit%vec(6)) * orbit%p0c
  m2(1,:) = [1/orbit%beta, -orbit%vec(5) * mc2**2 * orbit%p0c * orbit%beta / pc**3]
  m2(2,:) = [0.0_rp, orbit%beta]
  mat6(5:6,:) = matmul(m2, mat6(5:6,:))
endif

orbit%vec(5) = orbit%vec(5) / orbit%beta
orbit%vec(6) = (1 + orbit%vec(6)) / orbit%beta

end subroutine to_energy_coords

!---------------------------------------------------------------------------------------
! contains

subroutine to_momentum_coords(orbit, pz, mc2, mat6, make_mat)

type (coord_struct) orbit
real(rp), optional :: mat6(6,6)
real(rp) pz, mc2, m2(2,2), pc
logical make_mat

! Convert back from (c(t0-t), E) coords to (z, pz)
! Note: pz is passed in as an argument to eliminate round-off error if pz were
! to be calculated in this routine.

pc = (1 + pz) * orbit%p0c
orbit%beta = pc / (orbit%p0c * orbit%vec(6))

if (make_mat) then
  m2(1,:) = [orbit%beta, orbit%vec(5) * mc2**2 * orbit%p0c * orbit%beta**2 / pc**3]
  m2(2,:) = [0.0_rp, 1.0_rp / orbit%beta]
  mat6(5:6,:) = matmul(m2, mat6(5:6,:))
endif

orbit%vec(5) = orbit%vec(5) * orbit%beta
orbit%vec(6) = pz

end subroutine to_momentum_coords

end subroutine track_a_lcavity
