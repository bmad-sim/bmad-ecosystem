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
real(rp) length, s_now, s_end, kmat(6,6), phase, ds
real(rp) mc2, gradient_tot
real(rp) an(0:n_pole_maxx), bn(0:n_pole_maxx), an_elec(0:n_pole_maxx), bn_elec(0:n_pole_maxx)

integer ix_mag_max, ix_elec_max, ix_step_start, ix_step_end, n_steps, direction
integer ix_step

logical, optional :: make_matrix
logical make_mat

character(*), parameter :: r_name = 'track_a_lcavity'

! 

if (ele%value(rf_frequency$) == 0  .and. (ele%value(voltage$) /= 0 .or. ele%value(voltage_err$) /= 0)) then
  call out_io (s_error$, r_name, 'LCAVITY ELEMENT HAS ZERO RF_FREQUENCY: ' // ele%name)
  orbit%state = lost$
  return
endif

! The lord will have the step information.
! See the documentation for the rf_ele_struct for some details.

length = orbit%time_dir * ele%value(l$)
if (length == 0) return

make_mat = logic_option(.false., make_matrix)
lord => pointer_to_super_lord(ele)
mc2 = mass_of(orbit%species)
n_steps = ubound(lord%rf%steps, 1)

direction = orbit%time_dir * orbit%direction 
if (direction == 1) then
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

gradient_tot = direction * ele%value(gradient_tot$)

call multipole_ele_to_ab (ele, .false., ix_mag_max,  an,      bn,      magnetic$, include_kicks$)
call multipole_ele_to_ab (ele, .false., ix_elec_max, an_elec, bn_elec, electric$)

!

call offset_particle (ele, set$, orbit, mat6 = mat6, make_matrix = make_mat)

! Beginning Edge

phase = this_rf_phase(orbit, lord)
call rf_coupler_kick (ele, param, first_track_edge$, phase, orbit, mat6, make_mat)
if (fringe_here(ele, orbit, first_track_edge$)) call fringe_kick(orbit, lord, +1, phase, mc2, mat6, make_mat)

! Body

do ix_step = ix_step_start, ix_step_end, direction
  step => lord%rf%steps(ix_step)

  if (ix_step == ix_step_end) then
    ! Drift to end. The first and last steps have no drift section.
    if (ix_step == 0 .or. ix_step == n_steps) cycle
    ds = s_end - s_now
    call step_drift(orbit, ds, step, ele, param, mat6, make_mat)

  else
    ! Drift to edge of step and kick
    if (direction == 1) then
      ds = step%s - s_now
      call step_drift(orbit, ds, step, ele, param, mat6, make_mat)
      s_now = step%s
      call this_energy_kick(orbit, lord, step, direction, mat6, make_mat)
    else
      ds = lord%rf%steps(ix_step-1)%s - s_now
      call step_drift(orbit, ds, step, ele, param, mat6, make_mat)
      s_now = lord%rf%steps(ix_step-1)%s
      call this_energy_kick(orbit, lord, lord%rf%steps(ix_step-1), direction, mat6, make_mat)
    endif
  endif
enddo

! End Edge

phase = this_rf_phase(orbit, lord)
if (fringe_here(ele, orbit, second_track_edge$)) call fringe_kick(orbit, lord, -1, phase, mc2, mat6, make_mat)
call rf_coupler_kick (ele, param, second_track_edge$, phase, orbit, mat6, make_mat)

!

call offset_particle (ele, unset$, orbit, mat6 = mat6, make_matrix = make_mat)

!---------------------------------------------------------------------------------------
contains

subroutine step_drift(orbit, ds, step, ele, param, mat6, make_mat)

type (coord_struct) orbit
type (ele_struct) ele
type (lat_param_struct) param
type (rf_stair_step_struct) :: step

real(rp) ds, ks_rel
real(rp), optional :: mat6(6,6)
logical make_mat

!

ks_rel = ele%value(ks$) * ele%value(p0c$)

if (lord%value(ks$) == 0) then
  call track_a_drift(orbit, ds, mat6, make_mat, ele%orientation)
else
  call solenoid_track_and_mat (ele, ds, param, orbit, orbit, mat6, make_mat, ks_rel/orbit%p0c, step%p0c/step%E_tot0)
endif

end subroutine step_drift

!---------------------------------------------------------------------------------------
! contains

subroutine fringe_kick(orbit, lord, edge, phase, mc2, mat6, make_mat)

type (coord_struct) orbit
type (ele_struct) lord

real(rp) phase, mc2
real(rp), optional :: mat6(6,6)
real(rp) f, ff, dE, pc, ez_field, dez_dz_field, omega, pz_end

integer edge  ! +1 -> entrance end, -1 -> exit end.
logical make_mat

!

ff = edge * orbit%time_dir * charge_of(orbit%species) / (2.0_rp * charge_of(lord%ref_species))
f = ff / orbit%p0c
pc = orbit%p0c * (1 + orbit%vec(6))
ez_field = gradient_tot * cos(phase)
omega = twopi * ele%value(rf_frequency$) / c_light
dez_dz_field = gradient_tot * sin(phase) * omega
dE = -ff * 0.5_rp * dez_dz_field * (orbit%vec(1)**2 + orbit%vec(3)**2)
pz_end = orbit%vec(6) + dpc_given_dE(pc, mc2, dE) / orbit%p0c

call to_energy_coords(orbit, mc2, mat6, make_mat)

if (make_mat) then
  call mat_make_unit(kmat)
  kmat(2,1) = -f * ez_field
  kmat(2,5) = -f * dez_dz_field * orbit%vec(1)
  kmat(4,3) = -f * ez_field
  kmat(4,5) = -f * dez_dz_field * orbit%vec(3)
  kmat(6,1) = -ff * dez_dz_field * orbit%vec(1) / orbit%p0c
  kmat(6,3) = -ff * dez_dz_field * orbit%vec(3) / orbit%p0c
  kmat(6,5) = -ff * 0.5_rp * ez_field * (orbit%vec(1)**2 + orbit%vec(3)**2) * omega**2 / orbit%p0c
  mat6 = matmul(kmat, mat6)
endif

orbit%vec(2) = orbit%vec(2) - f * ez_field * orbit%vec(1)
orbit%vec(4) = orbit%vec(4) - f * ez_field * orbit%vec(3)
orbit%vec(6) = orbit%vec(6) + dE / orbit%p0c

call to_momentum_coords(orbit, pz_end, mc2, mat6, make_mat)

end subroutine fringe_kick

!---------------------------------------------------------------------------------------
! contains

subroutine this_energy_kick(orbit, lord, step, direction, mat6, make_mat)

type (coord_struct) orbit
type (ele_struct) lord
type (rf_stair_step_struct) :: step

real(rp), optional :: mat6(6,6)
real(rp) scale, t_ref, phase, rel_p, mc2, m2(2,2), dE, dE_amp, pz_end, p1c, dp0c
real(rp) pc_start, pc_end, om, r_pc, r2_pc, dp_dE

integer direction
logical make_mat

! Multipole half kicks

scale = 0.5_rp * step%scale

if (ix_mag_max > -1)  call ab_multipole_kicks (an,      bn,      ix_mag_max,  lord, orbit, magnetic$, rp8(orbit%time_dir)*scale,   mat6, make_mat)
if (ix_elec_max > -1) call ab_multipole_kicks (an_elec, bn_elec, ix_elec_max, lord, orbit, electric$, length*scale, mat6, make_mat)

!-------------------------------------------------
! Standing wave transverse half kick

if (nint(lord%value(cavity_type$)) == standing_wave$) call standing_wave_transverse_kick(orbit, lord, scale, mat6, make_mat)

!-------------------------------------------------
! Calc some stuff

if (direction == 1) then
  dp0c = step%dp0c
  p1c = step%p0c + dp0c
else
  dp0c = -step%dp0c
  p1c = step%p0c
endif

phase = this_rf_phase(orbit, lord)
dE_amp = direction * step%dE_amp 
dE = dE_amp * cos(phase)

rel_p = 1 + orbit%vec(6)
pc_start = rel_p * orbit%p0c
mc2 = mass_of(orbit%species)
pz_end = orbit%vec(6) + dpc_given_dE(orbit%p0c*rel_p, mc2, dE) / orbit%p0c
pc_end = (1 + pz_end) * orbit%p0c

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
call orbit_reference_energy_correction(orbit, dp0c, mat6, make_mat)

!-------------------------------------------------
! Standing wave transverse half kick

if (nint(lord%value(cavity_type$)) == standing_wave$) call standing_wave_transverse_kick(orbit, lord, scale, mat6, make_mat)

!-------------------------------------------------
! Multipole half kicks

if (ix_mag_max > -1)  call ab_multipole_kicks (an,      bn,      ix_mag_max,  lord, orbit, magnetic$, rp8(orbit%time_dir)*scale,   mat6, make_mat)
if (ix_elec_max > -1) call ab_multipole_kicks (an_elec, bn_elec, ix_elec_max, lord, orbit, electric$, length*scale, mat6, make_mat)

end subroutine this_energy_kick

!---------------------------------------------------------------------------------------
! contains

subroutine standing_wave_transverse_kick(orbit, lord, scale, mat6, make_mat)

type (coord_struct) orbit
type (ele_struct) lord

real(rp) scale, coef, kmat(6,6), rel_p
real(rp), optional :: mat6(6,6)

logical make_mat

!

rel_p = 1.0_rp + orbit%vec(6)
coef = scale * orbit%time_dir * lord%value(l$) * lord%value(gradient_tot$)**2 / (8.0_rp * orbit%p0c**2 * rel_p)

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

end subroutine standing_wave_transverse_kick

!---------------------------------------------------------------------------------------
! contains

function this_rf_phase(orbit, lord) result (phase)

type (coord_struct) orbit
type (ele_struct) lord

real(rp) phase

!

phase = twopi * (lord%value(phi0_err$) + lord%value(phi0$) + lord%value(phi0_multipass$) + &
           (particle_rf_time (orbit, lord, .false.) - rf_ref_time_offset(lord)) * lord%value(rf_frequency$))
if (bmad_com%absolute_time_tracking .and. lord%orientation*orbit%time_dir*orbit%direction == -1) then
  phase = phase - twopi * lord%value(rf_frequency$) * lord%value(delta_ref_time$)
endif
phase = modulo2(phase, pi)

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
