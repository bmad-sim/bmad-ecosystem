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
type (em_field_struct) field

real(rp), optional :: mat6(6,6)
real(rp) length, pc, E_ref_start, E_ref_end, s_now, s_end, kmat(6,6), phase, ds, f
real(rp) p0c_start, p0c_end
real(rp) an(0:n_pole_maxx), bn(0:n_pole_maxx), an_elec(0:n_pole_maxx), bn_elec(0:n_pole_maxx)

integer ix_mag_max, ix_elec_max, ix_step_start, ix_step_end, n_steps, direction
integer ix_step

logical, optional :: make_matrix

character(*), parameter :: r_name = 'track_a_lcavity'

! 

if (ele%value(rf_frequency$) == 0  .and. (ele%value(voltage$) /= 0 .or. ele%value(voltage_err$) /= 0)) then
  call out_io (s_error$, r_name, 'LCAVITY ELEMENT HAS ZERO RF_FREQUENCY: ' // ele%name)
  orbit%state = lost$
  return
endif

length = orbit%time_dir * ele%value(l$)
if (length == 0) return

lord => pointer_to_super_lord(ele)

!

call multipole_ele_to_ab (ele, .false., ix_mag_max,  an,      bn,      magnetic$, include_kicks$)
call multipole_ele_to_ab (ele, .false., ix_elec_max, an_elec, bn_elec, electric$)

call offset_particle (ele, set$, orbit, mat6 = mat6, make_matrix = make_matrix)

direction = orbit%time_dir * orbit%direction 
if (direction == 1) then
  E_ref_start = ele%value(E_tot_start$)
  E_ref_end   = ele%value(E_tot$)
  p0c_start   = ele%value(p0c_start$)
  p0c_end     = ele%value(p0c$)
  s_now = ele%s_start - lord%s_start
  s_end = ele%s       - lord%s_start
else
  E_ref_start = ele%value(E_tot$)
  E_ref_end   = ele%value(E_tot_start$)
  p0c_start   = ele%value(p0c_start$)
  p0c_end     = ele%value(p0c$)
  s_now = ele%s       - lord%s_start
  s_end = ele%s_start - lord%s_start
endif

! lord will have the step information.
! See the documentation for the rf_ele_struct for some details.

n_steps = ubound(ele%rf%steps, 1)
ix_step_start = rf_step_index(E_ref_start, s_now, lord)
ix_step_end   = rf_step_index(E_ref_end, s_end, lord)

! Beginning Edge
! Traveling_wave fringe (standing_wave fringe is built-in to the body formulas)

if (fringe_here(ele, orbit, first_track_edge$)) then
  phase = this_rf_phase(ix_step_start, orbit, lord)
  call rf_coupler_kick (ele, param, first_track_edge$, phase, orbit, mat6, make_matrix)

  if (nint(ele%value(cavity_type$)) == traveling_wave$) then
    ds = bmad_com%significant_length / 10  ! Make sure inside field region
    call em_field_calc (ele, param, ds, orbit, .true., field, logic_option(.false., make_matrix))
    f = charge_of(orbit%species) / (2 * p0c_end)
    pc = orbit%p0c * (1 + orbit%vec(6))

    if (logic_option(.false., make_matrix)) then
      call mat_make_unit(kmat)
      kmat(2,1) = -f * (field%dE(3,1) * orbit%vec(1) + field%E(3))
      kmat(2,3) = -f * field%dE(3,2) * orbit%vec(1) 
      kmat(2,5) = -f * field%dE(3,3) * orbit%vec(1) * orbit%beta
      kmat(2,6) =  f * field%E(3) * orbit%vec(1) * f / pc
      kmat(4,1) = -f * field%dE(3,1) * orbit%vec(3)
      kmat(4,3) = -f * (field%dE(3,2) * orbit%vec(3) + field%E(3))
      kmat(4,5) = -f * field%dE(3,3) * orbit%vec(3) * orbit%beta
      kmat(4,6) =  f * field%E(3) * orbit%vec(1) * f / pc
      mat6 = matmul(kmat, mat6)
    endif

    orbit%vec(2) = orbit%vec(2) - f * field%E(3) * orbit%vec(1)
    orbit%vec(4) = orbit%vec(4) - f * field%E(3) * orbit%vec(3)
  endif
endif

! Body

do ix_step = ix_step_start, ix_step_end, direction
  if (ix_step == ix_step_end) then
    ! Drift to end. The first and last steps have no drift section.
    if (ix_step == 0 .or. ix_step == n_steps) cycle
    ds = s_end - s_now
    call track_a_drift(orbit, ds, mat6, make_matrix, ele%orientation)

  else
    ! Drift to edge of step and kick
    if (direction == 1) then
      ds = lord%rf%steps(ix_step)%s - s_now
      call track_a_drift(orbit, ds, mat6, make_matrix, ele%orientation)
      s_now = lord%rf%steps(ix_step)%s
      call this_energy_kick(orbit, lord, lord%rf%steps(ix_step), direction)
    else
      ds = lord%rf%steps(ix_step-1)%s - s_now
      call track_a_drift(orbit, ds, mat6, make_matrix, ele%orientation)
      s_now = lord%rf%steps(ix_step-1)%s
      call this_energy_kick(orbit, lord, lord%rf%steps(ix_step-1), direction)
    endif
  endif
enddo

! End Edge

if (fringe_here(ele, orbit, second_track_edge$)) then
  if (nint(ele%value(cavity_type$)) == traveling_wave$) then
    ds = bmad_com%significant_length / 10  ! Make sure inside field region
    call em_field_calc (ele, param, length - ds, orbit, .true., field, logic_option(.false., make_matrix))
    f = -charge_of(orbit%species) / (2 * p0c_end)

    if (logic_option(.false., make_matrix)) then
      call mat_make_unit(kmat)
      pc = orbit%p0c * (1 + orbit%vec(6))
      kmat(2,1) = -f * (field%dE(3,1) * orbit%vec(1) + field%E(3))
      kmat(2,3) = -f * field%dE(3,2) * orbit%vec(1) 
      kmat(2,5) = -f * field%dE(3,3) * orbit%vec(1) * orbit%beta
      kmat(2,6) =  f * field%E(3) * orbit%vec(1) * f / pc
      kmat(4,1) = -f * field%dE(3,1) * orbit%vec(3)
      kmat(4,3) = -f * (field%dE(3,2) * orbit%vec(3) + field%E(3))
      kmat(4,5) = -f * field%dE(3,3) * orbit%vec(3) * orbit%beta
      kmat(4,6) =  f * field%E(3) * orbit%vec(1) * f / pc
      mat6 = matmul(kmat, mat6)
    endif

    orbit%vec(2) = orbit%vec(2) - f * field%E(3) * orbit%vec(1)
    orbit%vec(4) = orbit%vec(4) - f * field%E(3) * orbit%vec(3)
  endif

  ! Coupler kick
  phase = this_rf_phase(ix_step_end, orbit, lord)
  call rf_coupler_kick (ele, param, second_track_edge$, phase, orbit, mat6, make_matrix)
endif

call offset_particle (ele, unset$, orbit, mat6 = mat6, make_matrix = make_matrix)

!---------------------------------------------------------------------------------------
contains

!+
! Function rf_step_index(E_ref, s_rel, lord) result (ix_step)
!
! Routine to return the step index at a particular s-position or referece energy.
!
! E_ref is used in the computation instead of s_rel since s_rel is ambiguous if
! s_rel corresponds to a slice boundary.
!
! Input:
!   E_ref         -- real(rp): Reference energy of step
!   s_rel         -- real(rp): S-position relative to the beginning of the element
!   lord          -- real(rp): RF cavity.
!
! Output:
!   ix_step       -- integer: Corresponding index in the lord%rf%steps(:) array.
!-

function rf_step_index(E_ref, s_rel, lord) result (ix_step)

type (ele_struct) :: lord
real(rp) E_ref, dE, dE_rel, s_rel
integer ix_step
integer n_slice

character(*), parameter :: r_name = 'rf_step_index'

! The step index corres

if (.not. associated(lord%rf)) then
  call out_io(s_error$, r_name, 'MISSING RF STEP BOOKKEEPING PARAMETERS. PLEASE REPORT THIS! FOR ELEMENT ' // ele_full_name(lord))
  return
endif

n_slice = ubound(lord%rf%steps, 1) - 1  ! Number of slices
dE = lord%value(E_tot$) - lord%value(E_tot_start$)

! dE == 0 case. Must use s_rel in this case.

if (dE == 0) then
  if (s_rel == 0) then
    ix_step = 0
  elseif (s_rel == lord%value(l$)) then
    ix_step = n_slice + 1
  else
    ix_step = nint(n_slice * s_rel / lord%value(l$)) + 1
  endif
  return
endif

! dE /= 0 case

ix_step = nint(2.0_rp * (E_ref - lord%value(E_tot_start$)) / dE)
if (ix_step == 2*n_slice) then
  ix_step = n_slice + 1
elseif (ix_step > 0) then
  ix_step = (ix_step + 1) / 2
endif

! Sanity check

if (abs(E_ref - lord%rf%steps(ix_step)%E_tot0) > 1.0e-10 * (lord%value(E_tot$) + lord%value(E_tot_start$))) then
  call out_io(s_error$, r_name, 'RF STEP BOOKKEEPING FAILURE. PLEASE REPORT THIS!  FOR ELEMENT ' // ele_full_name(lord))
  return
endif

end function rf_step_index

!---------------------------------------------------------------------------------------
! contains

subroutine this_energy_kick(orbit, lord, step, direction)

type (coord_struct) orbit
type (ele_struct) lord
type (rf_stair_step_struct) :: step

real(rp) dE, scale, t_ref, phase, rel_p, pc, mc2, m2(2,2), E_tot0, E_tot1, E0, E1
integer direction

!

if (direction == 1) then
  E_tot0 = step%E_tot0
  E_tot1 = step%E_tot1
else
  E_tot0 = step%E_tot1
  E_tot1 = step%E_tot0
endif

!

scale = 0.5_rp * dE / (lord%value(E_tot$) - lord%value(E_tot_start$))
phase = this_rf_phase(ix_step, orbit, lord)
rel_p = 1 + orbit%vec(6)
mc2 = mass_of(orbit%species)

! Half the multipole kicks
! TODO: REF MOMENTUM/ENERGY CORRECTION!
if (ix_mag_max > -1)  call ab_multipole_kicks (an,      bn,      ix_mag_max,  lord, orbit, magnetic$, rp8(orbit%time_dir)*scale,   mat6, make_matrix)
if (ix_elec_max > -1) call ab_multipole_kicks (an_elec, bn_elec, ix_elec_max, lord, orbit, electric$, length*scale, mat6, make_matrix)

! Convert to (x', y', c(t_ref-t), E) coords

orbit%vec(2) = orbit%vec(2) / rel_p    ! Convert to x'
orbit%vec(4) = orbit%vec(4) / rel_p    ! Convert to y'
orbit%vec(5) = orbit%vec(5) / orbit%beta
orbit%vec(6) = rel_p * orbit%p0c / orbit%beta
E0 = orbit%vec(6)

! Slice Kick
! Traveling wave
if (nint(lord%value(cavity_type$)) == traveling_wave$) then

! Standing wave
else

endif

! Update to new energy

E1 = orbit%vec(6)
call convert_total_energy_to(E1, orbit%species, pc = pc)
orbit%beta = pc / E1


! Convert back from (x', y', c(t_ref-t), E) coords

if (logic_option(.false., make_matrix)) then
  pc = orbit%p0c * (1 + orbit%vec(6))
  rel_p = pc / orbit%p0c
  mat6(2,:) = rel_p * mat6(2,:) + orbit%vec(2) * mat6(6,:) / (orbit%p0c * orbit%beta)
  mat6(4,:) = rel_p * mat6(4,:) + orbit%vec(4) * mat6(6,:) / (orbit%p0c * orbit%beta)

  m2(1,:) = [orbit%beta, orbit%vec(5) * mc2**2 / (pc * E1**2)]
  m2(2,:) = [0.0_rp, 1 / (p0c_end * orbit%beta)]

  mat6(5:6,:) = matmul(m2, mat6(5:6,:))
endif

orbit%vec(2) = orbit%vec(2) * (1 + orbit%vec(6))  ! Convert back to px
orbit%vec(4) = orbit%vec(4) * (1 + orbit%vec(6))  ! Convert back to py
orbit%vec(5) = orbit%vec(5) * orbit%beta

! Half the multipole kicks

if (ix_mag_max > -1)  call ab_multipole_kicks (an,      bn,      ix_mag_max,  lord, orbit, magnetic$, rp8(orbit%time_dir)*scale,   mat6, make_matrix)
if (ix_elec_max > -1) call ab_multipole_kicks (an_elec, bn_elec, ix_elec_max, lord, orbit, electric$, length*scale, mat6, make_matrix)

end subroutine this_energy_kick

!---------------------------------------------------------------------------------------
! contains

function this_rf_phase(ix_step, orbit, lord) result (phase)

type (coord_struct) orbit
type (ele_struct) lord

real(rp) phase
integer ix_step

!

phase = twopi * (lord%value(phi0_err$) + lord%value(phi0$) + lord%value(phi0_multipass$) + &
           (particle_rf_time (orbit, lord, .false.) - rf_ref_time_offset(lord)) * lord%value(rf_frequency$))
if (bmad_com%absolute_time_tracking .and. lord%orientation*orbit%time_dir*orbit%direction == -1) then
  phase = phase - twopi * lord%value(rf_frequency$) * lord%value(delta_ref_time$)
endif
phase = modulo2(phase, pi)

end function this_rf_phase

end subroutine track_a_lcavity
