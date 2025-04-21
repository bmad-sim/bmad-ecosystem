!+
! Subroutine rf_coupler_kick (ele, param, particle_at, phase, orbit, mat6, make_matrix)
! 
! Routine to add a RF cavity coupler kicks
!
! Input:
!   ele         -- ele_struct: Element being tracked through
!   param       -- lat_param_struct: branch parameters.
!   particle_at -- integer: first_track_edge$, or second_track_edge$.
!   phase       -- real(rp): phase of cavity
!   orbit       -- coord_struct: Position before kick.
!   mat6(6,6)   -- Real(rp), optional: Transfer matrix before the element.
!   make_matrix -- logical, optional: Propagate the transfer matrix? Default is false.
!
! Output:
!   orbit       -- coord_struct: Position after kick.
!   mat6(6,6)  -- real(rp), optional: Transfer matrix through the element.
!-

subroutine rf_coupler_kick (ele, param, particle_at, phase, orbit, mat6, make_matrix)

use bmad_routine_interface, dummy => rf_coupler_kick

implicit none

type (ele_struct) ele
type (lat_param_struct) param
type (coord_struct) orbit, old_orbit

real(rp), optional :: mat6(6,6)
real(rp) dp_coef, dp_x, dp_y, phase, ph, dE, mc2, p0c, pc, E, f, mc(6,6), coef, E_new, f2

integer particle_at, physical_end

logical, optional :: make_matrix

!

physical_end = physical_ele_end (particle_at, orbit, ele%orientation)
if (.not. at_this_ele_end (physical_end, nint(ele%value(coupler_at$)))) return
if (ele%value(coupler_strength$) == 0) return


ph = phase
if (ele%key == rfcavity$) ph = pi/2 - ph
ph = ph + twopi * ele%value(coupler_phase$)

mc2 = mass_of(orbit%species)
p0c = orbit%p0c
pc = p0c * (1 + orbit%vec(6))
E = pc / orbit%beta

f = twopi * ele%value(rf_frequency$) / c_light
dp_coef = e_accel_field(ele, gradient$) * ele%value(coupler_strength$)
dp_x = dp_coef * cos(twopi * ele%value(coupler_angle$))
dp_y = dp_coef * sin(twopi * ele%value(coupler_angle$))

if (nint(ele%value(coupler_at$)) == both_ends$) then
  dp_x = dp_x / 2
  dp_y = dp_y / 2
endif

! Track. 
! Since the transverse kick is time dependent there must be an associated energy kick
! to make things symplectic.


old_orbit = orbit

orbit%vec(2) = orbit%vec(2) + orbit%time_dir * dp_x * cos(ph) / orbit%p0c
orbit%vec(4) = orbit%vec(4) + orbit%time_dir * dp_y * cos(ph) / orbit%p0c

f2 = orbit%time_dir * sin(ph) * twopi * ele%value(rf_frequency$) / c_light
dE = (dp_x * orbit%vec(1) + dp_y * orbit%vec(3)) * f2
call apply_energy_kick (dE, orbit, [dp_x, dp_y] * f2)

! Matrix

if (logic_option(.false., make_matrix)) then
  call mat_make_unit (mc)

  mc(2,5) = dp_x * f * sin(ph) / (old_orbit%beta * p0c)
  mc(4,5) = dp_y * f * sin(ph) / (old_orbit%beta * p0c)

  mc(2,6) = -dp_x * f * sin(ph) * old_orbit%vec(5) * mc2**2 / (E * pc**2)
  mc(4,6) = -dp_y * f * sin(ph) * old_orbit%vec(5) * mc2**2 / (E * pc**2)

  coef = (dp_x * old_orbit%vec(1) + dp_y * old_orbit%vec(3)) * cos(ph) * f**2 
  mc(6,1) = dp_x * sin(ph) * f / (orbit%beta * p0c)
  mc(6,3) = dp_y * sin(ph) * f / (orbit%beta * p0c)
  mc(6,5) = -coef / (orbit%beta * old_orbit%beta * p0c) 
  mc(6,6) = old_orbit%beta/orbit%beta + coef * old_orbit%vec(5) * mc2**2 / (pc**2 * E * orbit%beta)

  f2 = old_orbit%vec(5) * mc2**2 / (pc * E**2 * p0c)
  E_new = p0c * (1 + orbit%vec(6)) / orbit%beta

  mc(5,1) = old_orbit%vec(5) * mc2**2 * p0c * mc(6,1) / (old_orbit%beta * E_new**3)
  mc(5,3) = old_orbit%vec(5) * mc2**2 * p0c * mc(6,3) / (old_orbit%beta * E_new**3)
  mc(5,5) = orbit%beta/old_orbit%beta + old_orbit%vec(5) * mc2**2 * p0c * mc(6,5) / (old_orbit%beta * E_new**3)
  mc(5,6) = old_orbit%vec(5) * mc2**2 * p0c * (mc(6,6) / (old_orbit%beta * E_new**3) - &
                                                     orbit%beta / (old_orbit%beta**2 * E**3))
  if (orbit%time_dir == -1) call mat_inverse(mc, mc)
  mat6 = matmul(mc, mat6)
endif

end subroutine rf_coupler_kick
