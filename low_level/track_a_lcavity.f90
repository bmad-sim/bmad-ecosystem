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
! Input:
!   orbit       -- Coord_struct: Starting position.
!   ele         -- ele_struct: Thick multipole element.
!   param       -- lat_param_struct: Lattice parameters.
!   make_matrix -- logical, optional: Propagate the transfer matrix? Default is false.
!
! Output:
!   orbit      -- coord_struct: End position.
!   mat6(6,6)  -- real(rp), optional: Transfer matrix through the element.
!-

subroutine track_a_lcavity (orbit, ele, param, mat6, make_matrix)

use track1_mod, except_dummy => track_a_lcavity

implicit none

type (coord_struct) :: orbit, start_orb
type (ele_struct), target :: ele
type (lat_param_struct) :: param

real(rp), optional :: mat6(6,6)
real(rp) length, pc_start, pc_end, gradient_ref, gradient_max, dz_factor, rel_p, coef, k2
real(rp) E_start_ref, E_end_ref, pc_start_ref, pc_end_ref, alpha, sin_a, cos_a, r_mat(2,2), volt_ref
real(rp) phase, cos_phi, gradient_net, e_start, e_end, e_ratio, voltage_max, dp_dg, sqrt_8, f, k1
real(rp) mc2, dE_start, dE_end, dE, beta_start, beta_end, beta_start_ref, beta_end_ref
real(rp) pxy2, xp1, xp2, yp1, yp2

logical, optional :: make_matrix

!

length = ele%value(l$)
if (length == 0) return

call offset_particle (ele, param, set$, orbit)

rel_p = 1 + orbit%vec(6)
E_start_ref  = ele%value(E_tot_start$)
E_end_ref    = ele%value(E_tot$)
gradient_ref = (E_end_ref - E_start_ref) / length
pc_start_ref = ele%value(p0c_start$)
pc_end_ref   = ele%value(p0c$)
beta_start_ref = pc_start_ref / E_start_ref
beta_end_ref   = pc_end_ref / E_end_ref

pc_start = pc_start_ref * rel_p
call convert_pc_to (pc_start, orbit%species, E_tot = E_start, beta = beta_start)

! The RF phase is defined with respect to the time at the beginning of the element.
! So if dealing with a slave element and absolute time tracking then need to correct.

phase = twopi * (ele%value(phi0_err$) + ele%value(phi0_autoscale$) + &
           ele%value(phi0$) + ele%value(phi0_multipass$) + &
           (particle_rf_time (orbit, ele, .false.) - rf_ref_time_offset(ele)) * ele%value(rf_frequency$))

gradient_max = e_accel_field(ele, gradient$)

cos_phi = cos(phase)
gradient_net = gradient_max * cos_phi + gradient_shift_sr_wake(ele, param)

dE = gradient_net * length
E_end = E_start + dE
if (E_end <= mass_of(orbit%species)) then
  orbit%state = lost_z_aperture$
  orbit%vec(6) = -1.01  ! Something less than -1
  return
endif

call convert_total_energy_to (E_end, orbit%species, pc = pc_end, beta = beta_end)
E_ratio = E_end / E_start
orbit%beta = beta_end
mc2 = mass_of(orbit%species)

! Coupler kick

call rf_coupler_kick (ele, param, first_track_edge$, phase, orbit)

! Body tracking longitudinal

orbit%vec(6) = (pc_end - pc_end_ref) / pc_end_ref 
orbit%p0c = pc_end_ref

if (abs(dE) <  1d-4*(pc_end+pc_start)) then
  dp_dg = length * (E_start / pc_start - mc2**2 * dE / (2 * pc_start**3) + (mc2 * dE)**2 * E_start / (2 * pc_start**5))
else
  dp_dg = (pc_end - pc_start) / gradient_net
endif

orbit%vec(5) = orbit%vec(5) * (beta_end / beta_start) - beta_end * (dp_dg - c_light * ele%value(delta_ref_time$))
orbit%t = orbit%t + dp_dg / c_light

! Body tracking transverse. Kick is only with standing wave cavities.

if (nint(ele%value(cavity_type$)) == traveling_wave$) then
  orbit%vec(2) = orbit%vec(2) / rel_p    ! Convert to x'
  orbit%vec(4) = orbit%vec(4) / rel_p    ! Convert to y'

  orbit%vec(1) = orbit%vec(1) + orbit%vec(2) * length
  orbit%vec(3) = orbit%vec(3) + orbit%vec(4) * length

  dz_factor = (orbit%vec(2)**2 + orbit%vec(4)**2) * dp_dg / 2
  orbit%vec(5) = orbit%vec(5) - beta_end * dz_factor
  orbit%t = orbit%t + dz_factor / c_light

  orbit%vec(2) = orbit%vec(2) * (1 + orbit%vec(6))  ! Convert back to px
  orbit%vec(4) = orbit%vec(4) * (1 + orbit%vec(6))  ! Convert back to py


else
  sqrt_8 = 2 * sqrt_2
  voltage_max = gradient_max * length

  if (abs(voltage_max * cos_phi) < 1d-5 * E_start) then
    f = voltage_max / E_start
    alpha = f * (1 + f * cos_phi / 2)  / sqrt_8
    coef = length * beta_start * (1 - voltage_max * cos_phi / (2 * E_start))
  else
    alpha = log(E_ratio) / (sqrt_8 * cos_phi)
    coef = sqrt_8 * pc_start * sin(alpha) / gradient_max
  endif

  cos_a = cos(alpha)
  sin_a = sin(alpha)

  r_mat(1,1) =  cos_a
  r_mat(1,2) =  coef 
  r_mat(2,1) = -sin_a * gradient_max / (sqrt_8 * pc_end)
  r_mat(2,2) =  cos_a * pc_start / pc_end

  orbit%vec(2) = orbit%vec(2) / rel_p    ! Convert to x'
  orbit%vec(4) = orbit%vec(4) / rel_p    ! Convert to y'

  k1 = -gradient_net / (2 * E_start)
  orbit%vec(2) = orbit%vec(2) + k1 * orbit%vec(1)    ! Entrance kick
  orbit%vec(4) = orbit%vec(4) + k1 * orbit%vec(3)    ! Entrance kick

  xp1 = orbit%vec(2)
  yp1 = orbit%vec(4)

  orbit%vec(1:2) = matmul(r_mat, orbit%vec(1:2))   ! R&S Eq 9.
  orbit%vec(3:4) = matmul(r_mat, orbit%vec(3:4))

  xp2 = orbit%vec(2)
  yp2 = orbit%vec(4)

  ! Correction of z for finite transverse velocity assumes a uniform change in slope.

  dz_factor = (xp1**2 + xp2**2 + xp1*xp2 + yp1**2 + yp2**2 + yp1*yp2) * dp_dg / 6
  orbit%vec(5) = orbit%vec(5) - beta_end * dz_factor
  orbit%t = orbit%t + dz_factor / c_light

  !

  k2 = gradient_net / (2 * E_end) 
  orbit%vec(2) = orbit%vec(2) + k2 * orbit%vec(1)         ! Exit kick
  orbit%vec(4) = orbit%vec(4) + k2 * orbit%vec(3)         ! Exit kick

  orbit%vec(2) = orbit%vec(2) * (1 + orbit%vec(6))  ! Convert back to px
  orbit%vec(4) = orbit%vec(4) * (1 + orbit%vec(6))  ! Convert back to py

endif

! Coupler kick

call rf_coupler_kick (ele, param, second_track_edge$, phase, orbit)

call offset_particle (ele, param, unset$, orbit)

end subroutine
