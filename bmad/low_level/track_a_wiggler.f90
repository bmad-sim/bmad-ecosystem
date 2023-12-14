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

use bmad_interface, except_dummy => track_a_wiggler

implicit none

type (coord_struct) :: orbit
type (ele_struct), target :: ele
type (ele_struct), pointer :: field_ele
type (lat_param_struct) :: param

real(rp), optional :: mat6(6,6)
real(rp) mat2(2,2), z_start, beta_ref, p_factor, k1x, k1y, k1yy, k1xx, k3l, length, ky2, kz, rel_p, t_start
real(rp) an(0:n_pole_maxx), bn(0:n_pole_maxx), an_elec(0:n_pole_maxx), bn_elec(0:n_pole_maxx)
real(rp) m43, m46, m52, m54, m56, mc2_rel, kmat(6,6), factor, gamma0
real(rp) dz_x(3), dz_y(3), ddz_x(3), ddz_y(3), r_step, step_len

integer i, ix_mag_max, ix_elec_max, n_step

logical, optional :: make_matrix

! For planar wigglers:
! In the horizontal plane the tracking looks like a drift.
! The tracking in the vertical plane is:
!   1) 1/2 the octupole kick at the entrance face.
!   2) Track as a quadrupole through the body
!   3) 1/2 the octupole kick at the exit face.

! For helical wigglers use the vertical plane tracking for the horizontal plane as well.

call multipole_ele_to_ab (ele, .false., ix_mag_max, an,      bn,      magnetic$, include_kicks$)
call multipole_ele_to_ab (ele, .false., ix_elec_max, an_elec, bn_elec, electric$)

field_ele => pointer_to_field_ele(ele, 1)

z_start = orbit%vec(5)
t_start = orbit%t

! ds_step is set for PTC tracking which needs a lot more steps. 
! Need only one step if there are no multipoles.

n_step = max(nint(ele%value(l$) / ele%value(ds_step$)), 1)
if (ix_mag_max  < 0 .and. ix_elec_max < 0) n_step = 1

r_step = real(orbit%time_dir, rp) / n_step
step_len = ele%value(l$) * r_step

length = orbit%time_dir * ele%value(l$)
mc2_rel = mass_of(orbit%species) / orbit%p0c

if (ele%value(l_period$) == 0) then
  kz = 1d100    ! Something large
  ky2 = 0
else
  kz = twopi / ele%value(l_period$)
  ky2 = kz**2 + ele%value(kx$)**2
endif

factor = abs(rel_tracking_charge_to_mass(orbit, param%particle)) * 0.5 * (c_light * ele%value(b_max$) / ele%value(p0c$))**2
if (field_ele%field_calc == helical_model$) then
  k1x = -factor
  k1y = -factor
else
  k1x =  factor * (ele%value(kx$) / kz)**2
  k1y = -factor * ky2 / kz**2
endif

!

call offset_particle (ele, set$, orbit, mat6 = mat6, make_matrix = make_matrix)

if (ix_mag_max > -1)  call ab_multipole_kicks (an,      bn,      ix_mag_max,  ele, orbit, magnetic$, r_step/2,   mat6, make_matrix)
if (ix_elec_max > -1) call ab_multipole_kicks (an_elec, bn_elec, ix_elec_max, ele, orbit, electric$, step_len/2, mat6, make_matrix)

! Body

do i = 1, n_step

  ! Notice that the averageed quadrupole and octupole kicks are independent of the sign of the particle charge and independent 
  ! of the longitudinal direction of travel of the particle!

  beta_ref = ele%value(p0c$) / ele%value(e_tot$)
  rel_p = 1 + orbit%vec(6)
  k1yy = k1y / rel_p**2

  ! 1/2 of the octupole kick.

  k3l = 2 * step_len * k1yy
  if (i == 1) k3l = k3l / 2

  if (logic_option(.false., make_matrix)) then
    m43 = k3l * rel_p * kz**2 * orbit%vec(3)**2
    m46 = -k3l * kz**2 * orbit%vec(3)**3 / 3
    mat6(4,:) = mat6(4,:) + m43 * mat6(3,:) + m46 * mat6(6,:)
    if (field_ele%field_calc == helical_model$) then
      mat6(2,:) = mat6(2,:) + m43 * mat6(1,:) + m46 * mat6(6,:)
    endif
  endif

  orbit%vec(4) = orbit%vec(4) + k3l * rel_p * kz**2 * orbit%vec(3)**3 / 3
  if (field_ele%field_calc == helical_model$) then
    orbit%vec(2) = orbit%vec(2) + k3l * rel_p * kz**2 * orbit%vec(1)**3 / 3
  endif

  ! Quadrupole body

  if (logic_option(.false., make_matrix)) call mat_make_unit (kmat)

  if (field_ele%field_calc == helical_model$) then
    call quad_mat2_calc (k1yy,   step_len, rel_p, kmat(1:2,1:2), dz_x, ddz_x)
  else
    call quad_mat2_calc (k1x / rel_p**2, step_len, rel_p, kmat(1:2,1:2), dz_x, ddz_x)
  endif

  call quad_mat2_calc (k1yy, step_len, rel_p, kmat(3:4,3:4), dz_y, ddz_y)

  ! The mat6(i,6) terms are constructed so that mat6 is sympelctic

  if (logic_option(.false., make_matrix)) then
    if (any(orbit%vec(1:4) /= 0)) then
      kmat(5,1) = 2 * orbit%vec(1) * dz_x(1) +     orbit%vec(2) * dz_x(2)
      kmat(5,2) =     orbit%vec(1) * dz_x(2) + 2 * orbit%vec(2) * dz_x(3)
      kmat(5,3) = 2 * orbit%vec(3) * dz_y(1) +     orbit%vec(4) * dz_y(2)
      kmat(5,4) =     orbit%vec(3) * dz_y(2) + 2 * orbit%vec(4) * dz_y(3)
      kmat(5,6) = orbit%vec(1)**2 * ddz_x(1) + orbit%vec(1)*orbit%vec(2) * ddz_x(2) + orbit%vec(2)**2 * ddz_x(3) + &
                  orbit%vec(3)**2 * ddz_y(1) + orbit%vec(3)*orbit%vec(4) * ddz_y(2) + orbit%vec(4)**2 * ddz_y(3)  
    endif

    if (any(kmat(5,1:4) /= 0)) then
      kmat(1,6) = kmat(5,2) * kmat(1,1) - kmat(5,1) * kmat(1,2)
      kmat(2,6) = kmat(5,2) * kmat(2,1) - kmat(5,1) * kmat(2,2)
      kmat(3,6) = kmat(5,4) * kmat(3,3) - kmat(5,3) * kmat(3,4)

      kmat(4,6) = kmat(5,4) * kmat(4,3) - kmat(5,3) * kmat(4,4)
    endif

    mat6 = matmul(kmat, mat6)
  endif

  !

  orbit%vec(5) = orbit%vec(5) + &
                  dz_x(1) * orbit%vec(1)**2 + dz_x(2) * orbit%vec(1) * orbit%vec(2) + dz_x(3) * orbit%vec(2)**2 + &
                  dz_y(1) * orbit%vec(3)**2 + dz_y(2) * orbit%vec(3) * orbit%vec(4) + dz_y(3) * orbit%vec(4)**2 

  orbit%vec(1:2) = matmul(kmat(1:2,1:2), orbit%vec(1:2))
  orbit%vec(3:4) = matmul(kmat(3:4,3:4), orbit%vec(3:4))

  orbit%vec(5) = orbit%vec(5) + low_energy_z_correction (orbit, ele, step_len, mat6, make_matrix)

  ! 1/2 of the octupole kick.

  k3l = 2 * step_len * k1yy
  if (i == n_step) k3l = k3l / 2

  if (logic_option(.false., make_matrix)) then
    m43 = k3l * rel_p * kz**2 * orbit%vec(3)**2
    m46 = -k3l * kz**2 * orbit%vec(3)**3 / 3
    mat6(4,:) = mat6(4,:) + m43 * mat6(3,:) + m46 * mat6(6,:)
    if (field_ele%field_calc == helical_model$) then
      mat6(2,:) = mat6(2,:) + m43 * mat6(1,:) + m46 * mat6(6,:)
    endif
  endif

  orbit%vec(4) = orbit%vec(4) + k3l * rel_p * kz**2 * orbit%vec(3)**3 / 3
  if (field_ele%field_calc == helical_model$) then
    orbit%vec(2) = orbit%vec(2) + k3l * rel_p * kz**2 * orbit%vec(1)**3 / 3
  endif

  !

  if (i == n_step) then
    if (ix_mag_max > -1)  call ab_multipole_kicks (an,      bn,      ix_mag_max,  ele, orbit, magnetic$, r_step/2,   mat6, make_matrix)
    if (ix_elec_max > -1) call ab_multipole_kicks (an_elec, bn_elec, ix_elec_max, ele, orbit, electric$, step_len/2, mat6, make_matrix)
  else
    if (ix_mag_max > -1)  call ab_multipole_kicks (an,      bn,      ix_mag_max,  ele, orbit, magnetic$, r_step,   mat6, make_matrix)
    if (ix_elec_max > -1) call ab_multipole_kicks (an_elec, bn_elec, ix_elec_max, ele, orbit, electric$, step_len, mat6, make_matrix)
  endif
enddo

call offset_particle (ele, unset$, orbit, mat6 = mat6, make_matrix = make_matrix)

! Add in term to take care of the fact that the particle's motion undulates

orbit%t = t_start + orbit%direction*length / (c_light * beta_ref) + (z_start - orbit%vec(5)) / (c_light * orbit%beta)

if (field_ele%field_calc == helical_model$) then
  factor = length * (kz * ele%value(osc_amplitude$))**2 / 2 
else
  factor = length * (kz * ele%value(osc_amplitude$))**2 / 4
endif

orbit%t = orbit%t + factor / (c_light * orbit%beta * rel_p**2)
orbit%vec(5) = orbit%vec(5) + factor * (orbit%beta / beta_ref - 1 / rel_p**2)
if (logic_option(.false., make_matrix)) then
  gamma0 = ele%value(E_tot$) / mass_of(orbit%species)
  mat6(5,:) = mat6(5,:) + factor * ((orbit%beta**3/gamma0**2 + 2.0_rp) / rel_p**3) * mat6(6,:)
endif

end subroutine
