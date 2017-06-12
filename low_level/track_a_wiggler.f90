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

use track1_mod, except_dummy => track_a_wiggler

implicit none

type (coord_struct) :: orbit
type (ele_struct), target :: ele
type (lat_param_struct) :: param

real(rp), optional :: mat6(6,6)
real(rp) mat2(2,2), z_start, beta_ref, p_factor, k1_factor, k1, k1l, length, k_z, rel_p, t_start
real(rp) an(0:n_pole_maxx), bn(0:n_pole_maxx), an_elec(0:n_pole_maxx), bn_elec(0:n_pole_maxx)
real(rp) m43, m46, m52, m54, m56, mc2_rel, kmat(6,6)
real(rp) dz_x(3), dz_y(3), ddz_x(3), ddz_y(3)

integer ix_pole_max, ix_elec_max

logical, optional :: make_matrix

! Only periodic type wigglers are handled here.
! In the horizontal plane the tracking looks like a drift.
! The tracking in the vertical plane is:
!   1) 1/2 the octupole kick at the entrance face.
!   2) Track as a quadrupole through the body
!   3) 1/2 the octupole kick at the exit face.

call multipole_ele_to_ab (ele, .false., ix_pole_max, an,      bn,      magnetic$, include_kicks = .true.)
call multipole_ele_to_ab (ele, .false., ix_elec_max, an_elec, bn_elec, electric$)

!

call offset_particle (ele, param, set$, orbit, mat6 = mat6, make_matrix = make_matrix)

if (ix_pole_max > -1) call ab_multipole_kicks (an,      bn,      param%particle, ele, orbit, magnetic$, 1.0_rp/2,   mat6, make_matrix)
if (ix_elec_max > -1) call ab_multipole_kicks (an_elec, bn_elec, param%particle, ele, orbit, electric$, length/2, mat6, make_matrix)

! Notice that the averageed quadrupole and octupole kicks are independent of the sign of the particle charge and independent 
! of the longitudinal direction of travel of the particle!

z_start = orbit%vec(5)
t_start = orbit%t
beta_ref = ele%value(p0c$) / ele%value(e_tot$)
length = ele%value(l$)
rel_p = 1 + orbit%vec(6)
mc2_rel = mass_of(orbit%species) / orbit%p0c

!

if (ele%value(l_pole$) == 0) then
  k_z = 1d100    ! Something large
else
  k_z = pi / ele%value(l_pole$)
endif

k1_factor = -abs(rel_tracking_charge_to_mass(orbit, param)) * 0.5 * (c_light * ele%value(b_max$) / ele%value(p0c$))**2
k1 = k1_factor / rel_p**2
k1l = length * k1

! 1/2 of the octupole kick at the entrance face.

if (logic_option(.false., make_matrix)) then
  m43 = k1l * rel_p * k_z**2 * orbit%vec(3)**2
  m46 = -k1l * k_z**2 * orbit%vec(3)**3 / 3
  mat6(4,:) = mat6(4,:) + m43 * mat6(3,:) + m46 * mat6(6,:)
endif

orbit%vec(4) = orbit%vec(4) + k1l * rel_p * k_z**2 * orbit%vec(3)**3 / 3

! Quadrupole body

if (logic_option(.false., make_matrix)) call mat_make_unit (kmat)

call quad_mat2_calc (0.0_rp, length, rel_p, kmat(1:2,1:2), dz_x, ddz_x)
call quad_mat2_calc (k1,     length, rel_p, kmat(3:4,3:4), dz_y, ddz_y)

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

orbit%vec(5) = orbit%vec(5) + low_energy_z_correction (orbit, ele, length, mat6, make_matrix)

! 1/2 of the octupole kick at the exit face.

if (logic_option(.false., make_matrix)) then
  m43 = k1l * rel_p * k_z**2 * orbit%vec(3)**2
  m46 = -k1l * k_z**2 * orbit%vec(3)**3 / 3
  mat6(4,:) = mat6(4,:) + m43 * mat6(3,:) + m46 * mat6(6,:)
endif

orbit%vec(4) = orbit%vec(4) + k1l * rel_p * k_z**2 * orbit%vec(3)**3 / 3

!

if (ix_pole_max > -1) call ab_multipole_kicks (an,      bn,      param%particle, ele, orbit, magnetic$, 1.0_rp/2,   mat6, make_matrix)
if (ix_elec_max > -1) call ab_multipole_kicks (an_elec, bn_elec, param%particle, ele, orbit, electric$, length/2, mat6, make_matrix)

call offset_particle (ele, param, unset$, orbit, mat6 = mat6, make_matrix = make_matrix)

orbit%t = t_start + length / (c_light * beta_ref) + (z_start - orbit%vec(5)) / (c_light * orbit%beta)

end subroutine
