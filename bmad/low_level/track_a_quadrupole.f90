!+
! Subroutine track_a_quadrupole (orbit, ele, param, mat6, make_matrix)
!
! Bmad_standard tracking through a quadrupole element. 
!
! Input:
!   orbit       -- Coord_struct: Starting position.
!   ele         -- ele_struct: Quadrupole element.
!   param       -- lat_param_struct: Lattice parameters.
!   make_matrix -- logical, optional: Propagate the transfer matrix? Default is false.
!
! Output:
!   orbit      -- coord_struct: End position.
!   mat6(6,6)  -- real(rp), optional: Transfer matrix through the element.
!-

subroutine track_a_quadrupole (orbit, ele, param, mat6, make_matrix)

use bmad_interface, except_dummy => track_a_quadrupole

implicit none

type (coord_struct) :: orbit, start_orb
type (ele_struct), target :: ele
type (lat_param_struct) :: param
type (fringe_field_info_struct) fringe_info

real(rp), optional :: mat6(6,6)
real(rp) kmat(6,6), rel_p, dz_x(3), dz_y(3), ddz_x(3), ddz_y(3), length
real(rp) k1, b1, rel_tracking_charge, charge_dir, r_step, step_len, s_off, mass, e_tot
real(rp) an(0:n_pole_maxx), bn(0:n_pole_maxx), an_elec(0:n_pole_maxx), bn_elec(0:n_pole_maxx)

integer i, n_step, orientation, ix_mag_max, ix_elec_max

logical, optional :: make_matrix
logical drifting

!

start_orb = orbit
orientation = ele%orientation * start_orb%direction * start_orb%time_dir
rel_tracking_charge = rel_tracking_charge_to_mass(start_orb, param%particle)
charge_dir = rel_tracking_charge * orientation
length = orbit%time_dir * ele%value(l$)

! The b1 component is folded into the body focusing. With zero length there is no body so b1 is
! left in bn and applied with the other multipoles.
if (ele%value(l$) == 0) then
  call multipole_ele_to_ab (ele, .false., ix_mag_max,  an,      bn,      magnetic$, include_kicks$)
  b1 = 0
else
  call multipole_ele_to_ab (ele, .false., ix_mag_max,  an,      bn,      magnetic$, include_kicks$, b1)
endif
call multipole_ele_to_ab (ele, .false., ix_elec_max, an_elec, bn_elec, electric$)

n_step = 1
if (ix_mag_max > -1 .or. ix_elec_max > -1) n_step = max(nint(abs(length) / ele%value(ds_step$)), 1)

r_step = real(orbit%time_dir, rp) / n_step
step_len = ele%value(l$) * r_step
if (length == 0) n_step = 0

! Entrance edge

if (ele%bookkeeping_state%has_misalign) call offset_particle (ele, set$, &
                                    orbit, set_hvkicks = .false., mat6 = mat6, make_matrix = make_matrix)

call init_fringe_info (fringe_info, ele)
if (fringe_info%has_fringe) then
  fringe_info%particle_at = first_track_edge$
  call apply_element_edge_kick(orbit, fringe_info, ele, param, .false., mat6, make_matrix)
  if (orbit%state /= alive$) return
endif

! Multipole kicks. Notice that the magnetic multipoles have already been normalized by the length.

if (ix_mag_max > -1)  call ab_multipole_kicks (an,      bn,      ix_mag_max,  ele, orbit, magnetic$, r_step/2,   mat6, make_matrix)
if (ix_elec_max > -1) call ab_multipole_kicks (an_elec, bn_elec, ix_elec_max, ele, orbit, electric$, step_len/2, mat6, make_matrix)

! Body
! The exact Hamiltonian is split into H = H_par + H_kin. H_par is the paraxial quadrupole Hamiltonian
!   H_par = (px^2 + py^2) / (2 * (1 + pz)) + k1 * (x^2 - y^2) / 2
! which is solved by quad_mat2_calc. H_kin is the kinematic remainder
!   H_kin = -sqrt((1 + pz)^2 - px^2 - py^2) + (1 + pz) - (px^2 + py^2) / (2 * (1 + pz))
! which depends only upon the momenta so it is also solved exactly. H_kin is fourth order in px and py
! and gives the terms cubic in px and py. The two are combined with Simpson weights:
!   H_kin(L/6) H_par(L/2) H_kin(2L/3) H_par(L/2) H_kin(L/6)
! To first order in H_kin this integrates H_kin along the paraxial orbit with Simpson's rule so the
! error falls off as (k1 L^2)^2 instead of k1 L^2 for a trapezoid split. With k1 = 0, H_par and H_kin
! commute and the result is an exact drift.

do i = 1, n_step

  rel_p = 1 + orbit%vec(6)  ! Can change when there are electric fields
  k1 = charge_dir * b1 / (length * rel_p)

  call quad_mat2_calc (-k1, step_len/2, rel_p, kmat(1:2,1:2), dz_x, ddz_x)
  call quad_mat2_calc ( k1, step_len/2, rel_p, kmat(3:4,3:4), dz_y, ddz_y)

  call kinematic_correction (orbit, step_len/6)
  if (orbit%state /= alive$) return
  call paraxial_quad_body (orbit)
  call kinematic_correction (orbit, 2 * step_len / 3)
  if (orbit%state /= alive$) return
  call paraxial_quad_body (orbit)
  call kinematic_correction (orbit, step_len/6)
  if (orbit%state /= alive$) return

  orbit%vec(5) = orbit%vec(5) + low_energy_z_correction (orbit, ele, step_len, mat6, make_matrix)

  if (orbit%direction == -1) then
    orbit%vec(5) = orbit%vec(5) - 2.0_rp * c_light * orbit%beta * ele%value(delta_ref_time$)
  endif

  if (i == n_step) then
    if (ix_mag_max > -1)  call ab_multipole_kicks (an,      bn,      ix_mag_max,  ele, orbit, magnetic$, r_step/2,   mat6, make_matrix)
    if (ix_elec_max > -1) call ab_multipole_kicks (an_elec, bn_elec, ix_elec_max, ele, orbit, electric$, step_len/2, mat6, make_matrix)
  else
    if (ix_mag_max > -1)  call ab_multipole_kicks (an,      bn,      ix_mag_max,  ele, orbit, magnetic$, r_step,   mat6, make_matrix)
    if (ix_elec_max > -1) call ab_multipole_kicks (an_elec, bn_elec, ix_elec_max, ele, orbit, electric$, step_len, mat6, make_matrix)
  endif

enddo

! Zero length: the body loop is skipped so the second half of the multipole kick is applied here.

if (n_step == 0) then
  if (ix_mag_max > -1)  call ab_multipole_kicks (an,      bn,      ix_mag_max,  ele, orbit, magnetic$, r_step/2,   mat6, make_matrix)
endif

! Exit edge

if (fringe_info%has_fringe) then
  fringe_info%particle_at = second_track_edge$
  call apply_element_edge_kick(orbit, fringe_info, ele, param, .false., mat6, make_matrix)
  if (orbit%state /= alive$) return
endif

if (ele%bookkeeping_state%has_misalign) call offset_particle (ele, unset$, &
                                  orbit, set_hvkicks = .false., mat6 = mat6, make_matrix = make_matrix)

orbit%t = start_orb%t + orbit%direction*orbit%time_dir*ele%value(delta_ref_time$) + (start_orb%vec(5) - orbit%vec(5)) / (orbit%beta * c_light)

!---------------------------------------------------------------------------------------------
contains

! Paraxial quadrupole map through the length used to make kmat(1:4,1:4), dz_x, etc.

subroutine paraxial_quad_body (orbit)

type (coord_struct) orbit
real(rp) r1(6), r3(6), r5(6)

! The mat6(i,6) terms are constructed so that mat6 is symplectic

if (logic_option(.false., make_matrix)) then
  kmat(5,1:4) = 0
  kmat(5,6) = 0
  kmat(1:4,6) = 0

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

  ! kmat is block diagonal in (1:2) and (3:4) plus row 5 and column 6. Row 6 is unity.
  ! For speed, apply directly to the rows of mat6 instead of using a full 6x6 matmul.
  r5 = mat6(5,:) + kmat(5,1) * mat6(1,:) + kmat(5,2) * mat6(2,:) + kmat(5,3) * mat6(3,:) + &
                   kmat(5,4) * mat6(4,:) + kmat(5,6) * mat6(6,:)   ! kmat(5,5) = 1
  r1 = mat6(1,:);  r3 = mat6(3,:)
  mat6(1,:) = kmat(1,1) * r1 + kmat(1,2) * mat6(2,:) + kmat(1,6) * mat6(6,:)
  mat6(2,:) = kmat(2,1) * r1 + kmat(2,2) * mat6(2,:) + kmat(2,6) * mat6(6,:)
  mat6(3,:) = kmat(3,3) * r3 + kmat(3,4) * mat6(4,:) + kmat(3,6) * mat6(6,:)
  mat6(4,:) = kmat(4,3) * r3 + kmat(4,4) * mat6(4,:) + kmat(4,6) * mat6(6,:)
  mat6(5,:) = r5
endif

!

orbit%vec(5) = orbit%vec(5) + &
               (dz_x(1) * orbit%vec(1)**2 + dz_x(2) * orbit%vec(1) * orbit%vec(2) + dz_x(3) * orbit%vec(2)**2 + &
                dz_y(1) * orbit%vec(3)**2 + dz_y(2) * orbit%vec(3) * orbit%vec(4) + dz_y(3) * orbit%vec(4)**2)

orbit%vec(1:2) = matmul(kmat(1:2,1:2), orbit%vec(1:2))
orbit%vec(3:4) = matmul(kmat(3:4,3:4), orbit%vec(3:4))

end subroutine paraxial_quad_body

!---------------------------------------------------------------------------------------------
! Exact map for H_kin (see above) through a length ds.
! H_kin depends only upon px, py, and pz so only x, y, and z change.

subroutine kinematic_correction (orbit, ds)

type (coord_struct) orbit
real(rp) ds, rel_p, px_rel, py_rel, pxy2, ps, g, f, h
real(rp) m12, m14, m16, m32, m34, m36, m52, m54, m56

!

rel_p = 1 + orbit%vec(6)
px_rel = orbit%vec(2) / rel_p
py_rel = orbit%vec(4) / rel_p
pxy2 = px_rel**2 + py_rel**2

if (pxy2 >= 1) then
  orbit%state = lost_pz$
  return
endif

ps = sqrt(1 - pxy2)
g = pxy2 / (ps * (1 + ps))     ! = 1/ps - 1 without round-off.

orbit%vec(1) = orbit%vec(1) + ds * px_rel * g
orbit%vec(3) = orbit%vec(3) + ds * py_rel * g
orbit%vec(5) = orbit%vec(5) - ds * pxy2**2 * (2 + ps) / (2 * ps * (1 + ps)**2)  ! = ds * (1 - 1/ps + pxy2/2)

if (logic_option(.false., make_matrix)) then
  f = 1 / ps**3
  h = -pxy2 * (1 + ps + ps**2) / ((1 + ps) * ps**3)   ! = 1 - 1/ps^3 without round-off.
  m12 =  ds * (g + px_rel**2 * f) / rel_p
  m14 =  ds * px_rel * py_rel * f / rel_p
  m16 = -ds * px_rel * (g + pxy2 * f) / rel_p
  m32 =  m14
  m34 =  ds * (g + py_rel**2 * f) / rel_p
  m36 = -ds * py_rel * (g + pxy2 * f) / rel_p
  m52 =  ds * px_rel * h / rel_p
  m54 =  ds * py_rel * h / rel_p
  m56 = -ds * pxy2 * h / rel_p

  ! As with a drift, only rows 1, 3, and 5 change and they depend only upon rows 2, 4, and 6.
  mat6(1,:) = mat6(1,:) + m12 * mat6(2,:) + m14 * mat6(4,:) + m16 * mat6(6,:)
  mat6(3,:) = mat6(3,:) + m32 * mat6(2,:) + m34 * mat6(4,:) + m36 * mat6(6,:)
  mat6(5,:) = mat6(5,:) + m52 * mat6(2,:) + m54 * mat6(4,:) + m56 * mat6(6,:)
endif

end subroutine kinematic_correction

end subroutine
