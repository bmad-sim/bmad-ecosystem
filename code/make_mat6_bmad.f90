!+
! Subroutine make_mat6_bmad (ele, param, orb_in, orb_out, end_in, err)
!
! Subroutine to make the 6x6 transfer matrix for an element. 
!
! Modules needed:
!   use bmad
!
! Input:
!   ele    -- Ele_struct: Element with transfer matrix
!   param  -- lat_param_struct: Parameters are needed for some elements.
!   orb_in -- Coord_struct: Coordinates at the beginning of element. 
!   end_in -- Logical, optional: If present and True then the end coords orb_out
!               will be taken as input. Not output as normal.
!
! Output:
!   ele       -- Ele_struct: Element with transfer matrix.
!     %vec0     -- 0th order map component
!     %mat6     -- 6x6 transfer matrix.
!   orb_out   -- Coord_struct: Coordinates at the end of element.
!   err       -- Logical, optional: Set True if there is an error. False otherwise.
!-

subroutine make_mat6_bmad (ele, param, orb_in, orb_out, end_in, err)

use track1_mod, dummy1 => make_mat6_bmad

implicit none

type (ele_struct), target :: ele
type (ele_struct) :: temp_ele1, temp_ele2
type (coord_struct) :: orb_in, orb_out
type (coord_struct) :: c00, orb_out1, c_int
type (coord_struct) orb, c0_off, orb_out_off
type (lat_param_struct)  param

real(rp), pointer :: mat6(:,:), v(:)

real(rp) mat6_pre(6,6), mat6_post(6,6), mat6_i(6,6)
real(rp) mat4(4,4), m2(2,2), kmat4(4,4), om_g, om, om_g2
real(rp) angle, k1, ks, length, e2, g, g_err, coef
real(rp) k2l, k3l, beta_ref, c_min, c_plu, dc_min, dc_plu
real(rp) t5_22, t5_33, t5_34, t5_44
real(rp) factor, kmat6(6,6), drift(6,6), ww(3,3)
real(rp) s_pos, s_pos_old, dr(3), axis(3), w_mat(3,3)
real(rp) knl(0:n_pole_maxx), tilt(0:n_pole_maxx)
real(rp) an_elec(0:n_pole_maxx), bn_elec(0:n_pole_maxx)
real(rp) c_e, c_m, gamma_old, gamma_new, voltage, sqrt_8
real(rp) arg, rel_p, rel_p2, dp_dg, dp_dg_dz1, dp_dg_dpz1
real(rp) cy, sy, k2, s_off, x_pitch, y_pitch, y_ave, k_z, stg, one_ctg
real(rp) dz_x(3), dz_y(3), xp_start, yp_start
real(rp) k, L, m55, m65, m66, new_pc, new_beta, dbeta_dpz
real(rp) cos_phi, sin_phi, cos_term, dcos_phi, gradient_net, e_start, e_end, e_ratio, pc, p0c
real(rp) alpha, sin_a, cos_a, fg, phase0, phase, t0, dt_ref, E, pxy2, dE
real(rp) g_tot, ct, st, x, px, y, py, z, pz, p_s, Dxy, Dy, px_t
real(rp) Dxy_t, dpx_t, df_dp, kx_1, ky_1, kx_2, ky_2
real(rp) mc2, pc_start, pc_end, pc_start_ref, pc_end_ref, gradient_max, voltage_max
real(rp) beta_start, beta_end, p_rel, beta_rel, xp, yp, s_ent, ds_ref
real(rp) drp1_dr0, drp1_drp0, drp2_dr0, drp2_drp0, xp1, xp2, yp1, yp2
real(rp) dp_long_dpx, dp_long_dpy, dp_long_dpz, dalpha_dpx, dalpha_dpy, dalpha_dpz
real(rp) Dy_dpy, Dy_dpz, dpx_t_dx, dpx_t_dpx, dpx_t_dpy, dpx_t_dpz, dp_ratio
real(rp) df_dx, df_dpx, df_dpz, deps_dx, deps_dpx, deps_dpy, deps_dpz
real(rp) ps, dps_dpx, dps_dpy, dps_dpz, dE_rel_dpz, dps_dx, sinh_c, cosh_c, ff 
real(rp) hk, vk, k_E, E_tot, E_rel, p_factor, sinh_k, cosh1_k, rel_tracking_charge, rel_charge

integer i, n_slice, key, dir, ix_pole_max, tm

real(rp) charge_dir, hkick, vkick, kick

logical, optional :: end_in, err
logical err_flag, fringe_here, drifting, do_track, set_tilt
character(16), parameter :: r_name = 'make_mat6_bmad'

!--------------------------------------------------------
! init

if (present(err)) err = .false.

mat6 => ele%mat6
v => ele%value

call mat_make_unit (mat6)

! If element is off.

key = ele%key

if (.not. ele%is_on) then
  select case (key)
  case (taylor$, match$, fiducial$, floor_shift$)
    if (.not. logic_option (.false., end_in)) call set_orb_out (orb_out, orb_in)
    return
  case (ab_multipole$, multipole$, lcavity$, sbend$, patch$)
    ! Nothing to do here
  case default
    key = drift$  
  end select
endif

if (key == sol_quad$ .and. v(k1$) == 0) key = solenoid$

!

select case (key)
case (sad_mult$, match$, beambeam$, sbend$, patch$, quadrupole$, drift$, &
      rcollimator$, ecollimator$, monitor$, instrument$, pipe$, kicker$, hkicker$, vkicker$, &
      elseparator$, rfcavity$, lcavity$)
  tm = ele%tracking_method
  if (key /= wiggler$ .or. ele%sub_key /= map_type$)   ele%tracking_method = bmad_standard$
  call track1 (orb_in, ele, param, c00, mat6 = mat6, make_matrix = .true.)
  ele%tracking_method = tm
  ele%vec0 = c00%vec - matmul(mat6, orb_in%vec)
  if (.not. logic_option (.false., end_in)) call set_orb_out (orb_out, c00)
  return
end select

!

ele%vec0 = 0

length = v(l$)
rel_p = 1 + orb_in%vec(6) 
rel_tracking_charge = rel_tracking_charge_to_mass(orb_in, param)
charge_dir = rel_tracking_charge * ele%orientation
c00 = orb_in
!! c00%direction = +1  ! Quad calc, for example, will not be correct if this is set.

! Note: sad_mult, match, etc. will handle the calc of orb_out if needed.

do_track = (.not. logic_option (.false., end_in))

if (do_track) then
  tm = ele%tracking_method
  if (key /= wiggler$ .or. ele%sub_key /= map_type$)   ele%tracking_method = bmad_standard$
  if (ele%tracking_method == linear$) then
    c00%state = alive$
    call track1_bmad (c00, ele, param, orb_out)
  else
    call track1 (c00, ele, param, orb_out)
  endif
  ele%tracking_method = tm

  if (orb_out%state /= alive$) then
    mat6 = 0
    if (present(err)) err = .true.
    call out_io (s_error$, r_name, 'PARTICLE LOST IN TRACKING AT: ' // trim(ele%name) // '  (\i0\) ', &
                 i_array = [ele%ix_ele] )
    return
  endif
endif

orb_out1 = orb_out

!--------------------------------------------------------
! Selection...

select case (key)

!--------------------------------------------------------
! Custom

case (custom$)

  if (present(err)) err = .true.
  call out_io (s_fatal$, r_name,  'MAT6_CALC_METHOD = BMAD_STANDARD IS NOT ALLOWED FOR A CUSTOM ELEMENT: ' // ele%name)
  if (global_com%exit_on_error) call err_exit
  return

!--------------------------------------------------------
! Marker, branch, photon_branch, etc.

case (marker$, detector$, fork$, photon_fork$, floor_shift$, fiducial$, mask$) 
  return

!--------------------------------------------------------
! Match

case (match$)
  call match_ele_to_mat6 (ele, err_flag)
  if (present(err)) err = err_flag
  if (err_flag) return
  if (.not. logic_option (.false., end_in)) then
    call track1_bmad (c00, ele, param, orb_out)

    ! If the particle is lost with a match element with match_end set to True,
    ! the problem is most likely that twiss_propagate_all has not yet
    ! been called (so the previous element's Twiss parameters are not yet set). 
    ! In this case, ignore a lost particle.

    if (orb_out%state /= alive$ .and. is_false(v(match_end$))) then
      call out_io (s_error$, r_name, 'PARTICLE LOST IN TRACKING AT: ' // trim(ele%name) // '  (\i0\) ', &
                   i_array = [ele%ix_ele] )
    endif
  endif
  return

!--------------------------------------------------------
! Multipole, AB_Multipole

case (multipole$, ab_multipole$)

  if (.not. ele%multipoles_on) return

  call offset_particle (ele, param, set$, c00, set_tilt = .false.)

  call multipole_ele_to_kt (ele, .true., ix_pole_max, knl, tilt)
  if (ix_pole_max > -1) then
    call multipole_kick_mat (knl, tilt, param%particle, ele, c00, 1.0_rp, ele%mat6)

    ! if knl(0) is non-zero then the reference orbit itself is bent
    ! and we need to account for this.

    if (knl(0) /= 0 .and. ele%key == multipole$) then
      ele%mat6(2,6) = knl(0) * cos(tilt(0))
      ele%mat6(4,6) = knl(0) * sin(tilt(0))
      ele%mat6(5,1) = -ele%mat6(2,6)
      ele%mat6(5,3) = -ele%mat6(4,6)
    endif
  endif

  ele%vec0 = orb_out%vec - matmul(mat6, orb_in%vec)

!--------------------------------------------------------
! Octupole
! the octupole is modeled as kick-drift-kick

case (octupole$)

  call offset_particle (ele, param, set$, c00) 
  n_slice = max(1, nint(length / v(ds_step$)))

  do i = 0, n_slice
    k3l = charge_dir * v(k3$) * length / n_slice
    if (i == 0 .or. i == n_slice) k3l = k3l / 2
    call mat4_multipole (k3l, 0.0_rp, 3, c00, kmat4)
    c00%vec(2) = c00%vec(2) + k3l * (3*c00%vec(1)*c00%vec(3)**2 - c00%vec(1)**3) / 6
    c00%vec(4) = c00%vec(4) + k3l * (3*c00%vec(3)*c00%vec(1)**2 - c00%vec(3)**3) / 6
    mat6(1:4,1:6) = matmul(kmat4, mat6(1:4,1:6))
    if (i /= n_slice) then
      call drift_mat6_calc (drift, length/n_slice, ele, param, c00)
      call track_a_drift (c00, length/n_slice)
      mat6 = matmul(drift,mat6)
    end if
  end do

  if (v(tilt_tot$) /= 0) then
    call tilt_mat6 (mat6, v(tilt_tot$))
  endif

  call add_multipoles_and_z_offset (.true.)
  ele%vec0 = orb_out%vec - matmul(mat6, orb_in%vec)

!--------------------------------------------------------
! rbends are not allowed internally

case (rbend$)

  if (present(err)) err = .true.
  call out_io (s_fatal$, r_name,  'RBEND ELEMENTS NOT ALLOWED INTERNALLY!')
  if (global_com%exit_on_error) call err_exit
  return

!--------------------------------------------------------
! Sextupole.
! the sextupole is modeled as kick-drift-kick

case (sextupole$)

  call offset_particle (ele, param, set$, c00)
  call hard_multipole_edge_kick (ele, param, first_track_edge$, c00, mat6, .true.)

  n_slice = max(1, nint(length / v(ds_step$)))
  
  do i = 0, n_slice
    k2l = charge_dir * v(k2$) * length / n_slice
    if (i == 0 .or. i == n_slice) k2l = k2l / 2
    call mat4_multipole (k2l, 0.0_rp, 2, c00, kmat4)
    c00%vec(2) = c00%vec(2) + k2l * (c00%vec(3)**2 - c00%vec(1)**2)/2
    c00%vec(4) = c00%vec(4) + k2l * c00%vec(1) * c00%vec(3)
    mat6(1:4,1:6) = matmul(kmat4,mat6(1:4,1:6))
    if (i /= n_slice) then
      call drift_mat6_calc (drift, length/n_slice, ele, param, c00)
      call track_a_drift (c00, length/n_slice)
      mat6 = matmul(drift,mat6)
    end if
  end do

  call hard_multipole_edge_kick (ele, param, second_track_edge$, c00, mat6, .true.)

  if (v(tilt_tot$) /= 0) then
    call tilt_mat6 (mat6, v(tilt_tot$))
  endif

  call add_multipoles_and_z_offset (.true.)
  ele%vec0 = orb_out%vec - matmul(mat6, orb_in%vec)

!--------------------------------------------------------
! solenoid

case (solenoid$)

  call offset_particle (ele, param, set$, c00)
  call solenoid_track_and_mat (ele, length, param, c00, c00, mat6)
  call offset_particle (ele, param, unset$, c00)

  call add_multipoles_and_z_offset (.true.)
  ele%vec0 = c00%vec - matmul(mat6, orb_in%vec)

!--------------------------------------------------------
! solenoid/quad

case (sol_quad$)

  call offset_particle (ele, param, set$, c00)

  call sol_quad_mat6_calc (v(ks$) * rel_tracking_charge, v(k1$) * charge_dir, length, c00%vec, mat6)

  if (v(tilt_tot$) /= 0) then
    call tilt_mat6 (mat6, v(tilt_tot$))
  endif

  call add_multipoles_and_z_offset (.true.)
  call add_M56_low_E_correction()
  ele%vec0 = orb_out%vec - matmul(mat6, orb_in%vec)

!--------------------------------------------------------
! taylor

case (taylor$)

  call make_mat6_taylor (ele, param, orb_in)

!--------------------------------------------------------
! wiggler

case (wiggler$, undulator$)

  call offset_particle (ele, param, set$, c00)
  call offset_particle (ele, param, set$, orb_out1, ds_pos = length)

  call mat_make_unit (mat6)     ! make a unit matrix

  if (length == 0) then
    call add_multipoles_and_z_offset (.true.)
  call add_M56_low_E_correction()
    return
  endif

  k1 = -0.5 * charge_dir * (c_light * v(b_max$) / (v(p0c$) * rel_p))**2

  ! octuple correction to k1

  y_ave = (c00%vec(3) + orb_out1%vec(3)) / 2
  if (v(l_pole$) == 0) then
    k_z = 0
  else
    k_z = pi / v(l_pole$)
  endif
  k1 = k1 * (1 + 2 * (k_z * y_ave)**2)

  !

  mat6(1, 1) = 1
  mat6(1, 2) = length / rel_p
  mat6(2, 1) = 0
  mat6(2, 2) = 1

  call quad_mat2_calc (k1, length, rel_p, mat6(3:4,3:4))

  cy = mat6(3, 3)
  sy = mat6(3, 4)

  t5_22 = -length / 2
  t5_33 =  k1 * (length - sy*cy) / 4
  t5_34 = -k1 * sy**2 / 2
  t5_44 = -(length + sy*cy) / 4

  ! the mat6(i,6) terms are constructed so that mat6 is sympelctic

  mat6(5,2) = 2 * c00%vec(2) * t5_22 / rel_p
  mat6(5,3) = 2 * c00%vec(3) * t5_33 +     c00%vec(4) * t5_34 / rel_p
  mat6(5,4) =     c00%vec(3) * t5_34 + 2 * c00%vec(4) * t5_44 / rel_p

  mat6(1,6) = mat6(5,2) * mat6(1,1)
  mat6(2,6) = mat6(5,2) * mat6(2,1)
  mat6(3,6) = mat6(5,4) * mat6(3,3) - mat6(5,3) * mat6(3,4)
  mat6(4,6) = mat6(5,4) * mat6(4,3) - mat6(5,3) * mat6(4,4)

  if (v(tilt_tot$) /= 0) then
    call tilt_mat6 (mat6, v(tilt_tot$))
  endif

  call add_multipoles_and_z_offset (.true.)
  call add_M56_low_E_correction()
  ele%vec0 = orb_out%vec - matmul(mat6, orb_in%vec)

!--------------------------------------------------------
! unrecognized element

case default

  if (present(err)) err = .true.
  call out_io (s_fatal$, r_name,  'UNKNOWN ELEMENT KEY: \i0\ ', &
                                  'FOR ELEMENT: ' // ele%name, i_array = [ele%key])
  if (global_com%exit_on_error) call err_exit
  return

end select

!--------------------------------------------------------
contains

subroutine set_orb_out (orb_out, c00)

type (coord_struct) orb_out, c00

orb_out = c00
if (orb_out%direction == 1) then
  orb_out%s = ele%s
else
  orb_out%s = ele%s_start
endif

end subroutine set_orb_out

!--------------------------------------------------------
! contains

! put in multipole components

subroutine add_multipoles_and_z_offset (add_pole)

real(rp) mat6_m(6,6)
integer ix_pole_max
logical add_pole

!

if (add_pole) then
  call multipole_ele_to_kt (ele, .true., ix_pole_max, knl, tilt)
  if (ix_pole_max > -1) then
    knl = knl * ele%orientation
    call multipole_kick_mat (knl, tilt, param%particle, ele, orb_in, 0.5_rp, mat6_m)
    mat6(:,1) = mat6(:,1) + mat6(:,2) * mat6_m(2,1) + mat6(:,4) * mat6_m(4,1)
    mat6(:,3) = mat6(:,3) + mat6(:,2) * mat6_m(2,3) + mat6(:,4) * mat6_m(4,3)
    call multipole_kick_mat (knl, tilt, param%particle, ele, orb_out, 0.5_rp, mat6_m)
    mat6(2,:) = mat6(2,:) + mat6_m(2,1) * mat6(1,:) + mat6_m(2,3) * mat6(3,:)
    mat6(4,:) = mat6(4,:) + mat6_m(4,1) * mat6(1,:) + mat6_m(4,3) * mat6(3,:)
  endif
endif

if (v(z_offset_tot$) /= 0) then
  s_off = v(z_offset_tot$) * ele%orientation
  mat6(1,:) = mat6(1,:) - s_off * mat6(2,:)
  mat6(3,:) = mat6(3,:) - s_off * mat6(4,:)
  mat6(:,2) = mat6(:,2) + mat6(:,1) * s_off
  mat6(:,4) = mat6(:,4) + mat6(:,3) * s_off
endif

! pitch corrections

call mat6_add_pitch (v(x_pitch_tot$), v(y_pitch_tot$), ele%orientation, ele%mat6)

end subroutine add_multipoles_and_z_offset

!----------------------------------------------------------------
! contains

subroutine add_M56_low_E_correction()

real(rp) mass, e_tot

! 1/gamma^2 m56 correction

mass = mass_of(orb_in%species)
e_tot = v(p0c$) * (1 + orb_in%vec(6)) / orb_in%beta
mat6(5,6) = mat6(5,6) + length * mass**2 * v(e_tot$) / e_tot**3

end subroutine add_M56_low_E_correction

end subroutine make_mat6_bmad
