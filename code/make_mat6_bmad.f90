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

use sad_mod, dummy => make_mat6_bmad
use mad_mod, dummy1 => make_mat6_bmad

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
real(rp) k2l, k3l, del_l, beta_ref, c_min, c_plu, dc_min, dc_plu
real(rp) t5_22, t5_33, t5_34, t5_44
real(rp) factor, kmat6(6,6), drift(6,6), ww(3,3)
real(rp) s_pos, s_pos_old, z_slice(100), dr(3), axis(3), w_mat(3,3)
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
real(rp) dbeta1_dE1, dbeta2_dE2, dalpha_dt1, dalpha_dE1, dcoef_dt1, dcoef_dE1, z21, z22
real(rp) drp1_dr0, drp1_drp0, drp2_dr0, drp2_drp0, xp1, xp2, yp1, yp2
real(rp) dp_long_dpx, dp_long_dpy, dp_long_dpz, dalpha_dpx, dalpha_dpy, dalpha_dpz
real(rp) Dy_dpy, Dy_dpz, dpx_t_dx, dpx_t_dpx, dpx_t_dpy, dpx_t_dpz, dp_ratio
real(rp) df_dx, df_dpx, df_dpz, deps_dx, deps_dpx, deps_dpy, deps_dpz
real(rp) ps, dps_dpx, dps_dpy, dps_dpz, dE_rel_dpz, dps_dx, sinh_c, cosh_c, ff 
real(rp) hk, vk, k_E, E_tot, E_rel, p_factor, sinh_k, cosh1_k, rel_tracking_charge, rel_charge

integer i, n_slice, key, dir

real(rp) charge_dir, hkick, vkick, kick

logical, optional :: end_in, err
logical err_flag, has_nonzero_pole, has_nonzero_elec, fringe_here, drifting, do_track
character(16), parameter :: r_name = 'make_mat6_bmad'

!--------------------------------------------------------
! init

if (present(err)) err = .false.

mat6 => ele%mat6
v => ele%value

call mat_make_unit (mat6)
ele%vec0 = 0

length = v(l$)
rel_p = 1 + orb_in%vec(6) 
key = ele%key
rel_tracking_charge = relative_tracking_charge(orb_in, param)
charge_dir = rel_tracking_charge * ele%orientation
c00 = orb_in
c00%direction = +1

! Note: sad_mult, match, etc. will handle the calc of orb_out if needed.

do_track = (.not. logic_option (.false., end_in))
if (ele%key == sad_mult$ .or. ele%key == match$) do_track = .false.

if (do_track) then
  if (ele%tracking_method == linear$) then
    c00%state = alive$
    call track1_bmad (c00, ele, param, orb_out)
  else
    call track1 (c00, ele, param, orb_out)
  endif

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
! If element is off.

if (.not. ele%is_on) then
  select case (key)
  case (taylor$, match$, fiducial$, floor_shift$)
    return
  case (ab_multipole$, multipole$, lcavity$, sbend$, patch$)
    ! Nothing to do here
  case default
    key = drift$  
  end select
endif

!--------------------------------------------------------
! Selection...

if (key == sol_quad$ .and. v(k1$) == 0) key = solenoid$

select case (key)

!--------------------------------------------------------
! beam-beam interaction

case (beambeam$)

  call offset_particle (ele, param, set$, c00)
  call canonical_to_angle_coords (c00)

  call offset_particle (ele, param, set$, orb_out1, ds_pos = length)
  call canonical_to_angle_coords (orb_out1)

  n_slice = nint(v(n_slice$))
  if (n_slice < 1) then
    if (present(err)) err = .true.
    call out_io (s_fatal$, r_name,  'N_SLICE FOR BEAMBEAM ELEMENT IS NEGATIVE')
    call type_ele (ele, .true., 0, .false., 0, .false.)
    return
  endif

  if (v(charge$) == 0 .or. param%n_part == 0) return

  ! factor of 2 in orb%vec(5) since relative motion of the two beams is 2*c_light

  if (n_slice == 1) then
    call bbi_kick_matrix (ele, param, c00, 0.0_rp, mat6)
  else
    call bbi_slice_calc (ele, n_slice, z_slice)

    s_pos = 0          ! start at IP
    orb = c00
    orb%vec(2) = c00%vec(2) - v(x_pitch_tot$)
    orb%vec(4) = c00%vec(4) - v(y_pitch_tot$)
    call mat_make_unit (mat4)

    do i = 1, n_slice + 1
      s_pos_old = s_pos  ! current position
      s_pos = (z_slice(i) + c00%vec(5)) / 2 ! position of slice relative to IP
      del_l = s_pos - s_pos_old
      mat4(1,1:4) = mat4(1,1:4) + del_l * mat4(2,1:4)
      mat4(3,1:4) = mat4(3,1:4) + del_l * mat4(4,1:4)
      if (i == n_slice + 1) exit
      orb%vec(1) = c00%vec(1) + s_pos * orb%vec(2)
      orb%vec(3) = c00%vec(3) + s_pos * orb%vec(4)
      call bbi_kick_matrix (ele, param, orb, s_pos, kmat6)
      mat4(2,1:4) = mat4(2,1:4) + kmat6(2,1) * mat4(1,1:4) + &
                                  kmat6(2,3) * mat4(3,1:4)
      mat4(4,1:4) = mat4(4,1:4) + kmat6(4,1) * mat4(1,1:4) + &
                                  kmat6(4,3) * mat4(3,1:4)
    enddo

    mat6(1:4,1:4) = mat4

  endif

  call add_multipoles_and_z_offset (.true.)
  ele%vec0 = orb_out%vec - matmul(mat6, orb_in%vec)

!--------------------------------------------------------
! Custom

case (custom$)

  if (present(err)) err = .true.
  call out_io (s_fatal$, r_name,  'MAT6_CALC_METHOD = BMAD_STANDARD IS NOT ALLOWED FOR A CUSTOM ELEMENT: ' // ele%name)
  if (global_com%exit_on_error) call err_exit
  return

!-----------------------------------------------
! Drift

case (drift$) 
  call drift_mat6_calc (mat6, length, ele, param, c00)
  ele%vec0 = orb_out%vec - matmul(mat6, orb_in%vec)
  return

!-----------------------------------------------
! Elseparator

case (elseparator$)

  call offset_particle (ele, param, set$, c00, set_hvkicks = .false.)

  ! Compute kick

  rel_charge = abs(rel_tracking_charge) * sign(1, charge_of(orb_in%species))
  hk = ele%value(hkick$) * rel_charge
  vk = ele%value(vkick$) * rel_charge

  if (length == 0) length = 1d-50  ! To avoid divide by zero
  k_E = sqrt(hk**2 + vk**2) / length

  ! Rotate (x, y) so that kick is in +x direction.

  angle = atan2(vk, hk)
  call tilt_coords (angle, c00%vec)

  !

  E_tot = ele%value(p0c$) * (1 + c00%vec(6)) / c00%beta 
  E_rel = E_tot / ele%value(p0c$)
  mc2 = mass_of(orb_in%species)

  x = c00%vec(1)
  px = c00%vec(2)
  p_factor = (mc2 / ele%value(p0c$))**2 + c00%vec(2)**2 + c00%vec(4)**2

  ps = sqrt((E_rel + k_E * x)**2 - p_factor)
  alpha = length / ps
  coef = k_E * length / ps

  dE_rel_dpz = c00%beta
  dps_dx  = k_E * (E_rel + k_E * x) / ps
  dps_dpx = -c00%vec(2) / ps
  dps_dpy = -c00%vec(4) / ps
  dps_dpz = dE_rel_dpz * (E_rel + k_E * x) / ps

  sinh_c = sinh(coef)
  cosh_c = cosh(coef)

  if (abs(coef) < 1d-3) then
    sinh_k = alpha * (1 + coef**2 / 6 + coef**4/120)
    cosh1_k = alpha * coef * (1.0_rp / 2 + coef**2 / 24 + coef**4 / 720)
  else
    sinh_k = sinh_c / k_E
    cosh1_k = (cosh_c - 1) / k_E
  endif

  ff = -x * coef * sinh_c / ps - E_rel * length * sinh_c / ps**2 - px * length * cosh_c / ps**2
  mat6(1,1) = ff * dps_dx + cosh_c
  mat6(1,2) = ff * dps_dpx + sinh_k
  mat6(1,4) = ff * dps_dpy
  mat6(1,6) = ff * dps_dpz + dE_rel_dpz * cosh1_k

  ff = -(k_E * x + E_rel) * coef * cosh_c / ps - px * coef * sinh_c / ps
  mat6(2,1) = ff * dps_dx + k_E * sinh_c
  mat6(2,2) = ff * dps_dpx + cosh_c
  mat6(2,4) = ff * dps_dpy
  mat6(2,6) = ff * dps_dpz + dE_rel_dpz * sinh_c

  ff = -length * c00%vec(4) / ps**2
  mat6(3,1) = ff * dps_dx
  mat6(3,2) = ff * dps_dpx
  mat6(3,4) = ff * dps_dpy + length / ps
  mat6(3,6) = ff * dps_dpz

  ff = -x * coef * cosh_c / ps - E_rel * length * cosh_c / ps**2 - px * length * sinh_c / ps**2
  dbeta_dpz = mc2**2 * ele%value(p0c$) / E_tot**3
  beta_ref = ele%value(p0c$) / ele%value(e_tot$)
  mat6(5,1) = -c00%beta * (ff * dps_dx + sinh_c)
  mat6(5,2) = -c00%beta * (ff * dps_dpx + cosh1_k)
  mat6(5,4) = -c00%beta * (ff * dps_dpy)
  mat6(5,6) = dbeta_dpz * (length / beta_ref - (x * sinh_c + E_rel * sinh_k + px * cosh1_k)) - &
              c00%beta * (ff * dps_dpz + dE_rel_dpz * sinh_k)

  if (v(tilt_tot$) + angle /= 0) then
    call tilt_mat6 (mat6, v(tilt_tot$) + angle)
  endif

  call add_multipoles_and_z_offset (.true.)
  ele%vec0 = orb_out%vec - matmul(mat6, orb_in%vec)

!--------------------------------------------------------
! Kicker, etc.

case (kicker$, hkicker$, vkicker$, rcollimator$, ecollimator$, monitor$, instrument$, pipe$)

  call offset_particle (ele, param, set$, c00, set_hvkicks = .false.)

  charge_dir = rel_tracking_charge * ele%orientation

  hkick = charge_dir * v(hkick$) 
  vkick = charge_dir * v(vkick$) 
  kick  = charge_dir * v(kick$) 
  
  n_slice = max(1, nint(length / v(ds_step$)))
  if (ele%key == hkicker$) then
     c00%vec(2) = c00%vec(2) + kick / (2 * n_slice)
  elseif (ele%key == vkicker$) then
     c00%vec(4) = c00%vec(4) + kick / (2 * n_slice)
  else
     c00%vec(2) = c00%vec(2) + hkick / (2 * n_slice)
     c00%vec(4) = c00%vec(4) + vkick / (2 * n_slice)
  endif

  do i = 1, n_slice 
     call track_a_drift (c00, length/n_slice)
     call drift_mat6_calc (drift, length/n_slice, ele, param, c00)
     mat6 = matmul(drift,mat6)
     if (i == n_slice) then
        if (ele%key == hkicker$) then
           c00%vec(2) = c00%vec(2) + kick / (2 * n_slice)
        elseif (ele%key == vkicker$) then
           c00%vec(4) = c00%vec(4) + kick / (2 * n_slice)
        else
           c00%vec(2) = c00%vec(2) + hkick / (2 * n_slice)
           c00%vec(4) = c00%vec(4) + vkick / (2 * n_slice)
        endif
     else 
        if (ele%key == hkicker$) then
           c00%vec(2) = c00%vec(2) + kick / n_slice
        elseif (ele%key == vkicker$) then
           c00%vec(4) = c00%vec(4) + kick / n_slice
        else
           c00%vec(2) = c00%vec(2) + hkick / n_slice
           c00%vec(4) = c00%vec(4) + vkick / n_slice
        endif
     endif
  end do

  if (v(tilt_tot$) /= 0) then
    call tilt_mat6 (mat6, v(tilt_tot$))
  endif

  call add_multipoles_and_z_offset (.true.)
  ele%vec0 = orb_out%vec - matmul(mat6, orb_in%vec)

!--------------------------------------------------------
! LCavity: Linac rf cavity.
! Modified version of the ultra-relativistic formalism from:
!       J. Rosenzweig and L. Serafini
!       Phys Rev E, Vol. 49, p. 1599, (1994)
! with b_0 = b_-1 = 1. See the Bmad manual for more details.
!
! One must keep in mind that we are NOT using good canonical coordinates since
!   the energy of the reference particle is changing.
! This means that the resulting matrix will NOT be symplectic.

case (lcavity$)

  if (length == 0) return

  !

  call offset_particle (ele, param, set$, c00)

  phase = twopi * (v(phi0$) + v(phi0_multipass$) + &
                   v(phi0_autoscale$) +  v(phi0_err$) + &
                   (particle_ref_time (c00, ele) - rf_ref_time_offset(ele)) * v(rf_frequency$))

  ! Coupler kick

  if (v(coupler_strength$) /= 0) call mat6_coupler_kick(ele, param, first_track_edge$, phase, c00, mat6)

  ! 

  cos_phi = cos(phase)
  sin_phi = sin(phase)
  gradient_max = rel_tracking_charge * e_accel_field (ele, gradient$)
  gradient_net = gradient_max * cos_phi + gradient_shift_sr_wake(ele, param) 
  dE = gradient_net * length

  mc2 = mass_of(orb_in%species)
  pc_start_ref = v(p0c_start$) 
  pc_start = pc_start_ref * (1 + c00%vec(6))
  beta_start = c00%beta
  E_start = pc_start / beta_start
  E_end = E_start + dE
  if (E_end <= 0) then
    if (present(err)) err = .true.
    call out_io (s_error$, r_name, 'END ENERGY IS NEGATIVE AT ELEMENT: ' // ele%name)
    mat6 = 0   ! garbage.
    return 
  endif

  pc_end_ref = v(p0c$)
  call convert_total_energy_to (E_end, orb_in%species, pc = pc_end, beta = beta_end)
  E_end = pc_end / beta_end
  E_ratio = E_end / E_start

  om = twopi * v(rf_frequency$) / c_light
  om_g = om * gradient_max * length
  dbeta1_dE1 = mc2**2 / (pc_start * E_start**2)
  dbeta2_dE2 = mc2**2 / (pc_end * E_end**2)

  ! First convert from (x, px, y, py, z, pz) to (x, x', y, y', c(t_ref-t), E) coords 

  rel_p = 1 + c00%vec(6)
  mat6(2,:) = mat6(2,:) / rel_p - c00%vec(2) * mat6(6,:) / rel_p**2
  mat6(4,:) = mat6(4,:) / rel_p - c00%vec(4) * mat6(6,:) / rel_p**2

  m2(1,:) = [1/c00%beta, -c00%vec(5) * mc2**2 * c00%p0c / (pc_start**2 * E_start)]
  m2(2,:) = [0.0_rp, c00%p0c * c00%beta]
  mat6(5:6,:) = matmul(m2, mat6(5:6,:))

  c00%vec(2) = c00%vec(2) / rel_p
  c00%vec(4) = c00%vec(4) / rel_p
  c00%vec(5) = c00%vec(5) / c00%beta 
  c00%vec(6) = rel_p * c00%p0c / c00%beta - 1

  ! Body tracking longitudinal

  kmat6 = 0
  kmat6(6,5) = om_g * sin_phi
  kmat6(6,6) = 1

  if (abs(dE) <  1d-4*(pc_end+pc_start)) then
    dp_dg = length * (1 / beta_start - mc2**2 * dE / (2 * pc_start**3) + (mc2 * dE)**2 * E_start / (2 * pc_start**5))
    kmat6(5,5) = 1 - length * (-mc2**2 * kmat6(6,5) / (2 * pc_start**3) + mc2**2 * dE * kmat6(6,5) * E_start / pc_start**5)
    kmat6(5,6) = -length * (-dbeta1_dE1 / beta_start**2 + 2 * mc2**2 * dE / pc_start**4 + &
                    (mc2 * dE)**2 / (2 * pc_start**5) - 5 * (mc2 * dE)**2 / (2 * pc_start**5))
  else
    dp_dg = (pc_end - pc_start) / gradient_net
    kmat6(5,5) = 1 - kmat6(6,5) / (beta_end * gradient_net) + kmat6(6,5) * (pc_end - pc_start) / (gradient_net**2 * length)
    kmat6(5,6) = -1 / (beta_end * gradient_net) + 1 / (beta_start * gradient_net)
  endif

  ! Body tracking transverse

  if (nint(ele%value(cavity_type$)) == traveling_wave$) then

    kmat6(1,1) = 1
    kmat6(1,2) = length
    kmat6(2,2) = 1

    kmat6(3,3) = 1
    kmat6(3,4) = length
    kmat6(4,4) = 1

    kmat6(5,2) = -length * c00%vec(4)
    kmat6(5,4) = -length * c00%vec(4)

    c00%vec(5) = c00%vec(5) - (c00%vec(2)**2 + c00%vec(4)**2) * dp_dg / 2

    mat6 = matmul(kmat6, mat6)

    c00%vec(1:2) = matmul(kmat6(1:2,1:2), c00%vec(1:2))
    c00%vec(3:4) = matmul(kmat6(3:4,3:4), c00%vec(3:4))
    c00%vec(5) = c00%vec(5) - (dp_dg - c_light * v(delta_ref_time$))

  else
    sqrt_8 = 2 * sqrt_2
    voltage_max = gradient_max * length

    if (abs(voltage_max * cos_phi) < 1d-5 * E_start) then
      g = voltage_max / E_start
      alpha = g * (1 + g * cos_phi / 2)  / sqrt_8
      coef = length * beta_start * (1 - voltage_max * cos_phi / (2 * E_start))
      dalpha_dt1 = g * g * om * sin_phi / (2 * sqrt_8) 
      dalpha_dE1 = -(voltage_max / E_start**2 + voltage_max**2 * cos_phi / E_start**3) / sqrt_8
      dcoef_dt1 = -length * beta_start * sin_phi * om_g / (2 * E_start)
      dcoef_dE1 = length * beta_start * voltage_max * cos_phi / (2 * E_start**2) + coef * dbeta1_dE1 / beta_start
    else
      alpha = log(E_ratio) / (sqrt_8 * cos_phi)
      coef = sqrt_8 * pc_start * sin(alpha) / gradient_max
      dalpha_dt1 = kmat6(6,5) / (E_end * sqrt_8 * cos_phi) - log(E_ratio) * om * sin_phi / (sqrt_8 * cos_phi**2)
      dalpha_dE1 = 1 / (E_end * sqrt_8 * cos_phi) - 1 / (E_start * sqrt_8 * cos_phi)
      dcoef_dt1 = sqrt_8 * pc_start * cos(alpha) * dalpha_dt1 / gradient_max
      dcoef_dE1 = coef / (beta_start * pc_start) + sqrt_8 * pc_start * cos(alpha) * dalpha_dE1 / gradient_max
    endif

    cos_a = cos(alpha)
    sin_a = sin(alpha)

    z21 = -gradient_max / (sqrt_8 * pc_end)
    z22 = pc_start / pc_end  

    c_min = cos_a - sqrt_2 * beta_start * sin_a * cos_phi
    c_plu = cos_a + sqrt_2 * beta_end * sin_a * cos_phi
    dc_min = -sin_a - sqrt_2 * beta_start * cos_a * cos_phi 
    dc_plu = -sin_a + sqrt_2 * beta_end * cos_a * cos_phi 

    cos_term = 1 + 2 * beta_start * beta_end * cos_phi**2
    dcos_phi = om * sin_phi

    kmat6(1,1) =  c_min
    kmat6(1,2) =  coef 
    kmat6(2,1) =  z21 * (sqrt_2 * (beta_start - beta_end) * cos_phi * cos_a + cos_term * sin_a)
    kmat6(2,2) =  c_plu * z22

    kmat6(1,5) = c00%vec(1) * (dc_min * dalpha_dt1 - sqrt_2 * beta_start * sin_a * dcos_phi) + & 
                 c00%vec(2) * (dcoef_dt1)

    kmat6(1,6) = c00%vec(1) * (dc_min * dalpha_dE1 - sqrt_2 * dbeta1_dE1 * sin_a * cos_phi) + &
                 c00%vec(2) * (dcoef_dE1)

    kmat6(2,5) = c00%vec(1) * z21 * (sqrt_2 * (beta_start - beta_end) * (dcos_phi * cos_a - cos_phi * sin_a * dalpha_dt1)) + &
                 c00%vec(1) * z21 * sqrt_2 * (-dbeta2_dE2 * kmat6(6,5)) * cos_phi * cos_a + &
                 c00%vec(1) * z21 * (4 * beta_start * beta_end *cos_phi * dcos_phi * sin_a + cos_term * cos_a * dalpha_dt1) + &
                 c00%vec(1) * z21 * (2 * beta_start * dbeta2_dE2 * kmat6(6,5) * sin_a) + &
                 c00%vec(1) * (-kmat6(2,1) * kmat6(6,5) / (beta_end * pc_end)) + &
                 c00%vec(2) * z22 * (dc_plu * dalpha_dt1 + sqrt_2 * sin_a * (beta_end * dcos_phi + dbeta2_dE2 * kmat6(6,5) * cos_phi)) + &
                 c00%vec(2) * z22 * (-c_plu * kmat6(6,5) / (beta_end * pc_end))

    kmat6(2,6) = c00%vec(1) * z21 * (sqrt_2 * cos_phi * ((dbeta1_dE1 - dbeta2_dE2) * cos_a - (beta_start - beta_end) * sin_a * dalpha_dE1)) + &
                 c00%vec(1) * z21 * (2 * cos_phi**2 * (dbeta1_dE1 * beta_end + beta_start * dbeta2_dE2) * sin_a + cos_term * cos_a * dalpha_dE1) + &
                 c00%vec(1) * (-kmat6(2,1) / (beta_end * pc_end)) + &
                 c00%vec(2) * z22 * (dc_plu * dalpha_dE1 + sqrt_2 * dbeta2_dE2 * sin_a * cos_phi) + &
                 c00%vec(2) * c_plu * (1 / (beta_start * pc_end) - pc_start / (beta_end * pc_end**2))

    kmat6(3:4,3:4) = kmat6(1:2,1:2)

    kmat6(3,5) = c00%vec(3) * (dc_min * dalpha_dt1 - sqrt_2 * beta_start * sin_a * dcos_phi) + & 
                 c00%vec(4) * (dcoef_dt1)

    kmat6(3,6) = c00%vec(3) * (dc_min * dalpha_dE1 - sqrt_2 * dbeta1_dE1 * sin_a * cos_phi) + &
                 c00%vec(4) * (dcoef_dE1)

    kmat6(4,5) = c00%vec(3) * z21 * (sqrt_2 * (beta_start - beta_end) * (dcos_phi * cos_a - cos_phi * sin_a * dalpha_dt1)) + &
                 c00%vec(3) * z21 * sqrt_2 * (-dbeta2_dE2 * kmat6(6,5)) * cos_phi * cos_a + &
                 c00%vec(3) * z21 * (4 * beta_start * beta_end *cos_phi * dcos_phi * sin_a + cos_term * cos_a * dalpha_dt1) + &
                 c00%vec(3) * z21 * (2 * beta_start * dbeta2_dE2 * kmat6(6,5) * sin_a) + &
                 c00%vec(3) * (-kmat6(2,1) * kmat6(6,5) / (beta_end * pc_end)) + &
                 c00%vec(4) * z22 * (dc_plu * dalpha_dt1 + sqrt_2 * sin_a * (beta_end * dcos_phi + dbeta2_dE2 * kmat6(6,5) * cos_phi)) + &
                 c00%vec(4) * z22 * (-c_plu * kmat6(6,5) / (beta_end * pc_end))

    kmat6(4,6) = c00%vec(3) * z21 * (sqrt_2 * cos_phi * ((dbeta1_dE1 - dbeta2_dE2) * cos_a - (beta_start - beta_end) * sin_a * dalpha_dE1)) + &
                 c00%vec(3) * z21 * (2 * cos_phi**2 * (dbeta1_dE1 * beta_end + beta_start * dbeta2_dE2) * sin_a + cos_term * cos_a * dalpha_dE1) + &
                 c00%vec(3) * (-kmat6(2,1) / (beta_end * pc_end)) + &
                 c00%vec(4) * z22 * (dc_plu * dalpha_dE1 + sqrt_2 * dbeta2_dE2 * sin_a * cos_phi) + &
                 c00%vec(4) * c_plu * (1 / (beta_start * pc_end) - pc_start / (beta_end * pc_end**2))


    ! Correction to z for finite x', y'
    ! Note: Corrections to kmat6(5,5) and kmat6(5,6) are ignored since these are small (quadratic
    ! in the transvers coords).

    c_plu = sqrt_2 * cos_phi * cos_a + sin_a

    drp1_dr0  = -gradient_net / (2 * E_start)
    drp1_drp0 = 1

    xp1 = drp1_dr0 * c00%vec(1) + drp1_drp0 * c00%vec(2)
    yp1 = drp1_dr0 * c00%vec(3) + drp1_drp0 * c00%vec(4)

    drp2_dr0  = (c_plu * z21)
    drp2_drp0 = (cos_a * z22)

    xp2 = drp2_dr0 * c00%vec(1) + drp2_drp0 * c00%vec(2)
    yp2 = drp2_dr0 * c00%vec(3) + drp2_drp0 * c00%vec(4)

    kmat6(5,1) = -(c00%vec(1) * (drp1_dr0**2 + drp1_dr0*drp2_dr0 + drp2_dr0**2) + &
                   c00%vec(2) * (drp1_dr0 * drp1_drp0 + drp2_dr0 * drp2_drp0 + &
                                (drp1_dr0 * drp2_drp0 + drp1_drp0 * drp2_dr0) / 2)) * dp_dg / 3

    kmat6(5,2) = -(c00%vec(2) * (drp1_drp0**2 + drp1_drp0*drp2_drp0 + drp2_drp0**2) + &
                   c00%vec(1) * (drp1_dr0 * drp1_drp0 + drp2_dr0 * drp2_drp0 + &
                                (drp1_dr0 * drp2_drp0 + drp1_drp0 * drp2_dr0) / 2)) * dp_dg / 3

    kmat6(5,3) = -(c00%vec(3) * (drp1_dr0**2 + drp1_dr0*drp2_dr0 + drp2_dr0**2) + &
                   c00%vec(4) * (drp1_dr0 * drp1_drp0 + drp2_dr0 * drp2_drp0 + &
                                (drp1_dr0 * drp2_drp0 + drp1_drp0 * drp2_dr0) / 2)) * dp_dg / 3

    kmat6(5,4) = -(c00%vec(4) * (drp1_drp0**2 + drp1_drp0*drp2_drp0 + drp2_drp0**2) + &
                   c00%vec(3) * (drp1_dr0 * drp1_drp0 + drp2_dr0 * drp2_drp0 + &
                                (drp1_dr0 * drp2_drp0 + drp1_drp0 * drp2_dr0) / 2)) * dp_dg / 3

    c00%vec(5) = c00%vec(5) - (xp1**2 + xp1*xp2 + xp2**2 + yp1**2 + yp1*yp2 + yp2**2) * dp_dg / 6

    !

    mat6 = matmul(kmat6, mat6)

    c00%vec(1:2) = matmul(kmat6(1:2,1:2), c00%vec(1:2))
    c00%vec(3:4) = matmul(kmat6(3:4,3:4), c00%vec(3:4))
    c00%vec(5) = c00%vec(5) - (dp_dg - c_light * v(delta_ref_time$))

  endif

  ! Convert back from (x, x', y, y', c(t-t_ref), E)  to (x, px, y, py, z, pz) coords
  ! Here the effective t used in calculating m2 is zero so m2(1,2) is zero.

  rel_p = pc_end / pc_end_ref
  mat6(2,:) = rel_p * mat6(2,:) + c00%vec(2) * mat6(6,:) / (pc_end_ref * beta_end)
  mat6(4,:) = rel_p * mat6(4,:) + c00%vec(4) * mat6(6,:) / (pc_end_ref * beta_end)

  m2(1,:) = [beta_end, c00%vec(5) * mc2**2 / (pc_end * E_end**2)]
  m2(2,:) = [0.0_rp, 1 / (pc_end_ref * beta_end)]

  mat6(5:6,:) = matmul(m2, mat6(5:6,:))

  c00%vec(2) = c00%vec(2) / rel_p
  c00%vec(4) = c00%vec(4) / rel_p
  c00%vec(6) = (pc_end - pc_end_ref) / pc_end_ref 
  c00%p0c = pc_end_ref
  c00%beta = beta_end

  ! Coupler kick

  if (v(coupler_strength$) /= 0) call mat6_coupler_kick(ele, param, second_track_edge$, phase, c00, mat6)

  ! multipoles and z_offset

  if (v(tilt_tot$) /= 0) call tilt_mat6 (mat6, v(tilt_tot$))

  call add_multipoles_and_z_offset (.true.)
  ele%vec0 = orb_out%vec - matmul(mat6, orb_in%vec)

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
! Mirror

case (mirror$)

  mat6(1, 1) = -1
  mat6(2, 1) =  0   ! 
  mat6(2, 2) = -1
  mat6(4, 3) =  0

  if (ele%photon%surface%has_curvature) then
    print *, 'MIRROR CURVATURE NOT YET IMPLEMENTED!'
    call err_exit
  endif

  ! Offsets?

  ele%vec0 = orb_out%vec - matmul(mat6, orb_in%vec)

!--------------------------------------------------------
! Multipole, AB_Multipole

case (multipole$, ab_multipole$)

  if (.not. ele%multipoles_on) return

  call offset_particle (ele, param, set$, c00, set_tilt = .false.)

  call multipole_ele_to_kt (ele, .true., has_nonzero_pole, knl, tilt)
  if (has_nonzero_pole) then
    call multipole_kick_mat (knl*charge_dir, tilt, c00, 1.0_rp, ele%mat6)

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
! Patch

case (patch$) 

  if (ele%field_calc == custom$) then
    call out_io (s_fatal$, r_name, 'MAT6_CALC_METHOD=BMAD_STANDARD CANNOT HANDLE FIELD_CALC=CUSTOM', &
                                   'FOR PATCH ELEMENT: ' // ele%name)
    if (global_com%exit_on_error) call err_exit
    return
  endif

  mc2 = mass_of(orb_in%species)
  c00%vec(5) = 0
  call track_a_patch (ele, c00, .false., s_ent, ds_ref, ww, p_s, track_spin = .false.)

  ! Coordinate transform before drift

  rel_p = 1 + orb_in%vec(6)
  if (ele%orientation*orb_in%direction == 1) then
    dir = ele%value(upstream_ele_dir$) * orb_in%direction
  else
    dir = ele%value(downstream_ele_dir$) * orb_in%direction
  endif
  pz = sqrt(rel_p**2 - orb_in%vec(2)**2 - orb_in%vec(4)**2) * dir
  beta_ref = v(p0c$) / v(e_tot$)

  dps_dpx = ww(3,1) - ww(3,3) * orb_in%vec(2) / pz
  dps_dpy = ww(3,2) - ww(3,3) * orb_in%vec(4) / pz
  dps_dpz = ww(3,3) * rel_p / pz

  mat6(1,1) = ww(1,1) - ww(3,1) * c00%vec(2) / p_s
  mat6(1,2) = -s_ent * ((ww(1,1) - ww(1,3) * orb_in%vec(2) / pz) / p_s - c00%vec(2) * dps_dpx / p_s**2) 
  mat6(1,3) = ww(1,2) - ww(3,2) * c00%vec(2) / p_s
  mat6(1,4) = -s_ent * ((ww(1,2) - ww(1,3) * orb_in%vec(4) / pz) / p_s - c00%vec(2) * dps_dpy / p_s**2)
  mat6(1,6) = -s_ent * (ww(1,3) * rel_p / (p_s * pz) - c00%vec(2) * dps_dpz / p_s**2)

  mat6(2,2) = ww(1,1) - ww(1,3) * orb_in%vec(2) / pz
  mat6(2,4) = ww(1,2) - ww(1,3) * orb_in%vec(4) / pz
  mat6(2,6) = ww(1,3) * rel_p / pz

  mat6(3,1) = ww(2,1) - ww(3,1) * c00%vec(4) / p_s
  mat6(3,2) = -s_ent * ((ww(2,1) - ww(2,3) * orb_in%vec(2) / pz) / p_s - c00%vec(4) * dps_dpx / p_s**2)
  mat6(3,3) = ww(2,2) - ww(3,2) * c00%vec(4) / p_s
  mat6(3,4) = -s_ent * ((ww(2,2) - ww(2,3) * orb_in%vec(4) / pz) / p_s - c00%vec(4) * dps_dpy / p_s**2) 
  mat6(3,6) = -s_ent * (ww(2,3) * rel_p / (p_s * pz) - c00%vec(4) * dps_dpz / p_s**2)

  mat6(4,2) = ww(2,1) - ww(2,3) * orb_in%vec(2) / pz
  mat6(4,4) = ww(2,2) - ww(2,3) * orb_in%vec(4) / pz
  mat6(4,6) = ww(2,3) * rel_p / pz

  mat6(5,1) = ww(3,1) * rel_p / p_s
  mat6(5,2) = -s_ent * rel_p * dps_dpx / p_s**2
  mat6(5,3) = ww(3,2) * rel_p / p_s
  mat6(5,4) = -s_ent * rel_p * dps_dpy / p_s**2
  mat6(5,6) = v(t_offset$) * c_light * mc2**2 * orb_in%beta**3 / (v(p0c_start$)**2 * rel_p**3) + &
              s_ent / p_s - s_ent * rel_p * dps_dpz / p_s**2 + &
              ds_ref * mc2**2 * c00%beta**3 / (rel_p**3 * v(p0c_start$)**2 * beta_ref)

  ! Energy offset

  if (v(p0c_start$) /= v(p0c$)) then
    dp_ratio = v(p0c_start$) / v(p0c$)
    mat6(2,:) = mat6(2,:) * dp_ratio
    mat6(4,:) = mat6(4,:) * dp_ratio
    mat6(6,:) = mat6(6,:) * dp_ratio
  endif

  ! veorb_in calc

  ele%vec0 = orb_out%vec - matmul(mat6, orb_in%vec)

!--------------------------------------------------------
! quadrupole

case (quadrupole$)

  call track_a_quadrupole (c00, ele, param, mat6, .true.)

!--------------------------------------------------------
! rbends are not allowed internally

case (rbend$)

  if (present(err)) err = .true.
  call out_io (s_fatal$, r_name,  'RBEND ELEMENTS NOT ALLOWED INTERNALLY!')
  if (global_com%exit_on_error) call err_exit
  return

!--------------------------------------------------------
! rf cavity
! Calculation Uses a 3rd order map assuming a linearized rf voltage vs time.

case (rfcavity$)

  mc2 = mass_of(orb_in%species)
  p0c = v(p0c$)
  beta_ref = p0c / v(e_tot$)
  n_slice = max(1, nint(length / v(ds_step$))) 
  dt_ref = length / (c_light * beta_ref)

  call offset_particle (ele, param, set$, c00, set_tilt = .false.)

  ! The cavity field is modeled as a standing wave antisymmetric wrt the center.
  ! Thus if the cavity is flipped (orientation = -1), the wave of interest, which is 
  ! always the accelerating wave, is the "backward" wave. And the phase of the backward 
  ! wave is different from the phase of the forward wave by a constant dt_ref * freq

  voltage = e_accel_field (ele, voltage$) * charge_dir

  phase0 = twopi * (v(phi0$) + v(phi0_multipass$) - v(phi0_autoscale$) - &
                  (particle_ref_time (c00, ele) - rf_ref_time_offset(ele)) * v(rf_frequency$))
  if (ele%orientation == -1) phase0 = phase0 + twopi * v(rf_frequency$) * dt_ref
  phase = phase0

  t0 = c00%t

  ! Track through slices.
  ! The phase of the accelerating wave traveling in the same direction as the particle is
  ! assumed to be traveling with a phase velocity the same speed as the reference velocity.

  if (v(coupler_strength$) /= 0) call mat6_coupler_kick(ele, param, first_track_edge$, phase, c00, mat6)

  do i = 0, n_slice

    factor = voltage / n_slice
    if (i == 0 .or. i == n_slice) factor = factor / 2

    dE = factor * sin(phase)
    pc = (1 + c00%vec(6)) * p0c 
    E = pc / c00%beta
    call convert_total_energy_to (E + dE, orb_in%species, pc = new_pc, beta = new_beta)
    ff = twopi * factor * v(rf_frequency$) * cos(phase) / (p0c * new_beta * c_light)

    m2(2,1) = ff / c00%beta
    m2(2,2) = c00%beta / new_beta - ff * c00%vec(5) *mc2**2 * p0c / (E * pc**2) 
    m2(1,1) = new_beta / c00%beta + c00%vec(5) * (mc2**2 * p0c * m2(2,1) / (E+dE)**3) / c00%beta
    m2(1,2) = c00%vec(5) * mc2**2 * p0c * (m2(2,2) / ((E+dE)**3 * c00%beta) - new_beta / (pc**2 * E))

    mat6(5:6, :) = matmul(m2, mat6(5:6, :))
  
    c00%vec(6) = (new_pc - p0c) / p0c
    c00%vec(5) = c00%vec(5) * new_beta / c00%beta
    c00%beta   = new_beta

    if (i /= n_slice) then
      call drift_mat6_calc (drift, length/n_slice, ele, param, c00)
      call track_a_drift (c00, length/n_slice)
      mat6 = matmul(drift, mat6)
      phase = phase0 + twopi * v(rf_frequency$) * ((i + 1) * dt_ref/n_slice - (c00%t - t0)) 
    endif

  enddo

  ! Coupler kick

  if (v(coupler_strength$) /= 0) call mat6_coupler_kick(ele, param, second_track_edge$, phase, c00, mat6)

  call offset_particle (ele, param, unset$, c00, set_tilt = .false.)

  !

  if (v(tilt_tot$) /= 0) call tilt_mat6 (mat6, v(tilt_tot$))

  call add_multipoles_and_z_offset (.true.)
  ele%vec0 = orb_out%vec - matmul(mat6, orb_in%vec)

!--------------------------------------------------------
! sad_mult

case (sad_mult$)

  call sad_mult_track_and_mat (ele, param, c00, orb_out1, ele%mat6)
  if (.not. logic_option (.false., end_in)) call set_orb_out (orb_out, c00)

!--------------------------------------------------------
! sbend

case (sbend$)

  call track_a_bend (c00, ele, param, mat6, .true.)

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
  call solenoid_track_and_mat (ele, param, c00, c00, mat6)
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

subroutine add_multipole_slice (knl, tilt, factor, orb, mat6)

type (coord_struct) orb
real(rp) knl(0:n_pole_maxx), tilt(0:n_pole_maxx)
real(rp) mat6(6,6), factor, mat6_m(6,6)

!

call multipole_kick_mat (knl, tilt, orb, factor, mat6_m)

mat6(2,:) = mat6(2,:) + mat6_m(2,1) * mat6(1,:) + mat6_m(2,3) * mat6(3,:)
mat6(4,:) = mat6(4,:) + mat6_m(4,1) * mat6(1,:) + mat6_m(4,3) * mat6(3,:)

call multipole_kicks (knl*factor, tilt, orb)

end subroutine

!--------------------------------------------------------
! contains

subroutine set_orb_out (orb_out, c00)

type (coord_struct) orb_out, c00
orb_out = c00
orb_out%location = downstream_end$

end subroutine set_orb_out

!--------------------------------------------------------
! contains

! put in multipole components

subroutine add_multipoles_and_z_offset (add_pole)

real(rp) mat6_m(6,6)
logical has_nonzero_pole, add_pole

!

if (add_pole) then
  call multipole_ele_to_kt (ele, .true., has_nonzero_pole, knl, tilt)
  if (has_nonzero_pole) then
    knl = knl * charge_dir
    call multipole_kick_mat (knl, tilt, orb_in, 0.5_rp, mat6_m)
    mat6(:,1) = mat6(:,1) + mat6(:,2) * mat6_m(2,1) + mat6(:,4) * mat6_m(4,1)
    mat6(:,3) = mat6(:,3) + mat6(:,2) * mat6_m(2,3) + mat6(:,4) * mat6_m(4,3)
    call multipole_kick_mat (knl, tilt, orb_out, 0.5_rp, mat6_m)
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

!----------------------------------------------------------------
!----------------------------------------------------------------
!----------------------------------------------------------------

subroutine mat6_coupler_kick(ele, param, particle_at, phase, orb, mat6)

use track1_mod

implicit none

type (ele_struct) ele
type (coord_struct) orb, old_orb
type (lat_param_struct) param
real(rp) phase, mat6(6,6), f, f2, coef, E_new
real(rp) dp_coef, dp_x, dp_y, ph, mc(6,6), E, pc, mc2, p0c
integer particle_at, physical_end

!

physical_end = physical_ele_end(particle_at, orb%direction, ele%orientation)
if (.not. at_this_ele_end (physical_end, nint(ele%value(coupler_at$)))) return

ph = phase
if (ele%key == rfcavity$) ph = pi/2 - ph
ph = ph + twopi * ele%value(coupler_phase$)

mc2 = mass_of(orb%species)
p0c = orb%p0c
pc = p0c * (1 + orb%vec(6))
E = pc / orb%beta

f = twopi * ele%value(rf_frequency$) / c_light
dp_coef = e_accel_field(ele, gradient$) * ele%value(coupler_strength$)
dp_x = dp_coef * cos(twopi * ele%value(coupler_angle$))
dp_y = dp_coef * sin(twopi * ele%value(coupler_angle$))

if (nint(ele%value(coupler_at$)) == both_ends$) then
  dp_x = dp_x / 2
  dp_y = dp_y / 2
endif

! Track

old_orb = orb
call rf_coupler_kick (ele, param, particle_at, phase, orb)

! Matrix

call mat_make_unit (mc)

mc(2,5) = dp_x * f * sin(ph) / (old_orb%beta * p0c)
mc(4,5) = dp_y * f * sin(ph) / (old_orb%beta * p0c)

mc(2,6) = -dp_x * f * sin(ph) * old_orb%vec(5) * mc2**2 / (E * pc**2)
mc(4,6) = -dp_y * f * sin(ph) * old_orb%vec(5) * mc2**2 / (E * pc**2)

coef = (dp_x * old_orb%vec(1) + dp_y * old_orb%vec(3)) * cos(ph) * f**2 
mc(6,1) = dp_x * sin(ph) * f / (orb%beta * p0c)
mc(6,3) = dp_y * sin(ph) * f / (orb%beta * p0c)
mc(6,5) = -coef / (orb%beta * old_orb%beta * p0c) 
mc(6,6) = old_orb%beta/orb%beta + coef * old_orb%vec(5) * mc2**2 / (pc**2 * E * orb%beta)

f2 = old_orb%vec(5) * mc2**2 / (pc * E**2 * p0c)
E_new = p0c * (1 + orb%vec(6)) / orb%beta

mc(5,1) = old_orb%vec(5) * mc2**2 * p0c * mc(6,1) / (old_orb%beta * E_new**3)
mc(5,3) = old_orb%vec(5) * mc2**2 * p0c * mc(6,3) / (old_orb%beta * E_new**3)
mc(5,5) = orb%beta/old_orb%beta + old_orb%vec(5) * mc2**2 * p0c * mc(6,5) / (old_orb%beta * E_new**3)
mc(5,6) = old_orb%vec(5) * mc2**2 * p0c * (mc(6,6) / (old_orb%beta * E_new**3) - &
                                     orb%beta / (old_orb%beta**2 * E**3))

mat6 = matmul(mc, mat6)

end subroutine mat6_coupler_kick

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

subroutine lcavity_edge_kick_matrix (ele, param, grad_max, phase, orb, mat6)

use bmad_interface

implicit none

type (ele_struct)  ele
type (coord_struct)  orb
type (lat_param_struct) param

real(rp) grad_max, phase, k1, mat6(6,6)
real(rp) f, mc2, E, pc

! Note that phase space here is (x, x', y, y', -c(t-t_ref), E) 

pc = (1 + orb%vec(6)) * orb%p0c
E = pc / orb%beta
k1 = grad_max * cos(phase)
f = grad_max * sin(phase) * twopi * ele%value(rf_frequency$) / c_light
mc2 = mass_of(orb%species)

mat6(2,:) = mat6(2,:) + k1 * mat6(1,:) + f * orb%vec(1) * mat6(5,:) 
mat6(4,:) = mat6(4,:) + k1 * mat6(3,:) - f * orb%vec(3) * mat6(5,:)

orb%vec(2) = orb%vec(2) + k1 * orb%vec(1)
orb%vec(4) = orb%vec(4) + k1 * orb%vec(3)

end subroutine lcavity_edge_kick_matrix

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

subroutine bbi_kick_matrix (ele, param, orb, s_pos, mat6)

use bmad_interface, except_dummy => bbi_kick_matrix

implicit none

type (ele_struct)  ele
type (coord_struct)  orb
type (lat_param_struct) param

real(rp) x_pos, y_pos, del, sig_x, sig_y, coef, garbage, s_pos
real(rp) ratio, k0_x, k1_x, k0_y, k1_y, mat6(6,6), beta, bbi_const
real(rp) beta_a0, alpha_a0, beta_b0, alpha_b0

!

call mat_make_unit (mat6)

sig_x = ele%value(sig_x$)
sig_y = ele%value(sig_y$)

if (sig_x == 0 .or. sig_y == 0) return

if (ele%value(beta_a$) == 0) then
  beta_a0 = ele%a%beta
  alpha_a0 = ele%a%alpha
else
  beta_a0 = ele%value(beta_a$)
  alpha_a0 = ele%value(alpha_a$)
endif

if (ele%value(beta_b$) == 0) then
  beta_b0 = ele%b%beta
  alpha_b0 = ele%b%alpha
else
  beta_b0 = ele%value(beta_b$)
  alpha_b0 = ele%value(alpha_b$)
endif

if (s_pos /= 0 .and. ele%a%beta /= 0) then
  beta = beta_a0 - 2 * alpha_a0 * s_pos + (1 + alpha_a0**2) * s_pos**2 / beta_a0
  sig_x = sig_x * sqrt(beta / beta_a0)
  beta = beta_b0 - 2 * alpha_b0 * s_pos + (1 + alpha_b0**2) * s_pos**2 / beta_b0
  sig_y = sig_y * sqrt(beta / beta_b0)
endif

x_pos = orb%vec(1) / sig_x  ! this has offset in it
y_pos = orb%vec(3) / sig_y

del = 0.001

ratio = sig_y / sig_x
call bbi_kick (x_pos, y_pos, ratio, k0_x, k0_y)
call bbi_kick (x_pos+del, y_pos, ratio, k1_x, garbage)
call bbi_kick (x_pos, y_pos+del, ratio, garbage, k1_y)

bbi_const = -param%n_part * ele%value(charge$) * classical_radius_factor /  &
                    (2 * pi * ele%value(p0c$) * (sig_x + sig_y))

coef = bbi_const / (ele%value(n_slice$) * del * (1 + orb%vec(6)))

mat6(2,1) = coef * (k1_x - k0_x) / sig_x
mat6(4,3) = coef * (k1_y - k0_y) / sig_y

end subroutine bbi_kick_matrix
