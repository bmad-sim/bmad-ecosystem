!+
! Subroutine track1_bmad (start_orb, ele, param, end_orb, err_flag)
!
! Particle tracking through a single element BMAD_standard style.
!
! Note: track1_bmad *never* relies on ele%mat6 for tracking excect for 
! hybrid elements.
! 
! Modules Needed:
!   use bmad
!
! Input:
!   start_orb  -- Coord_struct: Starting position
!   ele        -- Ele_struct: Element
!   param      -- lat_param_struct:
!     %particle     -- Particle type
!
! Output:
!   end_orb   -- Coord_struct: End position
!   err_flag  -- Logical, optional: Set true if there is an error. False otherwise.
!-

subroutine track1_bmad (start_orb, ele, param, end_orb, err_flag)

use sad_mod, dummy1 => track1_bmad
use track1_mod, dummy2 => track1_bmad
use mad_mod, dummy3 => track1_bmad
use geometry_mod, dummy4 => track1_bmad
use ptc_interface_mod, dummy5 => track1_bmad

implicit none

type (coord_struct) :: start_orb, start2_orb
type (coord_struct) :: end_orb, temp_orb
type (ele_struct) :: ele, temp_ele
type (ele_struct), pointer :: ele0
type (lat_param_struct) :: param
type (taylor_struct) taylor1(6), taylor2(6)

real(rp) k1, k2, k2l, k3l, length, phase0, phase, beta_start, beta_ref
real(rp) beta_end, beta_start_ref, beta_end_ref, hkick, vkick, kick
real(rp) e2, sig_x, sig_y, kx, ky, coef, bbi_const, voltage
real(rp) knl(0:n_pole_maxx), tilt(0:n_pole_maxx)
real(rp) ks, kss, ksr, sig_x0, sig_y0, beta, mat6(6,6), mat2(2,2), mat4(4,4)
real(rp) z_slice(100), s_pos, s_pos_old, vec0(6)
real(rp) rel_p, k_z, pc_start, pc_end, dt_ref, gradient_ref, gradient_max
real(rp) x_pos, y_pos, cos_phi, gradient_net, e_start, e_end, e_ratio, voltage_max
real(rp) alpha, sin_a, cos_a, f, r_mat(2,2), volt_ref, p_bend, p2_bend, p2xy
real(rp) x, y, z, px, py, pz, k, dE0, L, E, pxy2, xp1, xp2, yp1, yp2, px_end
real(rp) xp_start, yp_start, dz4_coef(4,4), dz_coef(3), sqrt_8, pr, kk
real(rp) dcos_phi, dgradient, dpz, sin_alpha_over_f
real(rp) mc2, dpc_start, dE_start, dE_end, dE, dp_dg, dp_dg_ref, g
real(rp) E_start_ref, E_end_ref, pc_start_ref, pc_end_ref
real(rp) new_pc, new_beta, len_slice, k0l, k1l, t0
real(rp) cosh1_k, sinh_k, hk, vk, dt, k_E, e_rel, beta_a0, beta_b0, alpha_a0, alpha_b0
real(rp) p_factor, sin_alpha, cos_alpha, sin_psi, cos_psi, wavelength
real(rp) cos_g, sin_g, cos_tc, sin_tc, angle, rel_tracking_charge, rtc
real(rp) k_in_norm(3), h_norm(3), k_out_norm(3), e_tot, pc, ps
real(rp) cap_gamma, gamma_0, gamma_h, b_err, dtheta_sin_2theta, b_eff
real(rp) m_in(3,3) , m_out(3,3), y_out(3), x_out(3), k_out(3)
real(rp) test, nn, mm, temp_vec(3), p_vec(3), r_vec(3), charge_dir

integer i, n, n_slice, key, orientation, ix_sec

logical, optional :: err_flag
logical err, has_nonzero_pole

character(*), parameter :: r_name = 'track1_bmad'

type (mad_map_struct) map
type (mad_energy_struct) energy

! initially set end_orb = start_orb

if (present(err_flag)) err_flag = .false.

start2_orb = start_orb ! In case start_orb and end_orb share the same memory.

end_orb = start_orb     ! transfer start to end
if (end_orb%species /= photon$) then
  end_orb%p0c = ele%value(p0c$)
endif
length = ele%value(l$)
rel_p = 1 + start_orb%vec(6)
orientation = ele%orientation * start_orb%direction
rel_tracking_charge = relative_tracking_charge(start_orb, param)
charge_dir = rel_tracking_charge * orientation

!-----------------------------------------------
! If element is off... 

key = ele%key
if (key == sol_quad$ .and. ele%value(k1$) == 0) key = solenoid$

if (.not. ele%is_on) then
  select case (key)
  case (taylor$, match$, fiducial$, floor_shift$)
    call set_end_orb_s()
    return
  case (ab_multipole$, multipole$, lcavity$, sbend$, patch$)
    ! Note: LCavities will still do wakefields.
  case default
    key = drift$  
  end select
endif

! Select.

select case (key)

!-----------------------------------------------
! beambeam
                        
case (beambeam$)

  if (ele%value(charge$) == 0 .or. param%n_part == 0) return

  sig_x0 = ele%value(sig_x$)
  sig_y0 = ele%value(sig_y$)
  if (sig_x0 == 0 .or. sig_y0 == 0) return

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

  call offset_particle (ele, param, set$, end_orb)
  call canonical_to_angle_coords (end_orb)

  n_slice = max(1, nint(ele%value(n_slice$)))
  call bbi_slice_calc (ele, n_slice, z_slice)
  s_pos = 0    ! end at the ip
  do i = 1, n_slice
    s_pos_old = s_pos
    s_pos = (end_orb%vec(5) + z_slice(i)) / 2
    end_orb%vec(1) = end_orb%vec(1) + end_orb%vec(2) * (s_pos - s_pos_old)
    end_orb%vec(3) = end_orb%vec(3) + end_orb%vec(4) * (s_pos - s_pos_old)
    if (beta_a0 == 0) then
      sig_x = sig_x0
      sig_y = sig_y0
    else
      beta = beta_a0 - 2 * alpha_a0 * s_pos + (1 + alpha_a0**2) * s_pos**2 / beta_a0
      sig_x = sig_x0 * sqrt(beta / beta_a0)
      beta = beta_b0 - 2 * alpha_b0 * s_pos + (1 + alpha_b0**2) * s_pos**2 / beta_b0
      sig_y = sig_y0 * sqrt(beta / beta_b0)
    endif

    call bbi_kick (end_orb%vec(1)/sig_x, end_orb%vec(3)/sig_y, sig_y/sig_x,  kx, ky)
    bbi_const = -param%n_part * ele%value(charge$) * classical_radius_factor /  &
                                          (2 * pi * ele%value(p0c$) * (sig_x + sig_y))
    coef = ele%value(bbi_const$) / (n_slice * rel_p)
    end_orb%vec(2) = end_orb%vec(2) + kx * coef
    end_orb%vec(4) = end_orb%vec(4) + ky * coef
  enddo
  end_orb%vec(1) = end_orb%vec(1) - end_orb%vec(2) * s_pos
  end_orb%vec(3) = end_orb%vec(3) - end_orb%vec(4) * s_pos

  call angle_to_canonical_coords (end_orb)
  call offset_particle (ele, param, unset$, end_orb)  

!-----------------------------------------------
! collimator

case (rcollimator$, ecollimator$, monitor$, instrument$, pipe$) 

  call offset_particle (ele, param, set$, end_orb, set_tilt = .false., set_hvkicks = .false.)

  kick = sqrt(ele%value(vkick$)**2 + ele%value(hkick$)**2) * charge_dir

  if (length == 0) then
    end_orb%vec(2) = end_orb%vec(2) + ele%value(hkick$) * charge_dir
    end_orb%vec(4) = end_orb%vec(4) + ele%value(vkick$) * charge_dir

  elseif (kick == 0) then
    call track_a_drift (end_orb, length)

  else
    angle = atan2(ele%value(vkick$), ele%value(hkick$))
    cos_a = cos(angle); sin_a = sin(angle)

    end_orb%vec(1:3:2) = [end_orb%vec(1) * cos_a + end_orb%vec(3) * sin_a, &
                         -end_orb%vec(1) * sin_a + end_orb%vec(3) * cos_a]
    end_orb%vec(2:4:2) = [end_orb%vec(2) * cos_a + end_orb%vec(4) * sin_a, &
                         -end_orb%vec(2) * sin_a + end_orb%vec(4) * cos_a]

    px = end_orb%vec(2)
    px_end = px + kick
    p2_bend = rel_p**2 - end_orb%vec(4)**2
    p_bend = sqrt(p2_bend)

    if (abs(px_end) > p_bend) then
      if (ele%value(hkick$) > abs(ele%value(vkick$))) then
        end_orb%state = lost_pos_x_aperture$
      elseif (ele%value(hkick$) < -abs(ele%value(vkick$))) then
        end_orb%state = lost_neg_x_aperture$
      elseif (ele%value(vkick$) > 0) then
        end_orb%state = lost_pos_y_aperture$
      else
        end_orb%state = lost_neg_y_aperture$
      endif
      return
    endif

    if (abs(kick) < 1d-3 * p_bend) then
      pr = px / p_bend
      kk = kick / p_bend
      sin_g = length * (1 + kk * pr / 2 + (1 + 2*pr**2) * kk**2 / 6 + pr * (3 + 2*pr**2) * kk**3 / 8) / (p_bend * sqrt(1 - pr**2))

      p2xy = p2_bend - px**2
      e = kick * (kick + 2*px) / p2xy
      end_orb%vec(1) = end_orb%vec(1) + length * (kick + 2*px) * (1 + e/4 + e**2/8) / (2 * sqrt(p2xy))
    else
      sin_g = length * (asin(px_end/p_bend) - asin(px/p_bend)) / kick
      end_orb%vec(1) = end_orb%vec(1) + length * (sqrt(p2_bend - px**2) - sqrt(p2_bend - px_end**2)) / kick
    endif

    end_orb%vec(2) = px_end
    end_orb%vec(3) = end_orb%vec(3) + end_orb%vec(4) * sin_g
    end_orb%vec(5) = end_orb%vec(5) + length * end_orb%beta * ele%value(e_tot$) / ele%value(p0c$) - rel_p * sin_g

    end_orb%vec(1:3:2) = [end_orb%vec(1) * cos_a - end_orb%vec(3) * sin_a, &
                          end_orb%vec(1) * sin_a + end_orb%vec(3) * cos_a]
    end_orb%vec(2:4:2) = [end_orb%vec(2) * cos_a - end_orb%vec(4) * sin_a, &
                          end_orb%vec(2) * sin_a + end_orb%vec(4) * cos_a]

    end_orb%t = end_orb%t + rel_p * sin_g / (end_orb%beta * c_light)
  endif

  call offset_particle (ele, param, unset$, end_orb, set_tilt = .false., set_hvkicks = .false.)
  call set_end_orb_s()

!-----------------------------------------------
! drift
 
case (drift$) 

  call offset_particle (ele, param, set$, end_orb)
  call track_a_drift (end_orb, length)
  call offset_particle (ele, param, unset$, end_orb)
  call set_end_orb_s()

!-----------------------------------------------
! elseparator

case (elseparator$)

  call offset_particle (ele, param, set$, end_orb, set_hvkicks = .false.) 

  ! Compute kick
  rtc = abs(rel_tracking_charge) * sign(1, charge_of(end_orb%species))
  hk = ele%value(hkick$) * rtc
  vk = ele%value(vkick$) * rtc

  ! Rotate (x, y) so that kick is in +x direction.

  angle = atan2(vk, hk)
  call tilt_coords (angle, end_orb%vec)

  ! Check if particle can make it though the separator.

  pc = ele%value(p0c$) * (1 + end_orb%vec(6))
  call convert_pc_to (pc, end_orb%species, E_tot = E_tot)
  E_rel = E_tot / ele%value(p0c$)
  mc2 = mass_of(end_orb%species)

  x = end_orb%vec(1)
  px = end_orb%vec(2)
  p_factor = (mc2 / ele%value(p0c$))**2 + end_orb%vec(2)**2 + end_orb%vec(4)**2
  if (length == 0) length = 1d-50  ! To avoid divide by zero
  k_E = sqrt(hk**2 + vk**2) / length

  if (x * k_E < sqrt(p_factor) - E_rel) then
    end_orb = start_orb
    end_orb%state = lost_z_aperture$
    return
  endif

  ! Track

  ps = sqrt((E_rel + k_E * x)**2 - p_factor)
  alpha = length / ps
  coef = k_E * length / ps

  if (abs(coef) > 10) then ! lost
    end_orb%state = lost$
    return
  endif

  if (abs(coef) < 1d-3) then
    sinh_k = alpha * (1 + coef**2 / 6 + coef**4/120)
    cosh1_k = alpha * coef * (1.0_rp / 2 + coef**2 / 24 + coef**4 / 720)
  else
    sinh_k = sinh(coef) / k_E
    cosh1_k = (cosh(coef) - 1) / k_E
  endif

  end_orb%vec(1) = x * cosh(coef) + E_rel * cosh1_k + px * sinh_k
  end_orb%vec(2) = (k_E * x  + E_rel) * sinh(coef) + px * cosh(coef)
  end_orb%vec(3) = end_orb%vec(3) + length * end_orb%vec(4) / ps

  dt = (x * sinh(coef) + E_rel * sinh_k + px * cosh1_k) / c_light
  end_orb%t = end_orb%t + dt
  beta_ref = ele%value(p0c$) / ele%value(e_tot$)
  end_orb%vec(5) = end_orb%vec(5) + end_orb%beta * (length / beta_ref - c_light * dt)

  ! Rotate back to lab coords
  
  call tilt_coords (-angle, end_orb%vec)

  call offset_particle (ele, param, unset$, end_orb, set_hvkicks = .false.) 
  call set_end_orb_s()

!-----------------------------------------------
! kicker
 
case (kicker$, hkicker$, vkicker$) 

  hkick = charge_dir * ele%value(hkick$) 
  vkick = charge_dir * ele%value(vkick$) 
  kick  = charge_dir * ele%value(kick$) 

  call offset_particle (ele, param, set$, end_orb, set_hvkicks = .false.)

  n_slice = max(1, nint(length / ele%value(ds_step$)))
  if (ele%key == hkicker$) then
     end_orb%vec(2) = end_orb%vec(2) + kick / (2 * n_slice)
  elseif (ele%key == vkicker$) then
     end_orb%vec(4) = end_orb%vec(4) + kick / (2 * n_slice)
  else
     end_orb%vec(2) = end_orb%vec(2) + hkick / (2 * n_slice)
     end_orb%vec(4) = end_orb%vec(4) + vkick / (2 * n_slice)
  endif
  do i = 1, n_slice
     call track_a_drift (end_orb, length/n_slice)
     if (i == n_slice) then
        if (ele%key == hkicker$) then
           end_orb%vec(2) = end_orb%vec(2) + kick / (2 * n_slice)
        elseif (ele%key == vkicker$) then
           end_orb%vec(4) = end_orb%vec(4) + kick / (2 * n_slice)
        else
           end_orb%vec(2) = end_orb%vec(2) + hkick / (2 * n_slice)
           end_orb%vec(4) = end_orb%vec(4) + vkick / (2 * n_slice)
        endif
     else 
        if (ele%key == hkicker$) then
           end_orb%vec(2) = end_orb%vec(2) + kick / n_slice
        elseif (ele%key == vkicker$) then
           end_orb%vec(4) = end_orb%vec(4) + kick / n_slice
        else
           end_orb%vec(2) = end_orb%vec(2) + hkick / n_slice
           end_orb%vec(4) = end_orb%vec(4) + vkick / n_slice
        endif
     endif
  end do

  call offset_particle (ele, param, unset$, end_orb, set_hvkicks = .false.)

  if (ele%key == kicker$) then
    end_orb%vec(1) = end_orb%vec(1) + ele%value(h_displace$)
    end_orb%vec(3) = end_orb%vec(3) + ele%value(v_displace$)
  endif

  call set_end_orb_s()

!-----------------------------------------------
! LCavity: Linac rf cavity.
! Modified version of the ultra-relativistic formalism from:
!       J. Rosenzweig and L. Serafini
!       Phys Rev E, Vol. 49, p. 1599, (1994)
! with b_0 = b_-1 = 1. See the Bmad manual for more details.

case (lcavity$)

  if (length == 0) return

  if (ele%value(E_tot_start$) == 0) then
    if (present(err_flag)) err_flag = .true.
    call out_io (s_fatal$, r_name, 'E_TOT_START IS 0 FOR A LCAVITY!' // ele%name)
    if (global_com%exit_on_error) call err_exit
    return
  endif

  call offset_particle (ele, param, set$, end_orb)

  E_start_ref  = ele%value(E_tot_start$)
  E_end_ref    = ele%value(E_tot$)
  gradient_ref = (E_end_ref - E_start_ref) / length
  pc_start_ref = ele%value(p0c_start$)
  pc_end_ref   = ele%value(p0c$)
  beta_start_ref = pc_start_ref / E_start_ref
  beta_end_ref   = pc_end_ref / E_end_ref

  pc_start = pc_start_ref * rel_p
  call convert_pc_to (pc_start, end_orb%species, E_tot = E_start, beta = beta_start)

  ! The RF phase is defined with respect to the time at the beginning of the element.
  ! So if dealing with a slave element and absolute time tracking then need to correct.

  phase = twopi * (ele%value(phi0_err$) + ele%value(phi0_autoscale$) + &
             ele%value(phi0$) + ele%value(phi0_multipass$) + &
             (particle_ref_time (end_orb, ele) - rf_ref_time_offset(ele)) * ele%value(rf_frequency$))

  gradient_max = e_accel_field(ele, gradient$)

  cos_phi = cos(phase)
  gradient_net = gradient_max * cos_phi + gradient_shift_sr_wake(ele, param)

  dE = gradient_net * length
  E_end = E_start + dE
  if (E_end <= mass_of(end_orb%species)) then
    end_orb%state = lost_z_aperture$
    end_orb%vec(6) = -1.01  ! Something less than -1
    return
  endif

  call convert_total_energy_to (E_end, end_orb%species, pc = pc_end, beta = beta_end)
  E_ratio = E_end / E_start
  end_orb%beta = beta_end
  mc2 = mass_of(end_orb%species)

  ! Coupler kick

  call rf_coupler_kick (ele, param, first_track_edge$, phase, end_orb)

  ! Body tracking longitudinal

  end_orb%vec(6) = (pc_end - pc_end_ref) / pc_end_ref 
  end_orb%p0c = pc_end_ref

  if (abs(dE) <  1d-4*(pc_end+pc_start)) then
    dp_dg = length * (E_start / pc_start - mc2**2 * dE / (2 * pc_start**3) + (mc2 * dE)**2 * E_start / (2 * pc_start**5))
  else
    dp_dg = (pc_end - pc_start) / gradient_net
  endif
  end_orb%vec(5) = end_orb%vec(5) * (beta_end / beta_start) - beta_end * (dp_dg - c_light * ele%value(delta_ref_time$))

  ! Body tracking transverse. Kick is only with standing wave cavities.

  if (nint(ele%value(cavity_type$)) == traveling_wave$) then
    end_orb%vec(2) = end_orb%vec(2) / rel_p    ! Convert to x'
    end_orb%vec(4) = end_orb%vec(4) / rel_p    ! Convert to y'

    end_orb%vec(1) = end_orb%vec(1) + end_orb%vec(2) * length
    end_orb%vec(3) = end_orb%vec(3) + end_orb%vec(4) * length

    end_orb%vec(5) = end_orb%vec(5) - (end_orb%vec(2)**2 + end_orb%vec(4)**2) * beta_end * dp_dg / 2

    end_orb%vec(2) = end_orb%vec(2) * (1 + end_orb%vec(6))  ! Convert back to px
    end_orb%vec(4) = end_orb%vec(4) * (1 + end_orb%vec(6))  ! Convert back to py


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

    end_orb%vec(2) = end_orb%vec(2) / rel_p    ! Convert to x'
    end_orb%vec(4) = end_orb%vec(4) / rel_p    ! Convert to y'

    k1 = -gradient_net / (2 * E_start)
    end_orb%vec(2) = end_orb%vec(2) + k1 * end_orb%vec(1)    ! Entrance kick
    end_orb%vec(4) = end_orb%vec(4) + k1 * end_orb%vec(3)    ! Entrance kick

    xp1 = end_orb%vec(2)
    yp1 = end_orb%vec(4)

    end_orb%vec(1:2) = matmul(r_mat, end_orb%vec(1:2))   ! R&S Eq 9.
    end_orb%vec(3:4) = matmul(r_mat, end_orb%vec(3:4))

    xp2 = end_orb%vec(2)
    yp2 = end_orb%vec(4)

    ! Correction of z for finite transverse velocity assumes a uniform change in slope.
    end_orb%vec(5) = end_orb%vec(5) - (xp1**2 + xp2**2 + xp1*xp2 + yp1**2 + yp2**2 + yp1*yp2) * beta_end * dp_dg / 6
    !

    k2 = gradient_net / (2 * E_end) 
    end_orb%vec(2) = end_orb%vec(2) + k2 * end_orb%vec(1)         ! Exit kick
    end_orb%vec(4) = end_orb%vec(4) + k2 * end_orb%vec(3)         ! Exit kick

    end_orb%vec(2) = end_orb%vec(2) * (1 + end_orb%vec(6))  ! Convert back to px
    end_orb%vec(4) = end_orb%vec(4) * (1 + end_orb%vec(6))  ! Convert back to py

  endif

  ! Coupler kick

  call rf_coupler_kick (ele, param, second_track_edge$, phase, end_orb)

  call offset_particle (ele, param, unset$, end_orb)

  ! Time & s calc

  f = gradient_net * length * mc2**2 / (pc_start**2 * E_start)

  if (abs(f) < 1d-6) then
    end_orb%t = start2_orb%t + length * (E_start / pc_start) * (1 - f/2) / c_light
  else
    end_orb%t = start2_orb%t + (pc_end - pc_start) / (gradient_net * c_light)
  endif

  call set_end_orb_s()

!-----------------------------------------------
! marker, etc.
! Note: floor_shift elements can have finite length in the case where it is a slice_slave of a taylor 
! element (the first slice is a taylor element and all other slices are floor_shifts).

case (marker$, fork$, photon_fork$, floor_shift$, fiducial$, detector$)

  end_orb%t = start2_orb%t + ele%value(delta_ref_time$)

  return

!-----------------------------------------------
! mask

case (mask$)

  ! If the plate/mask is turned off then all particles are simply transmitted through.

  if (.not. ele%is_on) return
  
  ! Particle is lost if in an opaque section

  temp_orb = end_orb
  call offset_particle (ele, param, set$, temp_orb, set_multipoles = .false., set_hvkicks = .false.)

  ix_sec = diffraction_plate_or_mask_hit_spot (ele, temp_orb)
  if (ix_sec == 0) end_orb%state = lost$

!-----------------------------------------------
! match

case (match$)

  if (is_true(ele%value(match_end_orbit$))) then
    ele%value(x0$)  = start2_orb%vec(1)
    ele%value(px0$) = start2_orb%vec(2)
    ele%value(y0$)  = start2_orb%vec(3)
    ele%value(py0$) = start2_orb%vec(4)
    ele%value(z0$)  = start2_orb%vec(5)
    ele%value(pz0$) = start2_orb%vec(6)
    end_orb%vec = [ ele%value(x1$), ele%value(px1$), &
                ele%value(y1$), ele%value(py1$), &
                ele%value(z1$), ele%value(pz1$) ]
    return
  endif

  ! Until match_end = False, use unit matrix.

  call match_ele_to_mat6 (ele, vec0, mat6, err)
  if (err) then
    ! Since there are cases where this error may be raised many 
    ! times, do not print an error message.
    if (present(err_flag)) err_flag = .true.
    end_orb%state = lost$
    return
  endif

  end_orb%vec = matmul (mat6, end_orb%vec) + vec0

  call time_and_s_calc ()

!-----------------------------------------------
! multipole, ab_multipole

case (multipole$, ab_multipole$) 

  call offset_particle (ele, param, set$, end_orb, set_multipoles = .false., set_tilt = .false.)

  call multipole_ele_to_kt(ele, .true., has_nonzero_pole, knl, tilt)
  if (has_nonzero_pole) call multipole_kicks (knl*charge_dir, tilt, end_orb, ref_orb_offset = (ele%key == multipole$))

  call offset_particle (ele, param, unset$, end_orb, set_multipoles = .false., set_tilt = .false.)

  call set_end_orb_s()

!-----------------------------------------------
! octupole
! The octupole is modeled using kick-drift.

case (octupole$)

  n_slice = max(1, nint(length / ele%value(ds_step$)))

  k3l = charge_dir * ele%value(k3$) * length / n_slice

  call offset_particle (ele, param, set$, end_orb)

  end_orb%vec(2) = end_orb%vec(2) + k3l *  (3*end_orb%vec(1)*end_orb%vec(3)**2 - end_orb%vec(1)**3) / 12
  end_orb%vec(4) = end_orb%vec(4) + k3l *  (3*end_orb%vec(3)*end_orb%vec(1)**2 - end_orb%vec(3)**3) / 12

  do i = 1, n_slice

    call track_a_drift (end_orb, length / n_slice)

    if (i == n_slice) then
      end_orb%vec(2) = end_orb%vec(2) + k3l *  (3*end_orb%vec(1)*end_orb%vec(3)**2 - end_orb%vec(1)**3) / 12
      end_orb%vec(4) = end_orb%vec(4) + k3l *  (3*end_orb%vec(3)*end_orb%vec(1)**2 - end_orb%vec(3)**3) / 12
    else
      end_orb%vec(2) = end_orb%vec(2) + k3l *  (3*end_orb%vec(1)*end_orb%vec(3)**2 - end_orb%vec(1)**3) / 6
      end_orb%vec(4) = end_orb%vec(4) + k3l *  (3*end_orb%vec(3)*end_orb%vec(1)**2 - end_orb%vec(3)**3) / 6
    endif

  enddo

  call offset_particle (ele, param, unset$, end_orb)

  call set_end_orb_s()

!-----------------------------------------------
! patch

case (patch$)

  call track_a_patch(ele, end_orb)

!-----------------------------------------------
! quadrupole

case (quadrupole$)

  k1 = charge_dir * ele%value(k1$) / rel_p

  call offset_particle (ele, param, set$, end_orb)

  ! Entrance edge

  call hard_multipole_edge_kick (ele, param, first_track_edge$, end_orb)
  call soft_quadrupole_edge_kick (ele, param, first_track_edge$, end_orb)

  ! Body

  call quad_mat2_calc (-k1, length, rel_p, mat2, dz_coef)
  end_orb%vec(5) = end_orb%vec(5) + dz_coef(1) * end_orb%vec(1)**2 + &
                      dz_coef(2) * end_orb%vec(1) * end_orb%vec(2) + dz_coef(3) * end_orb%vec(2)**2 

  end_orb%vec(1:2) = matmul(mat2, end_orb%vec(1:2))

  call quad_mat2_calc (k1, length, rel_p, mat2, dz_coef)
  end_orb%vec(5) = end_orb%vec(5) + dz_coef(1) * end_orb%vec(3)**2 + &
                      dz_coef(2) * end_orb%vec(3) * end_orb%vec(4) + dz_coef(3) * end_orb%vec(4)**2 

  end_orb%vec(3:4) = matmul(mat2, end_orb%vec(3:4))

  ! Exit edge

  call soft_quadrupole_edge_kick (ele, param, second_track_edge$, end_orb)
  call hard_multipole_edge_kick (ele, param, second_track_edge$, end_orb)

  call offset_particle (ele, param, unset$, end_orb)  

  call track1_low_energy_z_correction (end_orb, ele, param)
  call time_and_s_calc ()

!-----------------------------------------------
! rfcavity

case (rfcavity$)

  beta_ref = ele%value(p0c$) / ele%value(e_tot$)
  n_slice = max(1, nint(length / ele%value(ds_step$))) 
  dt_ref = length / (c_light * beta_ref)

  call offset_particle (ele, param, set$, end_orb)

  ! The cavity field is modeled as a standing wave antisymmetric wrt the center.
  ! Thus if the cavity is flipped (orientation = -1), the wave of interest, which is 
  ! always the accelerating wave, is the "backward" wave.

  voltage = e_accel_field(ele, voltage$) * charge_dir

  phase0 = twopi * (ele%value(phi0$) + ele%value(phi0_multipass$) + ele%value(phi0_autoscale$) - &
          (particle_ref_time (end_orb, ele) - rf_ref_time_offset(ele)) * ele%value(rf_frequency$))
  phase = phase0

  t0 = end_orb%t

  call rf_coupler_kick (ele, param, first_track_edge$, phase, end_orb)

  ! Track through slices.
  ! The phase of the accelerating wave traveling in the same direction as the particle is
  ! assumed to be traveling with a phase velocity the same speed as the reference velocity.

  do i = 0, n_slice

    dE = voltage * sin(phase) / n_slice
    if (i == 0 .or. i == n_slice) dE = dE / 2

    call apply_energy_kick (dE, end_orb)
    
    if (end_orb%vec(6) == -1) then
      end_orb%state = lost_z_aperture$
      return
    endif

    if (i /= n_slice) then
      call track_a_drift (end_orb, length/n_slice)
      phase = phase0 + twopi * ele%value(rf_frequency$) * ((i + 1) * dt_ref/n_slice - (end_orb%t - t0)) 
    endif

  enddo

  ! coupler kick

  call rf_coupler_kick (ele, param, second_track_edge$, phase, end_orb)

  call offset_particle (ele, param, unset$, end_orb)

  call set_end_orb_s()

!-----------------------------------------------
! sad_multipole

case (sad_mult$)

  call sad_mult_track_and_mat (ele, param, start_orb, end_orb)
  call set_end_orb_s()

!-----------------------------------------------
! sbend

case (sbend$)

  call track_a_bend (start_orb, ele, param, end_orb)
  call time_and_s_calc ()

!-----------------------------------------------
! sextupole
! The sextupole is modeled using kick-drift.

case (sextupole$)

  n_slice = max(1, nint(length / ele%value(ds_step$)))

  call offset_particle (ele, param, set$, end_orb)

  do i = 0, n_slice
    k2l = charge_dir * ele%value(k2$) * length / n_slice
    if (i == 0 .or. i == n_slice) k2l = k2l / 2
    end_orb%vec(2) = end_orb%vec(2) + k2l * (end_orb%vec(3)**2 - end_orb%vec(1)**2)/2
    end_orb%vec(4) = end_orb%vec(4) + k2l * end_orb%vec(1) * end_orb%vec(3)
    if (i /= n_slice) call track_a_drift (end_orb, length/n_slice)
  enddo

  call offset_particle (ele, param, unset$, end_orb)
  call set_end_orb_s()

!-----------------------------------------------
! Solenoid
! Notice that ks is independent of the ele orientation

case (solenoid$)

  call offset_particle (ele, param, set$, end_orb, set_hvkicks = .false.)
  call solenoid_track_and_mat (ele, param, end_orb, end_orb)
  call offset_particle (ele, param, unset$, end_orb, set_hvkicks = .false.)

  if (ele%value(hkick$) /= 0 .or. ele%value(vkick$) /= 0) then
    ks = rel_tracking_charge * ele%value(ks$)
    ksr = ks / rel_p
    kss = ksr * length
    if (abs(kss) < 1d-2) then
      cos_a = length * (1 - kss**2 / 12 + kss**4 / 360) / 2
      sin_a = (1 - kss**2 / 12 + kss**4 / 240)
      f = length * kss * (1 - kss**2 / 20 + kss**4 / 840) / 6
    else
      cos_a = (1 - cos(kss)) / (kss * ksr)
      sin_a = (kss + sin(kss)) / (2 * kss)
      f = (kss - sin(kss)) / (kss * ksr)
    endif

    end_orb%vec(1:4) = end_orb%vec(1:4) + ele%value(hkick$) * [cos_a, sin_a, -f, -cos_a * ksr / 2] + &
                                          ele%value(vkick$) * [f, cos_a * ksr / 2, cos_a, sin_a]
  endif

  call time_and_s_calc ()

!-----------------------------------------------
! sol_quad

case (sol_quad$)

  call offset_particle (ele, param, set$, end_orb)

  ks = rel_tracking_charge * ele%value(ks$)
  k1 = charge_dir * ele%value(k1$)
  vec0 = 0
  vec0(6) = end_orb%vec(6)
  call sol_quad_mat6_calc (ks, k1, length, vec0, mat6, dz4_coef)
  end_orb%vec(5) = end_orb%vec(5) + sum(end_orb%vec(1:4) * matmul(dz4_coef, end_orb%vec(1:4)))   
  end_orb%vec(1:4) = matmul (mat6(1:4,1:4), end_orb%vec(1:4))

  call offset_particle (ele, param, unset$, end_orb)

  call track1_low_energy_z_correction (end_orb, ele, param)
  call time_and_s_calc ()

!-----------------------------------------------
! Taylor

case (taylor$)

  if (orientation == 1) then
    call track1_taylor (start_orb, ele, param, end_orb)

  else
    call taylor_inverse (ele%taylor, taylor1)
    taylor2 = ele%taylor
    ele%taylor = taylor1
    call track1_taylor (start_orb, ele, param, end_orb)
    ele%taylor = taylor2
    call kill_taylor(taylor1)
  endif

  call time_and_s_calc ()

!-----------------------------------------------
! wiggler:
! Only periodic type wigglers are handled here.
! In the horizontal plane the tracking looks like a drift.
! The tracking in the vertical plane is:
!   1) 1/2 the octupole kick at the entrance face.
!   2) Track as a quadrupole through the body
!   3) 1/2 the octupole kick at the exit face.

case (wiggler$, undulator$)

  if (ele%sub_key == map_type$) then
    if (present(err_flag)) err_flag = .true.
    call out_io (s_fatal$, r_name, &
            'MAP_TYPE WIGGLER: ' // ele%name, &
            'HAS TRACKING_METHOD = BMAD_STANDARD.', &
            'THIS IS NOT A POSSIBLE OPTION FOR THE TRACKING_METHOD.')
    if (global_com%exit_on_error) call err_exit
    return
  endif

  call offset_particle (ele, param, set$, end_orb)

  if (ele%value(l_pole$) == 0) then
    k_z = 1d100    ! Something large
  else
    k_z = pi / ele%value(l_pole$)
  endif
  k1 = -charge_dir * 0.5 * (c_light * ele%value(b_max$) / (ele%value(p0c$) * rel_p))**2

  p_factor = 1 - (end_orb%vec(2) / rel_p**2)**2 - (end_orb%vec(4) / rel_p**2)**2
  if (p_factor < 0) then
    end_orb%state = lost_z_aperture$
    return
  endif

  end_orb%vec(5) = end_orb%vec(5) + 0.5 * (length * (end_orb%beta * ele%value(e_tot$) / ele%value(p0c$) - & 
                   1/sqrt(p_factor)) - 0.5*k1*length / k_z**2 * (1 - rel_p**2))

  ! 1/2 of the octupole octupole kick at the entrance face.

  end_orb%vec(4) = end_orb%vec(4) + k1 * length * rel_p * k_z**2 * end_orb%vec(3)**3 / 3

  ! Quadrupole body

  call quad_mat2_calc (k1, length, rel_p, mat2)
  end_orb%vec(1) = end_orb%vec(1) + length * end_orb%vec(2) / rel_p
  end_orb%vec(3:4) = matmul (mat2, end_orb%vec(3:4))

  ! 1/2 of the octupole octupole kick at the exit face.

  end_orb%vec(4) = end_orb%vec(4) + k1 * length * rel_p * k_z**2 * end_orb%vec(3)**3 / 3
  
  p_factor = 1 - (end_orb%vec(2) / rel_p**2)**2 - (end_orb%vec(4) / rel_p**2)**2
  if (p_factor < 0) then
    end_orb%state = lost_z_aperture$
    return
  endif

  end_orb%vec(5) = end_orb%vec(5) + 0.5 * (length * (end_orb%beta * ele%value(e_tot$) / ele%value(p0c$) & 
                   - 1/sqrt(p_factor)) - 0.5*k1*length / k_z**2 * (1 - rel_p**2))
  
  call offset_particle (ele, param, unset$, end_orb)
   
  call track1_low_energy_z_correction (end_orb, ele, param)

  end_orb%t = start2_orb%t + (ele%value(l$) - 0.5*k1*length / k_z**2 * rel_p**2) / (end_orb%beta * c_light)

  call set_end_orb_s()

!-----------------------------------------------
! unknown

case default

  if (present(err_flag)) err_flag = .true.
  call out_io (s_fatal$, r_name, &
          'BMAD_STANDARD TRACKING_METHOD NOT IMPLMENTED FOR: ' // key_name(ele%key), &
          'FOR ELEMENT: ' // ele%name)
  if (global_com%exit_on_error) call err_exit
  return

end select

!------------------------------------------
contains

subroutine time_and_s_calc ()

end_orb%t = start2_orb%t + ele%value(delta_ref_time$) + (start2_orb%vec(5) - end_orb%vec(5)) / (end_orb%beta * c_light)
call set_end_orb_s()

end subroutine time_and_s_calc

!------------------------------------------
! contains

subroutine set_end_orb_s

if (end_orb%direction == 1) then
  end_orb%s = ele%s
else
  end_orb%s = ele%s - ele%value(l$)
endif

end subroutine set_end_orb_s

!--------------------------------------------------------------
! contains

! Rough calculation for change in longitudinal position using:
!      dz = -L * (<x'^2> + <y'^2>)/ 2 
! where <...> means average.
! The formula below assumes a linear change in velocity between 
! the beginning and the end:

subroutine end_z_calc ()

implicit none

end_orb%vec(5) = start2_orb%vec(5) - (length / rel_p**2) * &
      (start2_orb%vec(2)**2 + end_orb%vec(2)**2 + start2_orb%vec(2) * end_orb%vec(2) + &
       start2_orb%vec(4)**2 + end_orb%vec(4)**2 + start2_orb%vec(4) * end_orb%vec(4)) / 6

end subroutine end_z_calc

end subroutine track1_bmad
