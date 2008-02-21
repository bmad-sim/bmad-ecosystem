!+
! Subroutine track1_bmad (start, ele, param, end)
!
! Particle tracking through a single element BMAD_standard style.
! This routine is NOT meant for long term tracking since it does not get 
! all the 2nd order terms for the longitudinal motion.
!
! It is assumed that HKICK and VKICK are the kicks in the horizontal
! and vertical kicks irregardless of the value for TILT.
!
! Note: track1_bmad *never* relies on ele%mat6 for tracking excect for 
! hybrid elements.
! 
! Note: end%vec(6) will be set < -1 (and param%lost will be set) if the 
! particle fails to make it through an lcavity
!
! Modules Needed:
!   use bmad
!
! Input:
!   start  -- Coord_struct: Starting position
!   ele    -- Ele_struct: Element
!   param  -- lat_param_struct:
!     %particle -- Particle type
!
! Output:
!   end   -- Coord_struct: End position
!   param  -- lat_param_struct:
!     %lost -- Set True if particle is lost. False otherwise.
!-

#include "CESR_platform.inc"

subroutine track1_bmad (start, ele, param, end)

  use bmad, except_dummy => track1_bmad

  implicit none

  type (coord_struct) :: start
  type (coord_struct) :: end
  type (ele_struct) :: ele
  type (lat_param_struct) :: param

  real(rp) k1, k2, k2l, k3l, length, phase, beta_start
  real(rp) beta_end, beta_start_ref, beta_end_ref
  real(rp) e2, sig_x, sig_y, kx, ky, coef, bbi_const
  real(rp) knl(0:n_pole_maxx), tilt(0:n_pole_maxx)
  real(rp) ks, sig_x0, sig_y0, beta, mat6(6,6), mat2(2,2), mat4(4,4)
  real(rp) z_slice(100), s_pos, s_pos_old, vec0(6)
  real(rp) rel_pc, ff, k_z, pc_start, pc_end, dt_ref, gradient_ref
  real(rp) x_pos, y_pos, cos_phi, gradient, e_start, e_end, e_ratio
  real(rp) alpha, sin_a, cos_a, f, r11, r12, r21, r22
  real(rp) x, y, z, px, py, pz, k, dE0, L, E, pxy2, xp0, xp1, yp0, yp1
  real(rp) xp_start, yp_start, dz4_coef(4,4), dz_coef(3)
  real(rp) dp_coupler, dp_x_coupler, dp_y_coupler
  real(rp) phase0, dphase, dcos_phi, dgradient, dpz, d_ct
  real(rp) E_start_ref, E_end_ref, pc_start_ref, pc_end_ref
  real(rp), parameter :: phase_cut = 1e-5

  integer i, n, n_slice, key

! initially set end = start

  end = start     ! transfer start to end
  length = ele%value(l$)
  rel_pc = 1 + start%vec(6)

!-----------------------------------------------
! Select
! If element is off looks like a drift. LCavities will still do wakefields.

  key = ele%key
  if (key == sol_quad$ .and. ele%value(k1$) == 0) key = solenoid$
  if (.not. ele%is_on .and. key /= lcavity$) key = drift$  

  select case (key)

!-----------------------------------------------
! marker

  case (marker$)

    return

!-----------------------------------------------
! drift

  case (drift$, rcollimator$, ecollimator$, monitor$, instrument$) 

    if (ele%is_on) call offset_particle (ele, param, end, set$, .false.)
    call track_a_drift (end%vec, length)
    if (ele%is_on) call offset_particle (ele, param, end, unset$, .false.)

!-----------------------------------------------
! patch

  case (patch$)

    rel_pc = 1 + end%vec(6)

    end%vec(2) = end%vec(2) - ele%value(x_pitch$) * rel_pc
    end%vec(4) = end%vec(4) - ele%value(y_pitch$) * rel_pc
    end%vec(5) = end%vec(5) + ele%value(x_pitch$) * end%vec(1) + &
                              ele%value(y_pitch$) * end%vec(3) 

    if (ele%value(tilt$) /= 0) call tilt_coords (ele%value(tilt$), end%vec, set$)

    end%vec(1) = end%vec(1) - ele%value(x_offset$)
    end%vec(3) = end%vec(3) - ele%value(y_offset$)
    end%vec(5) = end%vec(5) - ele%value(z_offset$)
    end%vec(6) = end%vec(6) - ele%value(pz_offset$)

!-----------------------------------------------
! match

  case (match$)

    call match_ele_to_mat6 (ele, vec0, mat6)
    end%vec = matmul (mat6, end%vec) + vec0

!-----------------------------------------------
! kicker, separator

  case (elseparator$, kicker$, hkicker$, vkicker$) 

    call offset_particle (ele, param, end, set$)

    end%vec(1) = end%vec(1) + length * end%vec(2)
    end%vec(3) = end%vec(3) + length * end%vec(4)

    if (ele%key == kicker$) then
      end%vec(1) = end%vec(1) + ele%value(h_displace$)
      end%vec(3) = end%vec(3) + ele%value(v_displace$)
    endif

    call offset_particle (ele, param, end, unset$)  
    call end_z_calc

!-----------------------------------------------
! beambeam
                          
  case (beambeam$)

    if (ele%value(charge$) == 0 .or. param%n_part == 0) return

    sig_x0 = ele%value(sig_x$)
    sig_y0 = ele%value(sig_y$)
    if (sig_x0 == 0 .or. sig_y0 == 0) return

    call offset_particle (ele, param, end, set$)

    n_slice = max(1, nint(ele%value(n_slice$)))
    call bbi_slice_calc (n_slice, ele%value(sig_z$), z_slice)
    s_pos = 0    ! end at the ip
    do i = 1, n_slice
      s_pos_old = s_pos
      s_pos = (end%vec(5) + z_slice(i)) / 2
      end%vec(1) = end%vec(1) + end%vec(2) * (s_pos - s_pos_old)
      end%vec(3) = end%vec(3) + end%vec(4) * (s_pos - s_pos_old)
      if (ele%a%beta == 0) then
        sig_x = sig_x0
        sig_y = sig_y0
      else
        beta = ele%a%beta - 2 * ele%a%alpha * s_pos + ele%a%gamma * s_pos**2
        sig_x = sig_x0 * sqrt(beta / ele%a%beta)
        beta = ele%b%beta - 2 * ele%b%alpha * s_pos + ele%b%gamma * s_pos**2
        sig_y = sig_y0 * sqrt(beta / ele%b%beta)
      endif

      call bbi_kick (end%vec(1)/sig_x, end%vec(3)/sig_y, sig_y/sig_x,  &
                                                                  kx, ky)
      bbi_const = -param%n_part * m_electron * ele%value(charge$) * r_e /  &
                      (2 * pi * ele%value(p0c$) * (sig_x + sig_y))
      coef = ele%value(bbi_const$) / (n_slice * rel_pc)
      end%vec(2) = end%vec(2) + kx * coef
      end%vec(4) = end%vec(4) + ky * coef
    enddo
    end%vec(1) = end%vec(1) - end%vec(2) * s_pos
    end%vec(3) = end%vec(3) - end%vec(4) * s_pos

    call offset_particle (ele, param, end, unset$)  

!-----------------------------------------------
! octupole
! The octupole is treated as a thin lens with a position dependent kick
! at the beginning and the end

  case (octupole$)

    k3l = ele%value(k3$) * length

    call offset_particle (ele, param, end, set$, set_canonical = .false.)

    end%vec(2) = end%vec(2) + k3l *  &
                    (3*end%vec(1)*end%vec(3)**2 - end%vec(1)**3) / 12
    end%vec(4) = end%vec(4) + k3l *  &
                    (3*end%vec(3)*end%vec(1)**2 - end%vec(3)**3) / 12

    call track_a_drift (end%vec, length)

    end%vec(2) = end%vec(2) + k3l *  &
                    (3*end%vec(1)*end%vec(3)**2 - end%vec(1)**3) / 12
    end%vec(4) = end%vec(4) + k3l *  &
                    (3*end%vec(3)*end%vec(1)**2 - end%vec(3)**3) / 12

    call offset_particle (ele, param, end, unset$, set_canonical = .false.)  

!-----------------------------------------------
! quadrupole

  case (quadrupole$)

    call offset_particle (ele, param, end, set$)

    k1 = ele%value(k1$) / rel_pc

    call quad_mat2_calc (-k1, length, mat2, dz_coef)
    end%vec(5) = end%vec(5) + dz_coef(1) * end%vec(1)**2 + &
                              dz_coef(2) * end%vec(1) * end%vec(2) + &
                              dz_coef(3) * end%vec(2)**2 

    end%vec(1:2) = matmul(mat2, end%vec(1:2))

    call quad_mat2_calc (k1, length, mat2, dz_coef)
    end%vec(5) = end%vec(5) + dz_coef(1) * end%vec(3)**2 + &
                              dz_coef(2) * end%vec(3) * end%vec(4) + &
                              dz_coef(3) * end%vec(4)**2 

    end%vec(3:4) = matmul(mat2, end%vec(3:4))

    call offset_particle (ele, param, end, unset$)  

!-----------------------------------------------
! sbend

  case (sbend$)

    call track_a_bend (start, ele, param, end)

!-----------------------------------------------
! rfcavity

  case (rfcavity$)

    call offset_particle (ele, param, end, set$, set_canonical = .false.)

    x = end%vec(1)
    y = end%vec(3)
    z = end%vec(5)

    px = end%vec(2)
    py = end%vec(4)
    pz = end%vec(6)

    if (ele%value(voltage$) == 0) then
      phase = 0
      k = 0
    else
      if (ele%value(RF_frequency$) == 0) then
        print *, 'ERROR IN TRACK1_BMAD: ', &
                   '"RF_FREQUENCY" ATTRIBUTE NOT SET FOR RF: ', trim(ele%name)
        print *, '      YOU NEED TO SET THIS OR THE "HARMON" ATTRIBUTE.'
        call err_exit
      endif
      ff = twopi * ele%value(rf_frequency$) / c_light
      phase = twopi * (ele%value(phi0$)+ele%value(dphi0$)) + ff * z
      k  =  ff * ele%value(voltage$) * cos(phase) / ele%value(p0c$)
    endif

    dE0 =  ele%value(voltage$) * sin(phase) / ele%value(E_tot$)
    L = ele%value(l$)
    E = 1 + pz
    E2 = E**2
    pxy2 = px**2 + py**2

!

    end = start
    end%vec(1) = x + px*L * (1/E - dE0/2 + pxy2*L/12 + pz*dE0 + dE0**2/3) 
    end%vec(3) = y + py*L * (1/E - dE0/2 + pxy2*L/12 + pz*dE0 + dE0**2/3)
    end%vec(5) = z + pxy2*L * (-1/(2*E2) + dE0/2)
    end%vec(6) = pz + dE0 + k*pxy2*L * (-1/(4*E2) + dE0/6) 
         
    call offset_particle (ele, param, end, unset$, set_canonical = .false.)

!-----------------------------------------------
! LCavity: Linac rf cavity
! Ultra-relativistic formalism from:
!       J. Rosenzweig and L. Serafini
!       Phys Rev E, Vol. 49, p. 1599, (1994)
! with b_0 = b_-1 = 1. The extra factors of beta are included to make the 
! transverse determinant (beta_i*gamma_i)/(beta_f*gamma_f) which it should
! be at low energies.

  case (lcavity$)

    if (length == 0) return

    if (ele%value(E_tot_start$) == 0) then
      print *, 'ERROR IN TRACK1_BMAD: E_TOT_START IS 0 FOR A LCAVITY!'
      call err_exit
    endif

    dphase = twopi * (ele%value(phi0_err$) - &
                              end%vec(5) * ele%value(rf_frequency$) / c_light)
    phase0 = twopi * (ele%value(phi0$) + ele%value(dphi0$)) 
    phase = phase0 + dphase

    cos_phi = cos(phase)
    gradient = (ele%value(gradient$) + ele%value(gradient_err$)) * cos_phi 

    if (.not. ele%is_on) then
      gradient = 0
    endif

    ! bp_com%grad_loss_sr_wake is an internal variable used with macroparticles.
    ! It accounts for the longitudinal short-range wakefields between macroparticles.
    ! Without macroparticles it should be zero.
    ! Since the reference energy is shifted by e_loss, dgradient is not affected here.

    if (bmad_com%sr_wakes_on) then
      if (bmad_com%grad_loss_sr_wake /= 0) then  
        ! use grad_loss_sr_wake and ignore e_loss
        gradient = gradient - bmad_com%grad_loss_sr_wake
      else
        gradient = gradient - ele%value(e_loss$) * param%n_part * &
                                                            e_charge / length
      endif
    endif

    E_start_ref  = ele%value(E_tot_start$)
    E_end_ref    = ele%value(E_tot$)
    gradient_ref = (E_end_ref - E_start_ref) / length
    pc_start_ref = ele%value(p0c_start$)
    pc_end_ref   = ele%value(p0c$)
    beta_start_ref = pc_start_ref / E_start_ref
    beta_end_ref   = pc_end_ref / E_end_ref


    ! If the cavity is off and the reference energy does not change then
    ! the tracking is simple.

    if (gradient == 0 .and. gradient_ref == 0) then
      if (ele%is_on) call offset_particle (ele, param, end, set$, .false.)
      call track_a_drift (end%vec, length)
      if (ele%is_on) call offset_particle (ele, param, end, unset$, .false.)
      return
    endif

    pc_start = pc_start_ref * rel_pc
    call convert_pc_to (pc_start, param%particle, &
                                      E_tot = E_start, beta = beta_start)
    E_end = E_start + gradient * length
    if (E_end <= mass_of(param%particle)) then
      param%lost = .true.
      end%vec(6) = -1.01  ! Something less than -1
      return
    endif

    call convert_total_energy_to (E_end, param%particle, &
                                             pc = pc_end, beta = beta_end)
    E_ratio = E_end / E_start

    call offset_particle (ele, param, end, set$)

    ! entrence kick

    k1 = -gradient / (2 * pc_start)
    end%vec(2) = end%vec(2) + k1 * end%vec(1)
    end%vec(4) = end%vec(4) + k1 * end%vec(3)

    ! coupler kick

    if (ele%value(coupler_strength$) /= 0) then

      dp_coupler = (ele%value(gradient$) * ele%value(gradient_err$)) * &
            ele%value(coupler_strength$) * cos(phase + twopi * ele%value(coupler_phase$))
      dp_x_coupler = dp_coupler * cos (twopi * ele%value(coupler_angle$))
      dp_y_coupler = dp_coupler * sin (twopi * ele%value(coupler_angle$))

      if (ele%coupler_at == both_ends$) then
        dp_x_coupler = dp_x_coupler / 2
        dp_y_coupler = dp_y_coupler / 2
      endif

      if (ele%coupler_at == entrance_end$ .or. ele%coupler_at == both_ends$) then
        end%vec(2) = end%vec(2) + dp_x_coupler / pc_start
        end%vec(4) = end%vec(4) + dp_y_coupler / pc_start
      endif

    endif

    ! track body

    if (gradient == 0) then
      r11 = 1
      r12 = length
      r21 = 0
      r22 = 1

    elseif (bmad_com%use_liar_lcavity) then  ! use liar formula
      r11 = 1
      r12 = E_start * log (E_ratio) / gradient
      r21 = 0
      r22 = 1 / E_ratio
    else
      alpha = log(E_ratio) / (2 * sqrt_2 * cos_phi)
      cos_a = cos(alpha)
      sin_a = sin(alpha)
      f = gradient / (2 * sqrt_2 * cos_phi)
      r11 =  cos_a
      r12 =  sin_a * beta_start * E_start / f
      r21 = -sin_a * f / (E_end * beta_end)
      r22 =  cos_a * beta_start * E_start / (E_end * beta_end)
    endif

    x_pos = end%vec(1)
    y_pos = end%vec(3)

    end%vec(1) = r11 * x_pos + r12 * end%vec(2)
    end%vec(2) = r21 * x_pos + r22 * end%vec(2)
    end%vec(3) = r11 * y_pos + r12 * end%vec(4)
    end%vec(4) = r21 * y_pos + r22 * end%vec(4)

    ! coupler kick

    if (ele%value(coupler_strength$) /= 0) then
      if (ele%coupler_at == exit_end$ .or. ele%coupler_at == both_ends$) then
        end%vec(2) = end%vec(2) + dp_x_coupler / pc_end
        end%vec(4) = end%vec(4) + dp_y_coupler / pc_end
      endif
    endif

    ! exit kick

    k2 = gradient / (2 * pc_end) 
    end%vec(2) = end%vec(2) + k2 * end%vec(1)
    end%vec(4) = end%vec(4) + k2 * end%vec(3)

    ! Final momentum.
    ! dgradient is the deviation of the gradient from the reference and is used to 
    ! avoid round-off errors when the change in pz is small.

    if (abs(dphase) < phase_cut .and. ele%is_on) then
      dcos_phi = -sin(phase0) * dphase - cos(phase0) * dphase**2 / 2
      dgradient = ele%value(gradient$) * dcos_phi + ele%value(gradient_err$) * cos(phase0) 
      dpz = (dgradient * length - start%vec(6) * pc_start_ref * &
                           (beta_end_ref - beta_start_ref)) / beta_end_ref 
      end%vec(6) = (dpz + start%vec(6) * pc_start_ref) / pc_end_ref
    else
      end%vec(6) = (pc_end - pc_end_ref) / pc_end_ref 
    endif
    call offset_particle (ele, param, end, unset$)

    ! correct z for change in velocity

    if (gradient_ref == 0) then
      d_ct = (pc_end - pc_start) / gradient - length / beta_end_ref
    elseif (gradient == 0) then
      d_ct = length / beta_end - (pc_end_ref - pc_start_ref) / gradient_ref
    elseif (abs(dphase) < phase_cut) then
      d_ct = dpz / gradient - dgradient * &
                      (pc_end_ref - pc_start_ref) / (gradient * gradient_ref)
    else
      d_ct = (pc_end - pc_start) / gradient - &
                    (pc_end_ref - pc_start_ref) / gradient_ref
    endif                      

    end%vec(5) = end%vec(5) * (beta_end / beta_start) - beta_end * d_ct

    ! This assumes a uniform change in slope.

    xp0 = start%vec(2) / rel_pc
    xp1 = end%vec(2) / (1 + end%vec(6))
    yp0 = start%vec(4) / rel_pc
    yp1 = end%vec(4) / (1 + end%vec(6))

    end%vec(5) = end%vec(5) - (length / 6) * &
            (xp0**2 + xp1**2 + xp0*xp1 + yp0**2 + yp1**2 + yp0*yp1)

!-----------------------------------------------
! sextupole
! The sextupole is treated as a drift with position dependent kick
! at the beginning and the end

  case (sextupole$)

    k2l = ele%value(k2$) * length 

    call offset_particle (ele, param, end, set$, set_canonical = .false.)

    end%vec(2) = end%vec(2) + k2l * (end%vec(3)**2 - end%vec(1)**2)/4
    end%vec(4) = end%vec(4) + k2l * end%vec(1) * end%vec(3) / 2
    call track_a_drift (end%vec, length)
    end%vec(2) = end%vec(2) + k2l * (end%vec(3)**2 - end%vec(1)**2)/4
    end%vec(4) = end%vec(4) + k2l * end%vec(1) * end%vec(3) / 2

    call offset_particle (ele, param, end, unset$, set_canonical = .false.)

!-----------------------------------------------
! solenoid

  case (solenoid$)

    call offset_particle (ele, param, end, set$)

    ks = ele%value(ks$) / rel_pc

    xp_start = end%vec(2) + ks * end%vec(3) / 2
    yp_start = end%vec(4) - ks * end%vec(1) / 2
    end%vec(5) = end%vec(5) - length * (xp_start**2 + yp_start**2 ) / 2

    call solenoid_mat_calc (ks, length, mat4)
    end%vec(1:4) = matmul (mat4, end%vec(1:4))

    call offset_particle (ele, param, end, unset$)

!-----------------------------------------------
! sol_quad

  case (sol_quad$)

    call offset_particle (ele, param, end, set$)

    ks = ele%value(ks$) / rel_pc
    k1 = ele%value(k1$) / rel_pc
    vec0 = 0
    call sol_quad_mat6_calc (ks, k1, length, mat6, vec0, dz4_coef)
    end%vec(5) = end%vec(5) + sum(end%vec(1:4) * matmul(dz4_coef, end%vec(1:4)))   
    end%vec(1:4) = matmul (mat6(1:4,1:4), end%vec(1:4))

    call offset_particle (ele, param, end, unset$)

!-----------------------------------------------
! wiggler:
! Only periodic type wigglers are handled here.
! In the horizontal plane the tracking looks like a drift.
! The tracking in the vertical plane is:
!   1) 1/2 the octupole kick at the entrance face.
!   2) Track as a quadrupole through the body
!   3) 1/2 the octupole kick at the exit face.

  case (wiggler$)

    if (ele%sub_key == map_type$) then
      print *, 'ERROR IN TRACK1_BMAD: NEW STYLE WIGGLER: ', ele%name
      print *, '       HAS TRACKING_METHOD = BMAD_STANDARD.'
      print *, '       THIS IS NOT A POSSIBLE OPTION FOR THE TRACKING_METHOD.'
      call err_exit
    endif

    call offset_particle (ele, param, end, set$)

    if (ele%value(l_pole$) == 0) then
      k_z = 0
    else
      k_z = pi / ele%value(l_pole$)
    endif
    k1 = -0.5 * (c_light * ele%value(b_max$) / &
                    (ele%value(p0c$) * rel_pc))**2

    ! 1/2 of the octupole octopole kick at the entrance face.

    end%vec(4) = end%vec(4) + k1 * length * k_z**2 * end%vec(3)**3 / 3

    ! Quadrupole body

    call quad_mat2_calc (k1, length, mat2)
    end%vec(1) = end%vec(1) + length * end%vec(2)
    end%vec(3:4) = matmul (mat2, end%vec(3:4))

    ! 1/2 of the octupole octopole kick at the exit face.

    end%vec(4) = end%vec(4) + k1 * length * k_z**2 * end%vec(3)**3 / 3

    call offset_particle (ele, param, end, unset$)
    call end_z_calc

!-----------------------------------------------
! multipole

  case (multipole$, ab_multipole$) 

    call offset_particle (ele, param, end, set$, &
                          set_canonical = .false., set_multipoles = .false.)

    call multipole_ele_to_kt(ele, param%particle, knl, tilt, .false.)
    do n = 0, n_pole_maxx
      call multipole_kick (knl(n), tilt(n), n, end, .true.)
    enddo

    call offset_particle (ele, param, end, unset$, &
                           set_canonical = .false., set_multipoles = .false.)

!-----------------------------------------------
! accel_sol
! look at former routine in track1_mod:track_a_accel_sol ()

  case (accel_sol$)

    print *, 'ERROR: ACCEL_SOL NOT IMPLEMENTED!' 
    call err_exit

!-----------------------------------------------
! Taylor

  case (taylor$)

    call track1_taylor (start, ele, param, end)

!-----------------------------------------------
! unknown

  case default

    print *, 'ERROR IN TRACK1_BMAD: UNKNOWN ELEMENT: ', &
                                        key_name(ele%key), ele%type
    call err_exit

  end select

contains

!--------------------------------------------------------------
! Rough calculation for change in longitudinal position using:
!      dz = -L * (<x'^2> + <y'^2>)/ 2 
! where <...> means average.
! The formula below assumes a linear change in velocity between 
! the beginning and the end:

subroutine end_z_calc

  end%vec(5) = end%vec(5) - (length / rel_pc**2) * &
        (start%vec(2)**2 + end%vec(2)**2 + start%vec(2) * end%vec(2) + &
         start%vec(4)**2 + end%vec(4)**2 + start%vec(4) * end%vec(4)) / 6

end subroutine

end subroutine
