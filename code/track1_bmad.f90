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
! Modules Needed:
!   use bmad
!
! Input:
!   start  -- Coord_struct: Starting position
!   ele    -- Ele_struct: Element
!   param  -- Param_struct:
!
! Output:
!   end   -- Coord_struct: End position
!-

#include "CESR_platform.inc"

subroutine track1_bmad (start, ele, param, end)

  use bmad_struct
  use bmad_interface
  use track1_mod
  use multipole_mod
  use bookkeeper_mod

  implicit none

  type (coord_struct), intent(in)  :: start
  type (coord_struct), intent(out) :: end
  type (ele_struct),   intent(inout)  :: ele
  type (param_struct), intent(inout) :: param

  type (coord_struct)  c0

  real(rp) x_kick, y_kick, k1, k2, k2l, k3l, length, phase
  real(rp) del, e1, e2, del_x_vel, del_y_vel, sig_x, sig_y, kx, ky, coef
  real(rp) knl(0:n_pole_maxx), tilt(0:n_pole_maxx)
  real(rp) ks, sig_x0, sig_y0, beta, mat6(6,6), mat2(2,2), mat4(4,4)
  real(rp) z_slice(100), s_pos, s_pos_old, vec0(6)
  real(rp) ave_x_vel2, ave_y_vel2, rel_E, dE, ff
  real(rp) x_pos, y_pos, cos_phi, gradient, e_start, e_end, e_ratio
  real(rp) alpha, sin_a, cos_a, f, z_ave, r11, r12, r21, r22
  real(rp) x, y, z, px, py, pz, k, dE0, L, E, pxy2
  real(rp) xp_start, xp_end, yp_start, yp_end, dz4_coef(4,4), dz_coef(3)

  integer i, j, n, n_slice, key

  logical init_needed / .true. /

! initially set end = start

  end = start     ! transfer start to end
  length = ele%value(l$)
  rel_E = 1 + start%vec(6)

!-----------------------------------------------
! select

  key = ele%key
  if (key == sol_quad$ .and. ele%value(k1$) == 0) key = solenoid$
  if (.not. ele%is_on) key = drift$  ! if element is off looks like a drift

  select case (key)

!-----------------------------------------------
! marker

  case (marker$)

    return

!-----------------------------------------------
! drift

  case (drift$, rcollimator$, ecollimator$, monitor$, instrument$) 

    call track_a_drift (end%vec, length)

!-----------------------------------------------
! patch

  case (patch$)

    rel_E = 1 + end%vec(6)

    end%vec(2) = end%vec(2) - ele%value(x_pitch$) * rel_E
    end%vec(4) = end%vec(4) - ele%value(y_pitch$) * rel_E
    end%vec(5) = end%vec(5) + ele%value(x_pitch$) * end%vec(1) + &
                              ele%value(y_pitch$) * end%vec(3) 

    if (ele%value(tilt$) /= 0) call tilt_coords (ele%value(tilt$), end%vec, set$)

    end%vec(1) = end%vec(1) - ele%value(x_offset$)
    end%vec(3) = end%vec(3) - ele%value(y_offset$)
    end%vec(5) = end%vec(5) - ele%value(z_offset$)
    end%vec(6) = end%vec(6) - ele%value(dE_offset$)

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
    call attribute_bookkeeper (ele, param)
    call offset_particle (ele, param, end, set$)

    sig_x0 = ele%value(sig_x$)
    sig_y0 = ele%value(sig_y$)
    n_slice = max(1, nint(ele%value(n_slice$)))
    call bbi_slice_calc (n_slice, ele%value(sig_z$), z_slice)
    s_pos = 0    ! end at the ip
    do i = 1, n_slice
      s_pos_old = s_pos
      s_pos = (end%vec(5) + z_slice(i)) / 2
      end%vec(1) = end%vec(1) + end%vec(2) * (s_pos - s_pos_old)
      end%vec(3) = end%vec(3) + end%vec(4) * (s_pos - s_pos_old)
      if (ele%x%beta == 0) then
        sig_x = sig_x0
        sig_y = sig_y0
      else
        beta = ele%x%beta - 2 * ele%x%alpha * s_pos + ele%x%gamma * s_pos**2
        sig_x = sig_x0 * sqrt(beta / ele%x%beta)
        beta = ele%y%beta - 2 * ele%y%alpha * s_pos + ele%y%gamma * s_pos**2
        sig_y = sig_y0 * sqrt(beta / ele%y%beta)
      endif

      call bbi_kick (end%vec(1)/sig_x, end%vec(3)/sig_y, sig_y/sig_x,  &
                                                                  kx, ky)
      coef = ele%value(bbi_const$) / (n_slice * rel_E)
      end%vec(2) = end%vec(2) + kx * coef
      end%vec(4) = end%vec(4) + ky * coef
    enddo
    end%vec(1) = end%vec(1) - end%vec(2) * s_pos
    end%vec(3) = end%vec(3) - end%vec(4) * s_pos

    call offset_particle (ele, param, end, unset$)  
    call end_z_calc

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
    call z_pitch_correction

!-----------------------------------------------
! quadrupole

  case (quadrupole$)

    call offset_particle (ele, param, end, set$)

    k1 = ele%value(k1$) / rel_E

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
    call z_pitch_correction

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
      phase = twopi * ele%value(phi0$) + ff * z
      k  =  ff * ele%value(voltage$) * cos(phase) / param%beam_energy
    endif

    dE0 =  ele%value(voltage$) * sin(phase) / param%beam_energy
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
    call z_pitch_correction

!-----------------------------------------------
! LCavity: Linac rf cavity
! Ultra-relativistic formalism from:
!       J. Rosenzweig and L. Serafini
!       Phys Rev E, Vol. 49, p. 1599, (1994)
! with b_0 = b_-1 = 1
! bp_com%k_loss is an internal variable used with macroparticles.
! It accounts for the short-range wakefields between macroparticles.
! Without macroparticles it should be zero.

  case (lcavity$)

    if (ele%value(energy_start$) == 0) then
      print *, 'ERROR IN TRACK1_BMAD: ENERGY_START IS 0 FOR A LCAVITY!'
      call err_exit
    endif

    phase = twopi * (ele%value(phi0$) - &
                        end%vec(5) * ele%value(rf_frequency$) / c_light)
    cos_phi = cos(phase)
    gradient = ele%value(gradient$) * cos_phi 
    if (bmad_com%sr_wakes_on) gradient = gradient - bmad_com%k_loss - &
                                    ele%value(e_loss$) * param%charge / length

    e_start = ele%value(energy_start$) * rel_E 
    e_end = e_start + gradient * ele%value(l$)
    e_ratio = e_end / e_start

    if (e_ratio < 0) then
      if (bmad_status%type_out) print *, &
                'ERROR IN TRACK1_BMAD: NEGATIVE BEAM ENERGY FOR: ', ele%name
      if (bmad_status%exit_on_error) call err_exit
    endif

    if (gradient == 0) then
      call track_a_drift (end%vec, length)
      return
    endif

    call offset_particle (ele, param, end, set$)

! entrence kick

    k1 = -gradient / (2 * e_start)            
    end%vec(2) = end%vec(2) + k1 * end%vec(1)
    end%vec(4) = end%vec(4) + k1 * end%vec(3)

! track body

    if (bmad_com%use_liar_lcavity) then  ! use liar formula
      r11 = 1
      r12 = e_start * log (e_ratio) / gradient
      r21 = 0
      r22 = 1 / e_ratio
    else
      alpha = log(e_ratio) / (2 * sqrt_2 * cos_phi)
      cos_a = cos(alpha)
      sin_a = sin(alpha)
      f = gradient / (2 * sqrt_2 * cos_phi)
      r11 =  cos_a
      r12 =  sin_a * e_start / f
      r21 = -sin_a * f / e_end
      r22 =  cos_a * e_start / e_end
    endif

    x_pos = end%vec(1)
    y_pos = end%vec(3)

    end%vec(1) = r11 * x_pos + r12 * end%vec(2)
    end%vec(2) = r21 * x_pos + r22 * end%vec(2)
    end%vec(3) = r11 * y_pos + r12 * end%vec(4)
    end%vec(4) = r21 * y_pos + r22 * end%vec(4)

! exit kick

    k2 = gradient / (2 * e_end) 
    end%vec(2) = end%vec(2) + k2 * end%vec(1)
    end%vec(4) = end%vec(4) + k2 * end%vec(3)

! cleanup

    end%vec(6) = (e_end - ele%value(beam_energy$)) / ele%value(beam_energy$) 
    call offset_particle (ele, param, end, unset$)
    call end_z_calc

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
    call z_pitch_correction

!-----------------------------------------------
! solenoid

  case (solenoid$)

    call offset_particle (ele, param, end, set$)

    ks = ele%value(ks$) / rel_E

    xp_start = end%vec(2) + ks * end%vec(3) / 2
    yp_start = end%vec(4) - ks * end%vec(1) / 2
    end%vec(5) = end%vec(5) - length * (xp_start**2 + yp_start**2 ) / 2

    call solenoid_mat_calc (ks, length, mat4)
    end%vec(1:4) = matmul (mat4, end%vec(1:4))

    call offset_particle (ele, param, end, unset$)
    call z_pitch_correction

!-----------------------------------------------
! sol_quad

  case (sol_quad$)

    call offset_particle (ele, param, end, set$)

    ks = ele%value(ks$) / rel_E
    k1 = ele%value(k1$) / rel_E
    vec0 = 0
    call sol_quad_mat6_calc (ks, k1, length, mat6, vec0, dz4_coef)
    end%vec(5) = end%vec(5) + sum(end%vec(1:4) * matmul(dz4_coef, end%vec(1:4)))   
    end%vec(1:4) = matmul (mat6(1:4,1:4), end%vec(1:4))

    call offset_particle (ele, param, end, unset$)
    call z_pitch_correction

!-----------------------------------------------
! wiggler:

  case (wiggler$)

    if (ele%sub_key == map_type$) then
      print *, 'ERROR IN TRACK1_BMAD: NEW STYLE WIGGLER: ', ele%name
      print *, '       HAS TRACKING_METHOD = BMAD_STANDARD.'
      print *, '       THIS IS NOT A POSSIBLE OPTION FOR THE TRACKING_METHOD.'
      call err_exit
    endif

    call offset_particle (ele, param, end, set$, set_multipoles=.false.)

    k1 = ele%value(k1$) / rel_E**2
    call quad_mat2_calc (k1, length, mat2)
    end%vec(1) = end%vec(1) + length * end%vec(2)
    end%vec(3:4) = matmul (mat2, end%vec(3:4))

    call offset_particle (ele, param, end, unset$, set_multipoles=.false.)
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

    call track_taylor (start%vec, ele%taylor, end%vec)

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

  end%vec(5) = end%vec(5) - (length / rel_E**2) * &
        (start%vec(2)**2 + end%vec(2)**2 + start%vec(2) * end%vec(2) + &
         start%vec(4)**2 + end%vec(4)**2 + start%vec(4) * end%vec(4)) / 6

end subroutine

!-------------------------------------------------------------------

subroutine z_pitch_correction

  real(rp) x_pitch, y_pitch

  x_pitch = ele%value(x_pitch_tot$)
  y_pitch = ele%value(y_pitch_tot$)

  if (x_pitch /= 0) end%vec(5) = end%vec(5) + length*x_pitch**2/2 - &
                   x_pitch * (end%vec(1) - start%vec(1))

  if (y_pitch /= 0) end%vec(5) = end%vec(5) + length*y_pitch**2/2 - &
                   y_pitch * (end%vec(3) - start%vec(3))

end subroutine

end subroutine
