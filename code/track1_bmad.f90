!+
! Subroutine track1_bmad (start, ele, param, end)
!
! Particle tracking through a single element BMAD_standard style.
! This routine is NOT ment for long term tracking since it does not get 
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
  type (ele_struct)  bend

  real(rp) x_kick, y_kick, k1, k2l, k3l, length, phase, mat2(2,2), mat4(4,4)
  real(rp) del, e1, e2, del_x_vel, del_y_vel, sig_x, sig_y, kx, ky, coef
  real(rp) knl(0:n_pole_maxx), tilt(0:n_pole_maxx)
  real(rp) ks, sig_x0, sig_y0, beta, mat6(6,6)
  real(rp) z_slice(100), s_pos, s_pos_old, vec0(6)
  real(rp) ave_x_vel2, ave_y_vel2, dE_E
  real(rp) x_pos, y_pos, cos_phi, gradient, e_start, e_end, e_ratio
  real(rp) alpha, sin_a, cos_a, f, z_ave
  real(rp) x, y, z, px, py, pz, k, dE0, L, E, pxy2

  integer i, j, n, n_slice, key

  logical init_needed / .true. /

! init bend element                    

  if (init_needed) then
    call init_ele (bend)
    bend%key = sbend$
    init_needed = .false.
  endif

! initially set end = start

  end = start     ! transfer start to end
  length = ele%value(l$)

!-----------------------------------------------
! select

  key = ele%key
  if (.not. ele%is_on) key = drift$  ! if element is off looks like a drift

  select case (key)

! marker

  case (marker$)

    return

! drift

  case (drift$, rcollimator$, ecollimator$, monitor$, instrument$) 

    call track_a_drift (end%vec, length)

! patch

  case (patch$)

    dE_E = 1 + end%vec(6)
    end%vec(1) = end%vec(1) - ele%value(x_offset$)
    end%vec(2) = end%vec(2) - ele%value(x_pitch$) * dE_E
    end%vec(3) = end%vec(3) - ele%value(y_offset$)
    end%vec(4) = end%vec(4) - ele%value(y_pitch$) * dE_E
    end%vec(5) = end%vec(5) - ele%value(z_offset$) + &
                              ele%value(x_pitch$) * end%vec(1) + &
                              ele%value(y_pitch$) * end%vec(3) 
    end%vec(6) = end%vec(6) - ele%value(dE_offset$)


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
      coef = ele%value(bbi_const$) / (n_slice * (1 + end%vec(6)))
      end%vec(2) = end%vec(2) + kx * coef
      end%vec(4) = end%vec(4) + ky * coef
    enddo
    end%vec(1) = end%vec(1) - end%vec(2) * s_pos
    end%vec(3) = end%vec(3) - end%vec(4) * s_pos

    call offset_particle (ele, param, end, unset$)  
    call end_z_calc

! octupole
! The octupole is treated as a thin lens with a position dependent kick
! at the beginning and the end

  case (octupole$)

    call offset_particle (ele, param, end, set$)

    k3l = ele%value(k3$) * length / (1 + end%vec(6))

    end%vec(2) = end%vec(2) + k3l *  &
                    (3*end%vec(1)*end%vec(3)**2 - end%vec(1)**3) / 12
    end%vec(4) = end%vec(4) + k3l *  &
                    (3*end%vec(3)*end%vec(1)**2 - end%vec(3)**3) / 12

    end%vec(1) = end%vec(1) + end%vec(2) * length
    end%vec(3) = end%vec(3) + end%vec(4) * length

    end%vec(2) = end%vec(2) + k3l *  &
                    (3*end%vec(1)*end%vec(3)**2 - end%vec(1)**3) / 12
    end%vec(4) = end%vec(4) + k3l *  &
                    (3*end%vec(3)*end%vec(1)**2 - end%vec(3)**3) / 12

    call offset_particle (ele, param, end, unset$)  
    call end_z_calc

! quadrupole

  case (quadrupole$)

    call offset_particle (ele, param, end, set$)

    k1 = ele%value(k1$) / (1 + end%vec(6))
    call quad_mat_calc (-k1, length, mat2)
    end%vec(1:2) = matmul(mat2, end%vec(1:2))
    call quad_mat_calc (k1, length, mat2)
    end%vec(3:4) = matmul(mat2, end%vec(3:4))

    call offset_particle (ele, param, end, unset$)  
    call end_z_calc

! sbend
! A non-zero roll has a zeroth order effect that must be included

  case (sbend$)

    call track_a_bend (end, ele, param, end)

! rfcavity

  case (rfcavity$)

    call offset_particle (ele, param, end, set$, set_canonical = .false.)

    x = end%vec(1)
    y = end%vec(3)
    z = end%vec(5)

    px = end%vec(2)
    py = end%vec(4)
    pz = end%vec(6)

    if (ele%value(volt$) == 0) then
      phase = 0
      k = 0
    else
      if (ele%value(RF_wavelength$) == 0) then
        print *, 'ERROR IN TRACK1_BMAD: ', &
                   '"RF_WAVELENGTH" ATTRIBUTE NOT SET FOR RF: ', trim(ele%name)
        print *, '      YOU NEED TO SET THIS OR THE "HARMON" ATTRIBUTE.'
        call err_exit
      endif
      phase = twopi * (ele%value(phi0$) + z / ele%value(rf_wavelength$))
      k  =  twopi * ele%value(volt$) * cos(phase) / &
                              (param%beam_energy * ele%value(rf_wavelength$))
    endif

    dE0 =  ele%value(volt$) * sin(phase) / param%beam_energy
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

! linac rf cavity
! assumes the particle is ultra-relativistic.
! This uses the formalism from:
!       J. Rosenzweig and L. Serafini
!       Phys Rev E, Vol. 49, p. 1599, (1994)
! with b_0 = b_-1 = 1

  case (lcavity$)

    if (ele%value(energy_start$) == 0) then
      print *, 'ERROR IN TRACK1_BMAD: ENERGY_START IS 0 FOR A LCAVITY!'
      call err_exit
    endif

    phase = twopi * (ele%value(phi0$) + &
                        end%vec(5) * ele%value(rf_frequency$) / c_light)
    cos_phi = cos(phase)
    gradient = ele%value(gradient$) * cos_phi
    e_start = ele%value(energy_start$) * (1 + end%vec(6))
    e_end = e_start + gradient * ele%value(l$)
    e_ratio = e_end / e_start
    if (e_ratio < 0) then
      if (bmad_status%type_out) print *, &
                'ERROR IN TRACK1_BMAD: NEGATIVE BEAM ENERGY FOR: ', ele%name
      if (bmad_status%exit_on_error) call err_exit
    endif
    alpha = log(e_ratio) / (2 * sqrt_2 * cos_phi)
    cos_a = cos(alpha)
    sin_a = sin(alpha)

    if (gradient == 0) then
      call track_a_drift (end%vec, length)
      return
    endif

    call offset_particle (ele, param, end, set$)

    f = gradient / (2 * e_start)            ! entrence kick
    end%vec(2) = end%vec(2) - f * end%vec(1)
    end%vec(4) = end%vec(4) - f * end%vec(3)

    f = gradient / (2 * sqrt_2 * cos_phi)
    x_pos = end%vec(1)
    y_pos = end%vec(3)
    end%vec(1) =  cos_a * end%vec(1)         + sin_a * end%vec(2) * e_start / f
    end%vec(2) = -sin_a * x_pos * f / e_end + cos_a * end%vec(2) * e_start / e_end
    end%vec(3) =  cos_a * end%vec(3)         + sin_a * end%vec(4) * e_start / f
    end%vec(4) = -sin_a * y_pos * f / e_end + cos_a * end%vec(4) * e_start / e_end

    f = gradient / (2 * e_end)              ! exit kick
    end%vec(2) = end%vec(2) + f * end%vec(1)
    end%vec(4) = end%vec(4) + f * end%vec(3)

    end%vec(6) = (e_end - ele%value(energy$)) / ele%value(energy$) 

    call offset_particle (ele, param, end, unset$)
    call end_z_calc

! sextupole
! The sextupole is treated as a drift with position dependent kick
! at the beginning and the end

  case (sextupole$)

    call offset_particle (ele, param, end, set$)

    k2l = ele%value(k2$) * length / (1 + end%vec(6))
    end%vec(2) = end%vec(2) + k2l * (end%vec(3)**2 - end%vec(1)**2)/4
    end%vec(4) = end%vec(4) + k2l * end%vec(1) * end%vec(3) / 2
    end%vec(1) = end%vec(1) + end%vec(2) * length
    end%vec(3) = end%vec(3) + end%vec(4) * length
    end%vec(2) = end%vec(2) + k2l * (end%vec(3)**2 - end%vec(1)**2)/4
    end%vec(4) = end%vec(4) + k2l * end%vec(1) * end%vec(3) / 2

    call offset_particle (ele, param, end, unset$)
    call end_z_calc

! solenoid

  case (solenoid$)

    call offset_particle (ele, param, end, set$)

    ks = ele%value(ks$) / (1 + end%vec(6))
    call solenoid_mat_calc (ks, length, mat4)
    end%vec(1:4) = matmul (mat4, end%vec(1:4))

    call offset_particle (ele, param, end, unset$)
    call end_z_calc

! sol_quad

  case (sol_quad$)

    call offset_particle (ele, param, end, set$)

    ks = ele%value(ks$) / (1 + end%vec(6))
    k1 = ele%value(k1$) / (1 + end%vec(6))
    vec0 = 0
    call sol_quad_mat6_calc (ks, k1, length, mat6, vec0)
    end%vec(1:4) = matmul (mat6(1:4,1:4), end%vec(1:4))

    call offset_particle (ele, param, end, unset$)
    call end_z_calc

! wiggler:

  case (wiggler$)

    if (ele%sub_key == map_type$) then
      print *, 'ERROR IN TRACK1_BMAD: NEW STYLE WIGGLER: ', ele%name
      print *, '       HAS TRACKING_METHOD = BMAD_STANDARD.'
      print *, '       THIS IS NOT A POSSIBLE OPTION FOR THE TRACKING_METHOD.'
      call err_exit
    endif

    call offset_particle (ele, param, end, set$, set_multipoles=.false.)

    k1 = ele%value(k1$) / (1 + end%vec(6))**2
    call quad_mat_calc (k1, length, mat2)
    end%vec(1) = end%vec(1) + length * end%vec(2)
    end%vec(3:4) = matmul (mat2, end%vec(3:4))

    call offset_particle (ele, param, end, unset$, set_multipoles=.false.)
    call end_z_calc

! multipole

  case (multipole$, ab_multipole$) 

    call offset_particle (ele, param, end, set$, &
                  set_canonical = .false., set_tilt = .false.)

    call multipole_ele_to_kt(ele, param%particle, knl, tilt, .true.)
    do n = 0, n_pole_maxx
      call multipole_kick (knl(n), tilt(n), n, end)
    enddo

    call offset_particle (ele, param, end, unset$, &
                  set_canonical = .false., set_tilt = .false.)

! accel_sol
! look at former routine in track1_mod:track_a_accel_sol ()

  case (accel_sol$)

    print *, 'ERROR: ACCEL_SOL MUST BE RESUSITATED!' 
    call err_exit

! unknown

  case default

    print *, 'ERROR IN TRACK1_BMAD: UNKNOWN ELEMENT: ', &
                                        key_name(ele%key), ele%type
    call err_exit

  end select


contains

!--------------------------------------------------------------
! Rough calculation for change in longitudinal position using:
!      dz = -L * (<p_x^2> + <p_y^2>)/ 2 
! where <...> means average.
! The formula below assumes a linear change in velocity between 
! the beginning and the end:

subroutine end_z_calc

  end%vec(5) = end%vec(5) - length * &
      (start%vec(2)**2 + end%vec(2)**2 + start%vec(2) * end%vec(2) + &
       start%vec(4)**2 + end%vec(4)**2 + start%vec(4) * end%vec(4)) / 6

end subroutine

end subroutine
