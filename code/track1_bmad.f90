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
  real(rp) ave_x_vel2, ave_y_vel2, rel_E
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

  case (drift$) 

    call track_a_drift (end%vec, length)
    return

! patch

  case (patch$)

    print *, 'ERROR IN TRACK1_BMAD: PATCH ELEMENT NOT YET IMPLEMENTED!'
    call err_exit

! kicker, separator

  case (elseparator$, kicker$) 

    call offset_particle (ele, param, end, set$)

    end%vec(1) = end%vec(1) + length * end%vec(2)
    end%vec(3) = end%vec(3) + length * end%vec(4)

    if (ele%key == kicker$) then
      end%x%pos = end%x%pos + ele%value(h_displace$)
      end%y%pos = end%y%pos + ele%value(v_displace$)
    endif

    call offset_particle (ele, param, end, unset$)  

! beambeam
                          
  case (beambeam$)

    if (ele%value(charge$) == 0 .or. param%n_part == 0) return

    call offset_particle (ele, param, end, set$)

    sig_x0 = ele%value(sig_x$)
    sig_y0 = ele%value(sig_y$)
    n_slice = max(1, nint(ele%value(n_slice$)))
    call bbi_slice_calc (n_slice, ele%value(sig_z$), z_slice)
    s_pos = 0    ! end at the ip
    do i = 1, n_slice
      s_pos_old = s_pos
      s_pos = (end%z%pos + z_slice(i)) / 2
      end%x%pos = end%x%pos + end%x%vel * (s_pos - s_pos_old)
      end%y%pos = end%y%pos + end%y%vel * (s_pos - s_pos_old)
      if (ele%x%beta == 0) then
        sig_x = sig_x0
        sig_y = sig_y0
      else
        beta = ele%x%beta - 2 * ele%x%alpha * s_pos + ele%x%gamma * s_pos**2
        sig_x = sig_x0 * sqrt(beta / ele%x%beta)
        beta = ele%y%beta - 2 * ele%y%alpha * s_pos + ele%y%gamma * s_pos**2
        sig_y = sig_y0 * sqrt(beta / ele%y%beta)
      endif

      call bbi_kick (end%x%pos/sig_x, end%y%pos/sig_y, sig_y/sig_x,  &
                                                                  kx, ky)
      coef = ele%value(bbi_const$) / (n_slice * (1 + end%z%vel))
      end%x%vel = end%x%vel + kx * coef
      end%y%vel = end%y%vel + ky * coef
    enddo
    end%x%pos = end%x%pos - end%x%vel * s_pos
    end%y%pos = end%y%pos - end%y%vel * s_pos

    call offset_particle (ele, param, end, unset$)  

! octupole
! The octupole is treated as a thin lens with a position dependent kick
! at the beginning and the end

  case (octupole$)

    call offset_particle (ele, param, end, set$)

    k3l = ele%value(k3$) * length / (1 + end%z%vel)

    end%x%vel = end%x%vel + k3l *  &
                    (3*end%x%pos*end%y%pos**2 - end%x%pos**3) / 12
    end%y%vel = end%y%vel + k3l *  &
                    (3*end%y%pos*end%x%pos**2 - end%y%pos**3) / 12

    end%x%pos = end%x%pos + end%x%vel * length
    end%y%pos = end%y%pos + end%y%vel * length

    end%x%vel = end%x%vel + k3l *  &
                    (3*end%x%pos*end%y%pos**2 - end%x%pos**3) / 12
    end%y%vel = end%y%vel + k3l *  &
                    (3*end%y%pos*end%x%pos**2 - end%y%pos**3) / 12

    call offset_particle (ele, param, end, unset$)  

! quadrupole

  case (quadrupole$)

    call offset_particle (ele, param, end, set$)

    k1 = ele%value(k1$) / (1 + end%z%vel)
    call quad_mat_calc (-k1, length, mat2)
    end%vec(1:2) = matmul(mat2, end%vec(1:2))
    call quad_mat_calc (k1, length, mat2)
    end%vec(3:4) = matmul(mat2, end%vec(3:4))

    call offset_particle (ele, param, end, unset$)  

! sbend
! A non-zero roll has a zeroth order effect that must be included

  case (sbend$)

    call track_a_bend (end, ele, param, end)
    return ! do not do z-calc at end of this routine

! rfcavity

  case (rfcavity$)

    call offset_particle (ele, param, end, set$, set_canonical = .false.)

    x = end%x%pos
    y = end%y%pos
    z = end%z%pos

    px = end%x%vel
    py = end%y%vel
    pz = end%z%vel


    phase = twopi * (ele%value(phi0$) + z / ele%value(rf_wavelength$))
    k  =  twopi * ele%value(volt$) * cos(phase) / &
                              (param%beam_energy * ele%value(rf_wavelength$))
    dE0 =  ele%value(volt$) * sin(phase) / param%beam_energy

    L = ele%value(l$)
    E = 1 + pz
    E2 = E**2
    pxy2 = px**2 + py**2

!

    end = start
    end%x%pos = x + px*L * (1/E - dE0/2 + pxy2*L/12 + pz*dE0 + dE0**2/3) 
    end%y%pos = y + py*L * (1/E - dE0/2 + pxy2*L/12 + pz*dE0 + dE0**2/3)
    end%z%pos = z + pxy2*L * (-1/(2*E2) + dE0/2)
    end%z%vel = pz + dE0 + k*pxy2*L * (-1/(4*E2) + dE0/6) 
         
    call offset_particle (ele, param, end, unset$, set_canonical = .false.)

! linac rf cavity
! assumes the particle is ultra-relativistic.
! This uses the formalism from:
!       J. Rosenzweig and L. Serafini
!       Phys Rev E, Vol. 49, p. 1599, (1994)
! with b_0 = b_-1 = 1

  case (lcavity$)

    if (ele%value(energy_start$) == 0) then
      print *, 'ERROR IN TRACK1_BMAD: REFERENCE BEAM ENERGY IS 0 FOR A LCAVITY!'
      call err_exit
    endif

    phase = twopi * (ele%value(phi0$) + &
                        end%z%pos * ele%value(rf_frequency$) / c_light)
    cos_phi = cos(phase)
    gradient = ele%value(gradient$) * cos_phi
    e_start = ele%value(energy_start$) * (1 + end%z%vel)
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
    end%x%vel = end%x%vel - f * end%x%pos
    end%y%vel = end%y%vel - f * end%y%pos

    f = gradient / (2 * sqrt_2 * cos_phi)
    x_pos = end%x%pos
    y_pos = end%y%pos
    end%x%pos =  cos_a * end%x%pos         + sin_a * end%x%vel * e_start / f
    end%x%vel = -sin_a * x_pos * f / e_end + cos_a * end%x%vel * e_start / e_end
    end%y%pos =  cos_a * end%y%pos         + sin_a * end%y%vel * e_start / f
    end%y%vel = -sin_a * y_pos * f / e_end + cos_a * end%y%vel * e_start / e_end

    f = gradient / (2 * e_end)              ! exit kick
    end%x%vel = end%x%vel + f * end%x%pos
    end%y%vel = end%y%vel + f * end%y%pos

    end%z%vel = (e_end - ele%value(energy$)) / ele%value(energy$) 

    call offset_particle (ele, param, end, unset$)

! sextupole
! The sextupole is treated as a drift with position dependent kick
! at the beginning and the end

  case (sextupole$)

    call offset_particle (ele, param, end, set$)

    k2l = ele%value(k2$) * length / (1 + end%z%vel)
    end%x%vel = end%x%vel + k2l * (end%y%pos**2 - end%x%pos**2)/4
    end%y%vel = end%y%vel + k2l * end%x%pos * end%y%pos / 2
    end%x%pos = end%x%pos + end%x%vel * length
    end%y%pos = end%y%pos + end%y%vel * length
    end%x%vel = end%x%vel + k2l * (end%y%pos**2 - end%x%pos**2)/4
    end%y%vel = end%y%vel + k2l * end%x%pos * end%y%pos / 2

    call offset_particle (ele, param, end, unset$)

! solenoid

  case (solenoid$)

    call offset_particle (ele, param, end, set$)

    ks = ele%value(ks$) / (1 + end%z%vel)
    call solenoid_mat_calc (ks, length, mat4)
    end%vec(1:4) = matmul (mat4, end%vec(1:4))

    call offset_particle (ele, param, end, unset$)

! sol_quad

  case (sol_quad$)

    call offset_particle (ele, param, end, set$)

    ks = ele%value(ks$) / (1 + end%z%vel)
    k1 = ele%value(k1$) / (1 + end%z%vel)
    vec0 = 0
    call sol_quad_mat6_calc (ks, k1, length, mat6, vec0)
    end%vec(1:4) = matmul (mat6(1:4,1:4), end%vec(1:4))

    call offset_particle (ele, param, end, unset$)

! wiggler:

  case (wiggler$)

    if (ele%sub_key == map_type$) then
      print *, 'ERROR IN TRACK1_BMAD: NEW STYLE WIGGLER: ', ele%name
      print *, '       HAS TRACKING_METHOD = BMAD_STANDARD.'
      print *, '       THIS IS NOT A POSSIBLE OPTION FOR THE TRACKING_METHOD.'
      call err_exit
    endif

    call offset_particle (ele, param, end, set$, set_multipoles=.false.)
    k1 = ele%value(k1$) / (1 + end%z%vel)**2
    call quad_mat_calc (k1, length, mat2)
    end%vec(1) = end%vec(1) + length * end%vec(2)
    end%vec(3:4) = matmul (mat2, end%vec(3:4))
    call offset_particle (ele, param, end, unset$, set_multipoles=.false.)

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

    type *, 'ERROR IN TRACK1_BMAD: UNKNOWN ELEMENT: ', &
                                        key_name(ele%key), ele%type
    call err_exit

  end select

!--------------------------------------------------------------
! Rough calculation for change in longitudinal position using:
!      dz = -L * (<p_x^2> + <p_y^2>)/ 2 
! where <...> means average.
! The formula below assumes a linear change in velocity between 
! the beginning and the end:

  end%z%pos = end%z%pos - length * &
      (start%x%vel**2 + end%x%vel**2 + start%x%vel * end%x%vel + &
       start%y%vel**2 + end%y%vel**2 + start%y%vel * end%y%vel) / 6

end subroutine
