!+
! Subroutine track1_bmad (start, ele, param, end)
!
! Particle tracking through a single element BMAD_standard style.
! This routine is NOT ment for long term tracking since it does not get 
! all the 2nd order terms for the longitudinal motion.
!
! Modules Needed:
!   use bmad
!
! Input:
!   start  -- Coord_struct: Starting position
!   ele    -- Ele_struct: Element
!   param  -- Param_struct:
!     %aperture_limit_on -- If .true. then %LOST will be set if the
!                 particle is outsile the aperture.
!
! Output:
!   end   -- Coord_struct: End position
!   param
!     %lost -- Set .true. If the particle is outside the aperture and
!                %aperture_limit_on is set. Also: %lost is set .true. if
!                the particle does not make it through a bend irregardless
!                of the the setting of %aperture_limit_on.
!
! Notes:
!
! It is assumed that HKICK and VKICK are the kicks in the horizontal
! and vertical kicks irregardless of the value for TILT.
!
! TRACK1_BMAD *never* relies on ELE%MAT6 for tracking excect for 
! hybrid elements.
!-

#include "CESR_platform.inc"

subroutine track1_bmad (start, ele, param, end)

  use bmad

  implicit none

  type (coord_struct), intent(in)  :: start
  type (coord_struct), intent(out) :: end
  type (ele_struct),   intent(inout)  :: ele
  type (param_struct), intent(inout) :: param

  type (coord_struct)  c0
  type (ele_struct)  bend

  real(rdef) x_kick, y_kick, k1, k2l, k3l, length, phase, mat2(2,2), mat4(4,4)
  real(rdef) del, e1, e2, del_x_vel, del_y_vel, sig_x, sig_y, kx, ky, coef
  real(rdef) knl(0:n_pole_maxx), tilt(0:n_pole_maxx)
  real(rdef) ks, sig_x0, sig_y0, beta, mat6(6,6)
  real(rdef) z_slice(100), s_pos, s_pos_old, vec0(6)
  real(rdef) ave_x_vel2, ave_y_vel2

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

    call offset_particle (ele, param, end, set$)
    end%vec(1) = end%vec(1) + length * end%vec(2)
    end%vec(3) = end%vec(3) + length * end%vec(4)
    call offset_particle (ele, param, end, unset$)

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


    if (ele%value(k1$) /= 0) then
      call offset_particle (ele, param, end, set$)
      e1 = ele%value(e1$)
      e2 = ele%value(e2$)
      if (e1 /= 0) then
        del = tan(e1) * ele%value(g$) / (1 + end%z%vel)
        end%x%vel = end%x%vel + del * end%x%pos
        end%y%vel = end%y%vel - del * end%y%pos
      endif
      bend%value(k1$)       = ele%value(k1$) / (1 + end%z%vel)
      bend%value(g$)        = ele%value(g$)
      bend%value(g_design$) = ele%value(g_design$)
      bend%value(l$)        = ele%value(l$)
      call make_mat6(bend, param, c0, c0)
      end%vec = matmul(bend%mat6, end%vec)
      if (e2 /= 0) then
        del = tan(e2) * ele%value(g$) / (1 + end%z%vel)
        end%x%vel = end%x%vel + del * end%x%pos
        end%y%vel = end%y%vel - del * end%y%pos
      endif
      call offset_particle (ele, param, end, unset$)
    else
      call track_bend (end, ele, param, end)
      if (param%lost) return
    endif

! rfcavity

  case (rfcavity$)

    call offset_particle (ele, param, end, set$)

    end%vec(1) = end%vec(1) + length * end%vec(2) 
    end%vec(3) = end%vec(3) + length * end%vec(4) 

    if (ele%value(volt$) /= 0) then
      phase = ele%value(lag$) + end%z%pos / ele%value(rf_wavelength$)
      end%z%vel = end%z%vel + ele%value(volt$) * sin (twopi * phase) / &
                                                       (1e9 * param%energy)
    endif
         
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

  case (accel_sol$)

    print *, 'ERROR: ACCEL_SOL MUST BE RESUSITATED!' ! call track1_accel_sol ()
    call err_exit

! unknown

  case default

    type *, 'ERROR IN TRACK1_BMAD: UNKNOWN ELEMENT: ', &
                                        key_name(ele%key), ele%type
    call err_exit

  end select

!--------------------------------------------------------------
! Very crude calculation for change in longitudinal position

  if (ele%key /= sbend$ .and. param%lattice_type /= linac_lattice$ ) then
    if (start%x%vel*end%x%vel > 0) then  ! if same sign
      ave_x_vel2 = (start%x%vel**2 + end%x%vel**2) / 2
    else
      ave_x_vel2 = (start%x%vel**2 + end%x%vel**2) / 4
    endif
      
    if (start%y%vel*end%y%vel > 0) then  ! if same sign
      ave_y_vel2 = (start%y%vel**2 + end%y%vel**2) / 2
    else
      ave_y_vel2 = (start%y%vel**2 + end%y%vel**2) / 4
    endif

    end%z%pos = end%z%pos - length * (ave_x_vel2 + ave_y_vel2) / 2
  endif

end subroutine
