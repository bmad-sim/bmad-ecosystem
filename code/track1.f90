!+
! Subroutine TRACK1 (START, ELE, PARAM, END)
!
! Particle tracking through a single element. This routine is NOT ment for
! long term tracking since it does get all the 2nd order terms for 
! the longitudinal motion. For long term tracking see TRACK_LONG. 
!
! Modules Needed:
!   use bmad
!
! Input:
!   START  -- Coord_struct: Starting position
!   ELE    -- Ele_struct: Element
!   PARAM  -- Param_struct:
!     %APERTURE_LIMIT_ON -- If .true. then %LOST will be set if the
!                 particle is outsile the aperture.
!
! Output:
!   END   -- Coord_struct: End position
!   PARAM
!     %LOST -- Set .true. If the particle is outside the aperture and
!                %APERTURE_LIMIT_ON is set. Also: %LOST is set .true. if
!                the particle does not make it through a bend irregardless
!                of the the setting of %APERTURE_LIMIT_ON.
!
! Notes:
!
! It is assumed that HKICK and VKICK are the kicks in the horizontal
! and vertical kicks irregardless of the value for TILT.
!
! TRACK1 *never* relies on ELE%MAT6 for tracking excect for hybrid elements.
!-

!$Id$
!$Log$
!Revision 1.4  2002/02/23 20:32:25  dcs
!Double/Single Real toggle added
!
!Revision 1.3  2002/01/08 21:44:43  dcs
!Aligned with VMS version  -- DCS
!
!Revision 1.2  2001/09/27 18:31:58  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"

subroutine track1 (start, ele, param, end)

  use bmad

  implicit none

  type (coord_struct)  start, end, c0
  type (ele_struct)  ele, bend
  type (param_struct)  param

  real(rdef) x_kick, y_kick, k1, k2l, k3l, length, phase, mat2(2,2), mat4(4,4)
  real(rdef) del, e1, e2, del_x_vel, del_y_vel, sig_x, sig_y, kx, ky, coef
  real(rdef) knl(0:n_pole_maxx), tilt(0:n_pole_maxx)
  real(rdef) ks, sig_x0, sig_y0, beta, mat6(6,6)
  real(rdef) ave_x_vel2, ave_y_vel2, x_lim, y_lim
  real(rdef) z_slice(100), s_pos, s_pos_old, vec0(6)

  integer i, j, n, n_slice, key

  logical init_needed / .true. /

! init bend element                    

  if (init_needed) then
    call init_ele (bend)
    bend%key = sbend$
    init_needed = .false.
  endif

! runge_kutta tracking.
! track1_runge_kutta must be supplied by you if you want to do runge kutta
! tracking. 
! See the HTML bmad programming notes for more details.

  if (ele%tracking_method == runge_kutta$) then
    call track1_runge_kutta (start, ele, param, end)
    return
  elseif (ele%tracking_method == custom_calc$) then
    call custom_track1 (start, ele, param, end)
    return
  endif

! custom element

  if (ele%key == custom$) then
    call custom_track1 (start, ele, param, end)
    return
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

    call offset_particle (ele, param, end, set$)

    if (ele%value(k1$) /= 0) then
      e1 = ele%value(e1$)
      e2 = ele%value(e2$)
      if (e1 /= 0) then
        del = tan(e1) / (ele%value(rho$) * (1 + end%z%vel))
        end%x%vel = end%x%vel + del * end%x%pos
        end%y%vel = end%y%vel - del * end%y%pos
      endif
      bend%value(k1$)    = ele%value(k1$) / (1 + end%z%vel)
      bend%value(rho$) = ele%value(rho$)
      bend%value(rho_design$) = ele%value(rho_design$)
      bend%value(angle$) = ele%value(angle$)
      bend%value(l$)     = ele%value(l$)
      call make_mat6(bend, param, c0, c0)
      end%vec = matmul(bend%mat6, end%vec)
      if (e2 /= 0) then
        del = tan(e2) / (ele%value(rho$) * (1 + end%z%vel))
        end%x%vel = end%x%vel + del * end%x%pos
        end%y%vel = end%y%vel - del * end%y%pos
      endif
    else
      call track_bend (end, ele, end, param%lost)
      if (param%lost) return
    endif

    call offset_particle (ele, param, end, unset$)

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

! wiggler
! Note: k1 varies with 1/E^2, not 1/E as for a quad.
! Note: wiggler multipoles are per pole and thus cancel in the linear model.

  case (wiggler$)

    if (ele%tracking_method == linear$) then
      call offset_particle (ele, param, end, set$, set_multipoles=.false.)
      k1 = ele%value(k1$) / (1 + end%z%vel)**2
      call quad_mat_calc (k1, length, mat2)
      end%vec(1) = end%vec(1) + length * end%vec(2)
      end%vec(3:4) = matmul (mat2, end%vec(3:4))
      call offset_particle (ele, param, end, unset$, set_multipoles=.false.)
    else
      call track_wiggler (end, ele, param, end, param%lost)
      if (param%lost) return
    endif

! hybrid

  case (hybrid$)

    end%vec = matmul (ele%mat6, end%vec)

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

    print *, 'ERROR: ACCEL_SOL MUST BE RESUSITATED!' ! call track_accel_sol ()
    call err_exit

! redefinition of energy

  case (define_energy$)
    end%z%vel = param%energy / ele%value(energy$) * (end%z%vel + 1) - 1
    call make_mat6(ele, param)

! unknown

  case default

    type *, 'ERROR IN TRACK1: UNKNOWN ELEMENT: ', key_name(ele%key), ele%type
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

! check for particles outside aperture

  if (param%aperture_limit_on) then

    x_lim = ele%value(x_limit$)
    if (x_lim <= 0) x_lim = 1e10
    if (abs(end%x%pos) > x_lim) param%lost = .true.

    y_lim = ele%value(y_limit$)
    if (y_lim <= 0) y_lim = 1e10
    if (abs(end%y%pos) > y_lim) param%lost = .true.

  endif

end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine track_accel_sol (start, ele, param, end)
!
! Subroutine to track through an accel_sol element
!
! Modules Needed:
!   use bmad
!
! Input:
!   START  -- Coord_struct: Starting position
!   ELE    -- Ele_struct: Element
!   PARAM  -- Param_struct:
!     %APERTURE_LIMIT_ON -- If .true. then %LOST will be set if the
!                 particle is outsile the aperture.
!
! Output:
!   END   -- Coord_struct: End position
!   PARAM
!     %LOST -- Set .true. If the particle is outside the aperture and
!                %APERTURE_LIMIT_ON is set. Also: %LOST is set .true. if
!                the particle does not make it through a bend irregardless
!                of the the setting of %APERTURE_LIMIT_ON.
!-
                
subroutine track_accel_sol (start, ele, param, end) 

  use bmad

  implicit none

  type (coord_struct)  start, end, c0
  type (ele_struct)  ele
  type (param_struct)  param

  real(rdef) gamma_b, gamma_new, gamma_old, l_over_gamma, ll, ls(5), s_cumul
  real(rdef) s_grand_cum, vec_st(4), x_beg_lim, y_beg_lim, x_lim, y_lim
  real(rdef) x_lim_chng_rate, y_lim_chng_rate, phase, mat4(4,4), length
  real(rdef) b_x(5), b_y(5), beta_b, beta_s, c_e, c_m, en_gain, gam_inv2_b

  real(rdef), parameter :: beta_crit$ = 0.999

  integer i, j

  logical check_x_ap, check_y_ap, inside_segment

! linac stuff
! beta_b is the total speed in units of c_light (before entering the element)
! beta_s is the longitudinal speed in units of c_light

  gamma_b = param%energy * (end%z%vel + 1) / e_mass
  beta_b = sqrt(1 - 1 / gamma_b**2)
  gam_inv2_b = 1.0 / gamma_b**2
  if (gam_inv2_b <= 0.001) then
    beta_s = 1/sqrt(1 + end%x%vel**2 + end%y%vel**2) *(1-gam_inv2_b/2)
  else
    beta_s = sqrt((1 - gam_inv2_b) / (1 + end%x%vel**2 + end%y%vel**2))
  endif

! Calculation of end.z.vel

  end = start
  length = ele%value(l$)

  call offset_particle (ele, param, end, set$)

  if (ele%value(volt$) /= 0) then
    phase = ele%value(lag$) + end%z%pos  &
                      / (ele%value(rf_wavelength$) * beta_s)
    en_gain = ele%value(volt$) * sin(twopi * phase) / 1.e9
    if ((en_gain + gamma_b * e_mass) <= e_mass) then
      param%lost = .true.
      return
    else
      end%z%vel = end%z%vel + en_gain / param%energy
      c_e = en_gain / (e_mass * length)
    endif
  else
    c_e = 0.0
  endif

! Beginning fringe

  c_m = param%particle * c_light * ele%value(b_z$) / (e_mass * 1.e9)
  call mat_unit(mat4, 4, 4)
  mat4(2,3) = c_m / 2 *  &
                  sqrt((1 + end%x%vel**2 + end%y%vel**2) / (gamma_b**2 - 1))
  mat4(4,1) = -mat4(2,3)
  end%vec(1:4) = matmul(mat4, end%vec(1:4))

! Segment before first steerings:

  ls(1) = ele%value(s_st1$)
  b_x(1) = 0
  b_y(1) = 0

! Segment with first steerings:

  ls(2) = ele%value(l_st1$)
  b_x(2) = ele%value(b_x1$)
  b_y(2) = ele%value(b_y1$)

! Segment between steerings:

  ls(3) = ele%value(s_st2$) - (ele%value(s_st1$) + ele%value(l_st1$))
  b_x(3) = 0
  b_y(3) = 0

! Segment with second steerings:

  ls(4)= ele%value(l_st2$)
  b_x(4) = ele%value(b_x2$)
  b_y(4) = ele%value(b_y2$)

! Segment after second steerings:

  ls(5) = length - (ele%value(s_st2$) + ele%value(l_st2$))
  b_x(5) = 0
  b_y(5) = 0

  gamma_old = gamma_b
  gamma_new = gamma_old

  s_grand_cum = 0
  x_beg_lim = ele%value(x_beg_limit$)
  y_beg_lim = ele%value(y_beg_limit$)
  x_lim_chng_rate = (ele%value(x_limit$) - x_beg_lim) / length
  y_lim_chng_rate = (ele%value(y_limit$) - y_beg_lim) / length
  check_x_ap = .false.
  check_y_ap = .false.
  if (abs(c_m) > 0.001) then
    l_over_gamma = pi / (10 * abs(c_m))
    if (ele%value(x_limit$) * x_beg_lim /= 0) check_x_ap = .true.
    if (ele%value(y_limit$) * y_beg_lim /= 0) check_y_ap = .true.
  else
    l_over_gamma = 100 * pi
  endif

  do j = 1, 5

    s_cumul = 0
    inside_segment = .true.
    do while (inside_segment)
      gamma_old = gamma_new
      ll = l_over_gamma * gamma_old
      if (s_cumul + ll >= ls(j)) then
        ll = ls(j) - s_cumul
        inside_segment = .false.
      endif
      s_cumul = s_cumul + ll
      s_grand_cum = s_grand_cum + ll

      if (ll > 1.e-5) then
        gamma_new = gamma_old + c_e * ll
        call accel_sol_mat_calc (ll, c_m, c_e, gamma_old, gamma_new,  &
        b_x(j), b_y(j), end, mat4, vec_st)
        end%vec(1:4) = matmul(mat4, end%vec(1:4))
        do i = 1, 4
          end%vec(i) = end%vec(i) + param%particle * vec_st(i)
        enddo

! Calculation of end.z.pos

        if (abs(c_e) > 0.001) then
          end%z%pos = end%z%pos - ll + (sqrt(1 + (c_e * ll  &
                          + sqrt(gamma_old**2 - 1))**2) - gamma_old) / c_e
        else
          end%z%pos = end%z%pos + ll * (sqrt(1 - 1/gamma_old**2) - 1  &
              - c_e * ll * (gamma_old**2 - 0.5)/gamma_old**3)
        endif

        if (param%aperture_limit_on) then
          if (check_x_ap) then
            x_lim = x_lim_chng_rate * s_grand_cum + x_beg_lim
            if (abs(end%x%pos) > x_lim) param%lost = .true.
          endif
          if (check_y_ap) then
            y_lim = y_lim_chng_rate * s_grand_cum + y_beg_lim
            if (abs(end%y%pos) > y_lim) param%lost = .true.
          endif
          if (param%lost) return
        endif

      endif
    enddo

  enddo

! Ending fringe

  call mat_unit(mat4, 4, 4)
  mat4(4,1) = c_m / 2 *  &
                sqrt((1 + end%x%vel**2 + end%y%vel**2) / (gamma_new**2 - 1))
  mat4(2,3) = -mat4(4,1)
  end%vec(1:4) = matmul(mat4, end%vec(1:4))

  call offset_particle (ele, param, end, unset$)

end subroutine


!---------------------------------------------------------------------------
! old lattice stuff

!  if (param%lattice_type==linac_lattice$ .or. ele%key==accel_sol$) then
!
!    if (.not. ele%is_on .or. ele%key == drift$ .or. ele%key == elseparator$  &
!      .or. ele%key == kicker$) then
!      if (beta_b < beta_crit$) end%z%pos= end%z%pos + (beta_s - 1)*length
!    elseif (ele%key == solenoid$) then
!      if (beta_b < beta_crit$) end%z%pos= end%z%pos + (beta_s - 1)*length
!      ks = ele%value(ks$) / sqrt(1 + (end%z%vel**2  &
!             + 2 * end%z%vel)/(1 - (e_mass/param%energy)**2))
!    elseif (ele%key == sol_quad$) then
!      if (beta_b < beta_crit$) end%z%pos= end%z%pos + (beta_s - 1)*length
!      ks = ele%value(ks$) / sqrt(1 + (end%z%vel**2  &
!             + 2 * end%z%vel)/(1 - (e_mass/param%energy)**2))
!      k1 = ele%value(k1$) / sqrt(1 + (end%z%vel**2  &
!             + 2 * end%z%vel)/(1 - (e_mass/param%energy)**2))
!    elseif (ele%key == quadrupole$) then
!      if (beta_b < beta_crit$) end%z%pos= end%z%pos + (beta_s - 1)*length
!      k1 = ele%value(k1$) / sqrt(1 + (end%z%vel**2  &
!             + 2 * end%z%vel)/(1 - (e_mass/param%energy)**2))
!    elseif (ele%key == sextupole$) then
!      if (beta_b < beta_crit$) end%z%pos= end%z%pos + (beta_s - 1)*length
!      k2l = ele%value(k2$) * length / sqrt(1 + (end%z%vel**2  &
!             + 2 * end%z%vel)/(1 - (e_mass/param%energy)**2))
!    elseif (ele%key == octupole$) then
!      if (beta_b < beta_crit$) end%z%pos= end%z%pos + (beta_s - 1)*length
!      k3l = ele%value(k3$) * length / sqrt(1 + (end%z%vel**2  &
!              + 2 * end%z%vel)/(1 - (e_mass/param%energy)**2))
!    elseif (ele%key == wiggler$) then
!      if (beta_b < beta_crit$) end%z%pos= end%z%pos + (beta_s - 1)*length
!      k1 = ele%value(k1$) / sqrt(1 + (end%z%vel**2  &
!             + 2 * end%z%vel)/(1 - (e_mass/param%energy)**2))
!    elseif (ele%key == multipole$ .or. ele%key == ab_multipole$) then
!      if (beta_b < beta_crit$) end%z%pos= end%z%pos + (beta_s - 1)*length
!    endif
!
!  endif

