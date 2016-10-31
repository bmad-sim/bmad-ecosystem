!+
! Subroutine track_a_bend (orbit, ele, param, mat6, make_matrix)
!
! Particle tracking through a bend element. 
!
! Modules Needed:
!   use bmad
!
! Input:
!   orbit       -- Coord_struct: Starting position.
!   ele         -- Ele_struct: Bend element.
!   param       -- lat_param_struct: Lattice parameters.
!   make_matrix -- logical, optional: Propagate the transfer matrix? Default is false.
!
! Output:
!   orbit      -- Coord_struct: End position.
!   mat6(6,6)  -- Real(rp), optional: Transfer matrix through the element.
!-

subroutine track_a_bend (orbit, ele, param, mat6, make_matrix)

use track1_mod, except_dummy => track_a_bend

implicit none

type (coord_struct) :: orbit, c0_off, c1_off, start_orb
type (ele_struct), target :: ele
type (lat_param_struct) :: param
type (fringe_edge_info_struct) fringe_info

real(rp), optional :: mat6(6,6)
real(rp) mat6_i(6,6)

real(rp) dp_long_dpx, dp_long_dpy, dp_long_dpz, dalpha_dpx, dalpha_dpy, dalpha_dpz
real(rp) Dy_dpy, Dy_dpz, dpx_t_dx, dpx_t_dpx, dpx_t_dpy, dpx_t_dpz
real(rp) df_dx, df_dpx, df_dpz, deps_dx, deps_dpx, deps_dpy, deps_dpz
real(rp) df_dpy, dbeta_dx, dbeta_dpx, dbeta_dpy, dbeta_dpz
real(rp) dfactor_dx, dfactor_dpx, dfactor_dpy, dfactor_dpz, factor1, factor2
real(rp) r_step, axis(3), w_mat(3,3), dr(3), s_off, mass, e_tot

real(rp) angle, ct, st, x, px, y, py, z, pz, dpx_t, p_long
real(rp) rel_p, rel_p2, Dy, px_t, factor, g, g_err, c_dir, stg
real(rp) step_len, g_tot, eps, pxy2, ff, fg, alpha, beta, one_ct
real(rp) k_1, k_x, x_c, om_x, om_y, tau_x, tau_y, arg, s_x, c_x, z_2, s_y, c_y, r(6)
real(rp) an(0:n_pole_maxx), bn(0:n_pole_maxx), an_elec(0:n_pole_maxx), bn_elec(0:n_pole_maxx)

integer n, n_step

logical, optional :: make_matrix
logical has_nonzero_pole, has_nonzero_elec, drifting

!-----------------------------------------------------------------------

start_orb = orbit

call offset_particle (ele, param, set$, orbit, set_multipoles = .false., set_hvkicks = .false.)
c0_off = orbit

! Entrance edge kick

c_dir = ele%orientation * orbit%direction * relative_tracking_charge(orbit, param)
if (logic_option(.false., make_matrix)) call mat_make_unit (mat6)

nullify(fringe_info%hard_ele)
fringe_info%particle_at = first_track_edge$
call apply_element_edge_kick(orbit, fringe_info, 0.0_rp, ele, param, .false., mat6, make_matrix)

! If we have a sextupole component then step through in steps of length ds_step

n_step = 1

if (is_true(ele%value(exact_multipoles$))) then
  call init_exact_bend_multipole_coefs (ele, param, .true., has_nonzero_pole) 
else
  call multipole_ele_to_ab(ele, .false., has_nonzero_pole, an,      bn,      magnetic$, include_kicks = .true.)
  call multipole_ele_to_ab(ele, .false., has_nonzero_elec, an_elec, bn_elec, electric$)
endif

if (ele%value(k2$) /= 0 .and. is_false(ele%value(exact_multipoles$))) then
  bn(2) = bn(2) + ele%value(k2$) * ele%value(l$) / 2
  has_nonzero_pole = .true.
endif

if (has_nonzero_pole .or. has_nonzero_elec) n_step = max(nint(ele%value(l$) / ele%value(ds_step$)), 1)

! Set some parameters

r_step = 1.0_rp / n_step
step_len = ele%value(l$) * r_step
g = ele%value(g$)
g_tot = (g + ele%value(g_err$)) * c_dir
g_err = g_tot - g
angle = g * step_len
pz = orbit%vec(6)
rel_p  = 1 + pz
rel_p2 = rel_p**2
k_1 = ele%value(k1$) * c_dir
if (is_true(ele%value(exact_multipoles$))) k_1 = 0  ! Is folded in with multipoles.
drifting = .false.

if (.not. ele%is_on) then
  g_err = -g
  g_tot = 0
  k_1 = 0
endif

! multipole kick at the beginning.

if (has_nonzero_pole .or. has_nonzero_elec) call apply_multipole_kicks (0.5_rp)

! And track with n_step steps

do n = 1, n_step

  ! with k1 /= 0 use small angle approximation

  if (k_1 /= 0) then

    call sbend_body_with_k1_map (ele, param, n_step, orbit, mat6, make_matrix)

  elseif (g == 0 .and. g_err == 0) then
    call track_a_drift (orbit, step_len, mat6, make_matrix)
    drifting = .true.

  !-----------------------------------------------------------------------
  ! Track through main body...
  ! Use Eqs (12.18) from Etienne Forest: Beam Dynamics.

  else
    ct = cos(angle)
    st = sin(angle)
    if (abs(angle) < 1d-7) then
      stg = step_len * (1 - angle**2 / 6)
      one_ct = step_len * angle / 2
    else
      stg = sin(angle) / g
      one_ct = (1 - ct) / g
    endif

    x  = orbit%vec(1)
    px = orbit%vec(2)
    y  = orbit%vec(3)
    py = orbit%vec(4)
    z  = orbit%vec(5)
   
    pxy2 = px**2 + py**2
    if (rel_p2 - pxy2 < 0.01) then  ! somewhat arbitrary cutoff
      orbit%state = lost$
      orbit%vec(1) = 2 * bmad_com%max_aperture_limit
      orbit%vec(3) = 2 * bmad_com%max_aperture_limit
      return
    endif 

    p_long = sqrt(rel_p2 - pxy2)

    ! The following is to make sure that a beam entering on-axis remains 
    ! *exactly* on-axis.

    if (pxy2 < 1d-5) then  
      ff = pxy2 / (2 * rel_p)
      fg = g * (pz - ff - ff*ff/2 - g_tot*x) - g_err
    else
      fg = g * p_long - g_tot * (1 + x * g)
    endif

    Dy  = sqrt(rel_p2 - py**2)
    px_t = px*ct + fg*stg
    dpx_t = -px*st*g + fg*ct

    if (abs(px) > Dy .or. abs(px_t) > Dy) then
      if (max(px, px_t) > 0) then; orbit%state = lost_pos_x_aperture$
      else;                        orbit%state = lost_neg_x_aperture$
      endif
      return
    endif

    ! Matrix

    if (logic_option(.false., make_matrix)) then

      call mat_make_unit (mat6_i)

      dp_long_dpx = -px/p_long
      dp_long_dpy = -py/p_long
      dp_long_dpz = rel_p/p_long

      if (pxy2 < 1d-5) then
         df_dx  = -g_tot
         df_dpx = -px * pxy2 / (2 * rel_p2) - px/rel_p
         df_dpy = -py * pxy2 / (2 * rel_p2) - py/rel_p
         df_dpz = 1 + pxy2**2 / (4 * rel_p**3) + pxy2 / (2 * rel_p2)
      else
         df_dx  = -g_tot
         df_dpx = dp_long_dpx
         df_dpy = dp_long_dpy
         df_dpz = dp_long_dpz
      endif

      Dy_dpy = -py/Dy
      Dy_dpz = rel_p/Dy

      dpx_t_dx  = ct*df_dx
      dpx_t_dpx = -st + ct*df_dpx
      dpx_t_dpy = ct*df_dpy
      dpx_t_dpz = ct*df_dpz

      if (abs(angle) < 1d-5 .and. abs(g_tot * step_len) < 1d-5) then
        mat6_i(1,1) = 1
        mat6_i(1,2) = step_len / p_long + step_len * px**2 / p_long**3 - 3 * g_tot * px * (step_len * Dy)**2 / (2 * p_long**5) + &
                      g * step_len * (step_len *px + x * (p_long - px**2 / p_long)) / p_long**2 + &
                      g * step_len * px * (step_len * (rel_p2 + px**2 - py**2) + 2 * x * px * p_long) / p_long**4
        mat6_i(1,3) = 0
        mat6_i(1,4) = step_len * px *py / p_long**3 + &
                      g_tot * step_len**2 * (py / p_long**3 - 3 * py * Dy**2 / (2 * p_long**5)) + &
                      g * step_len * (-step_len * py - x * px * py / p_long) / p_long**2 + &
                      g * step_len * (step_len * (rel_p2 + px**2 - py**2) + 2 * x * px * p_long) * py / p_long**4
        mat6_i(1,5) = 0
        mat6_i(1,6) = -step_len * px * rel_p / p_long**3 + &
                      g_tot * step_len**2 * (3 * rel_p * Dy**2 / (2 * p_long**5) - rel_p / p_long**3) + &
                      g * step_len * (step_len * rel_p + x * px * rel_p / p_long) / p_long**2 - &
                      g * step_len * (step_len * (rel_p2 + px**2 - py**2) + 2 * x * px * p_long) * rel_p / p_long**4

      elseif (abs(g_tot) < 1d-5 * abs(g)) then
        alpha = p_long * ct - px * st
        dalpha_dpx = dp_long_dpx * ct - st
        dalpha_dpy = dp_long_dpy * ct
        dalpha_dpz = dp_long_dpz * ct
        mat6_i(1,1) = -(g_tot*st**2*(1+g*x)*Dy**2)/(g*alpha**3) + p_long/alpha &
                      +(3*g_tot**2*st**3*(1+g*x)**2*Dy**2*(ct*px+st*p_long))/(2*g**2*alpha**5)
        mat6_i(1,2) = (3*g_tot*st**2*(1+g*x)**2*Dy**2*dalpha_dpx)/(2*g**2*alpha**4) &
                      -(5*g_tot**2*st**3*(1+g*x)**3*Dy**2*(ct*px+st*p_long)*dalpha_dpx)/(2*g**3*alpha**6) &
                      -((-alpha+(1+g*x)*p_long)*dalpha_dpx)/(g*alpha**2) &
                      +(g_tot**2*st**3*(1+g*x)**3*Dy**2*(ct+st*dp_long_dpx))/(2*g**3*alpha**5) &
                      +(-dalpha_dpx+(1+g*x)*dp_long_dpx)/(g*alpha)
        mat6_i(1,4) = (3*g_tot*st**2*(1+g*x)**2*Dy**2*dalpha_dpy)/(2*g**2*alpha**4) &
                      -(5*g_tot**2*st**3*(1+g*x)**3*Dy**2*(ct*px+st*p_long)*dalpha_dpy)/(2*g**3*alpha**6) &
                      -((-alpha+(1+g*x)*p_long)*dalpha_dpy)/(g*alpha**2) &
                      +(g_tot**2*st**4*(1+g*x)**3*Dy**2*dp_long_dpy)/(2*g**3*alpha**5) &
                      +(-dalpha_dpy+(1+g*x)*dp_long_dpy)/(g*alpha) &
                      -(g_tot*st**2*(1+g*x)**2*Dy*Dy_dpy)/(g**2*alpha**3) &
                      +(g_tot**2*st**3*(1+g*x)**3*Dy*(ct*px+st*p_long)*Dy_dpy)/(g**3*alpha**5)
        mat6_i(1,6) = (3*g_tot*st**2*(1+g*x)**2*Dy**2*dalpha_dpz)/(2*g**2*alpha**4) &
                      -(5*g_tot**2*st**3*(1+g*x)**3*Dy**2*(ct*px+st*p_long)*dalpha_dpz)/(2*g**3*alpha**6) &
                      -((-alpha+(1+g*x)*p_long)*dalpha_dpz)/(g*alpha**2) &
                      +(g_tot**2*st**4*(1+g*x)**3*Dy**2*dp_long_dpz)/(2*g**3*alpha**5) &
                      +(-dalpha_dpz+(1+g*x)*dp_long_dpz)/(g*alpha) &
                      -(g_tot*st**2*(1+g*x)**2*Dy*Dy_dpz)/(g**2*alpha**3) &
                      +(g_tot**2*st**3*(1+g*x)**3*Dy*(ct*px+st*p_long)*Dy_dpz)/(g**3*alpha**5)
      else
        eps = px_t**2 + py**2
        deps_dx  = 2*px_t*st*df_dx
        deps_dpx = 2*px_t*(ct+st*df_dpx)
        deps_dpy = 2*px_t*st*df_dpy + 2*py
        deps_dpz = 2*px_t*st*df_dpz
!        if (eps < 1d-5 * rel_p2 ) then  ! use small angle approximation
!          eps = eps / (2 * rel_p)
!          deps_dx  = deps_dx / (2 * rel_p)
!          deps_dpx = deps_dpx / (2 * rel_p)
!          deps_dpy = deps_dpy / (2 * rel_p)
!          deps_dpz = deps_dpz / (2 * rel_p) - (px_t**2 + py**2) / (2*rel_p2) 
!          mat6_i(1,1) = (-dpx_t_dx + (eps/(2*rel_p) - 1)*deps_dx + eps*deps_dx/(2*rel_p))/g_tot
!          mat6_i(1,2) = (-dpx_t_dpx + (eps/(2*rel_p) - 1)*deps_dpx + eps*deps_dpx/(2*rel_p))/g_tot
!          mat6_i(1,4) = (-dpx_t_dpy + (eps/(2*rel_p) - 1)*deps_dpy + eps*deps_dpy/(2*rel_p))/g_tot
!          mat6_i(1,6) = (1 - dpx_t_dpz + (eps/(2*rel_p) - 1)*deps_dpz + eps*(deps_dpz/(2*rel_p) - eps/(2*rel_p2)))/g_tot
!        else
          mat6_i(1,1) = (-dpx_t_dx - deps_dx/(2*sqrt(rel_p2 - eps)))/g_tot
          mat6_i(1,2) = (-dpx_t_dpx - deps_dpx/(2*sqrt(rel_p2 - eps)))/g_tot
          mat6_i(1,4) = (-dpx_t_dpy - deps_dpy/(2*sqrt(rel_p2 - eps)))/g_tot
          mat6_i(1,6) = (-dpx_t_dpz + (2*rel_p - deps_dpz)/(2*sqrt(rel_p2 - eps)))/g_tot
!        endif
      endif
      
      mat6_i(2,1) = -g_tot * st
      mat6_i(2,2) = ct - px * st / p_long
      mat6_i(2,4) = -py * st / p_long
      mat6_i(2,6) = rel_p * st / p_long

      if (abs(g_tot) < 1d-5 * abs(g)) then
        beta = (1 + g * x) * st / (g * alpha) - &
               g_tot * (px * ct + p_long * st) * (st * (1 + g * x))**2 / (2 * g**2 * alpha**3)
        dbeta_dx  = st/alpha - (g_tot*st**2*(1+g*x)*(ct*px+st*p_long))/(g*alpha**3)
        dbeta_dpx = -(st*(1+g*x)*dalpha_dpx)/(g*alpha**2)-(g_tot*st**2*(1+g*x)**2*(ct+st*dp_long_dpx))/(2*g**2*alpha**3)
        dbeta_dpy = -(st*(1+g*x)*dalpha_dpy)/(g*alpha**2)-(g_tot*st**3*(1+g*x)**2*dp_long_dpy)/(2*g**2*alpha**3)
        dbeta_dpz = -(st*(1+g*x)*dalpha_dpz)/(g*alpha**2)-(g_tot*st**3*(1+g*x)**2*dp_long_dpz)/(2*g**2*alpha**3)
        mat6_i(3,1) = py*dbeta_dx
        mat6_i(3,2) = py*dbeta_dpx
        mat6_i(3,4) = beta + py*dbeta_dpy
        mat6_i(3,6) = py*dbeta_dpz
        mat6_i(5,1) = -rel_p*dbeta_dx
        mat6_i(5,2) = -rel_p*dbeta_dpx
        mat6_i(5,4) = -rel_p*dbeta_dpy
        mat6_i(5,6) = -beta - rel_p*dbeta_dpz
      else
        factor = (asin(px/Dy) - asin(px_t/Dy)) / g_tot
        factor1 = sqrt(1-(px/Dy)**2)
        factor2 = sqrt(1-(px_t/Dy)**2)
        dfactor_dx  = -st*df_dx/(Dy*factor2*g_tot)
        dfactor_dpx = (1/(factor1*Dy)-(ct+st*df_dpx)/(factor2*Dy))/g_tot
        dfactor_dpy = (-px*Dy_dpy/(factor1*Dy**2)-(-px_t*Dy_dpy/Dy**2 + st*df_dpy/Dy)/factor2)/g_tot
        dfactor_dpz = (-px*Dy_dpz/(factor1*Dy**2)-(-px_t*Dy_dpz/Dy**2 + st*df_dpz/Dy)/factor2)/g_tot
        mat6_i(3,1) = py*dfactor_dx
        mat6_i(3,2) = py*dfactor_dpx
        mat6_i(3,4) = angle/g_tot + factor + py*dfactor_dpy
        mat6_i(3,6) = py*dfactor_dpz
        mat6_i(5,1) = -rel_p*dfactor_dx
        mat6_i(5,2) = -rel_p*dfactor_dpx
        mat6_i(5,4) = -rel_p*dfactor_dpy
        mat6_i(5,6) = -angle/g_tot - factor - rel_p*dfactor_dpz
      endif

      mat6 = matmul(mat6_i, mat6)
    endif

    !

    if (abs(angle) < 1d-5 .and. abs(g_tot * step_len) < 1d-5) then
      orbit%vec(1) = orbit%vec(1) + step_len * px / p_long - &
                       g_tot * (step_len * Dy)**2 / (2 * p_long**3) + &
                       g * step_len * (step_len * (rel_p2 + px**2 - py**2) + 2 * x * px * p_long) / (2 * p_long**2)
    elseif (abs(g_tot) < 1d-5 * abs(g)) then
      alpha = p_long * ct - px * st
      orbit%vec(1) = (p_long * (1 + g * x) - alpha) / (g * alpha) - &
                   g_tot * (Dy * (1 + g * x) * st)**2 / (2 * alpha**3 * g**2) + &
                   g_tot**2 * Dy**2 * ((1 + g * x) * st)**3 * (px * ct + p_long * st) / (2 * alpha**5 * g**3)
    else
      eps = px_t**2 + py**2
!      if (eps < 1d-5 * rel_p2) then  ! use small angle approximation
!        eps = eps / (2 * rel_p)
!        orbit%vec(1) = (pz + px * st - ct * p_long + g_tot * x * ct + ct - g_err * one_ct + &
!                                                                 eps * (eps / (2 * rel_p) - 1)) / g_tot
!      else
        orbit%vec(1) = (sqrt(rel_p2 - eps) + px*st + g_tot*x*ct - p_long*ct) / g_tot - one_ct
!      endif
    endif

    orbit%vec(2) = px_t
    orbit%vec(4) = py

    if (abs(g_tot) < 1d-5 * abs(g)) then
      beta = (1 + g * x) * st / (g * alpha) - &
             g_tot * (px * ct + p_long * st) * (st * (1 + g * x))**2 / (2 * g**2 * alpha**3)
      orbit%vec(3) = y + py * beta
      orbit%vec(5) = z + step_len  - rel_p * beta 
    else
      factor = (asin(px/Dy) - asin(px_t/Dy)) / g_tot
      orbit%vec(3) = y + py * (angle/g_tot + factor)
      orbit%vec(5) = z + step_len * (g_err - g*pz) / g_tot - rel_p * factor
    endif

  endif

  ! multipole kick

  if (has_nonzero_pole .or. has_nonzero_elec) then
    if (n == n_step) then
      call apply_multipole_kicks (0.5_rp)
    else
      call apply_multipole_kicks (1.0_rp)
    endif
  endif

enddo

! Track through the exit face. Treat as thin lens.
! Need low energy z correction except when using track_a_drift.

if (orbit_too_large(orbit, param)) return

fringe_info%particle_at = second_track_edge$
call apply_element_edge_kick(orbit, fringe_info, 0.0_rp, ele, param, .false., mat6, make_matrix)

c1_off = orbit
call offset_particle (ele, param, unset$, orbit, set_multipoles = .false., set_hvkicks = .false.)

if (.not. drifting) call track1_low_energy_z_correction (orbit, ele, param)

orbit%t = start_orb%t + ele%value(delta_ref_time$) + (start_orb%vec(5) - orbit%vec(5)) / (orbit%beta * c_light)

if (orbit%direction == 1) then
  orbit%s = ele%s
else
  orbit%s = ele%s - ele%value(l$)
endif

! matrix

if (logic_option(.false., make_matrix)) then

  ! Roll
  ! c0_off is the coordinates *after* the roll at the entrance end
  ! So get the reverse roll matrix and take the inverse.

  if (ele%value(roll$) /= 0) then
    dr = 0
    if (ele%value(angle$) < 1d-20) then
      axis = [ele%value(angle$)/2, 0.0_rp, 1.0_rp]
    else
      axis = [cos(ele%value(angle$)) - 1, 0.0_rp, sin(ele%value(angle$))]
    endif
    call axis_angle_to_w_mat (axis, -ele%value(roll$), w_mat)
    call mat6_coord_transformation (mat6_i, ele, param, c0_off, dr, w_mat)
    mat6_i = mat_symp_conj(mat6_i)   ! Inverse
    mat6 = matmul(mat6, mat6_i)

    ! c1_off is the coordinates before the roll so this is what is needed
    axis(1) = -axis(1)  ! Axis in exit coordinates
    call axis_angle_to_w_mat (axis, -ele%value(roll$), w_mat)
    call mat6_coord_transformation (mat6_i, ele, param, c1_off, dr, w_mat)
    mat6 = matmul(mat6_i, mat6)

  endif

  !

  if (ele%value(ref_tilt_tot$) /= 0) call tilt_mat6 (mat6, ele%value(ref_tilt_tot$))

  if (ele%value(z_offset_tot$) /= 0) then
    s_off = ele%value(z_offset_tot$) * ele%orientation
    mat6(1,:) = mat6(1,:) - s_off * mat6(2,:)
    mat6(3,:) = mat6(3,:) - s_off * mat6(4,:)
    mat6(:,2) = mat6(:,2) + mat6(:,1) * s_off
    mat6(:,4) = mat6(:,4) + mat6(:,3) * s_off
  endif

  call mat6_add_pitch (ele%value(x_pitch_tot$), ele%value(y_pitch_tot$), ele%orientation, ele%mat6)

  ! 1/gamma^2 m56 correction

  if (.not. drifting) then
    mass = mass_of(orbit%species)
    e_tot = ele%value(p0c$) * (1 + orbit%vec(6)) / orbit%beta
    mat6(5,6) = mat6(5,6) + ele%value(l$) * mass**2 * ele%value(e_tot$) / e_tot**3
  endif

  ele%vec0 = orbit%vec - matmul(mat6, start_orb%vec)
endif

!-------------------------------------------------------------------------------------------------------
contains

subroutine apply_multipole_kicks (coef)

type (em_field_struct) field
type (em_potential_struct) potential

real(rp) coef, dmat6(6,6)
real(rp) rel_p, beta0, dt_ds_ref, p0, e_tot, direction, charge
real(rp) db_coef, dz_coef, vel(3), v2, E_force(3), B_force(3), f_bend, dt_ds, dp_ds, dbeta_ds, pz_p0
real(rp) dbeta_dpz, dE_dpz, dvel_dpx(3), dvel_dpy(3), dvel_dpz(3), dt_ds_dx, dt_ds_dpx, dt_ds_dpy, dt_ds_dpz
real(rp) dp_ds_dx, dp_ds_dpx, dp_ds_dy, dp_ds_dpy, dp_ds_dpz, cc, beta

integer i

!

if (is_true(ele%value(exact_multipoles$)) .and. ele%value(g$) /= 0) then

  rel_p = 1 + orbit%vec(6)
  beta = orbit%beta
  beta0 = ele%value(p0c$) / ele%value(e_tot$) 
  dt_ds_ref = orbit%direction / (beta0 * c_light)
  p0 = ele%value(p0c$) / c_light
  e_tot = orbit%p0c * rel_p / beta
  direction = ele%orientation * orbit%direction
  charge = charge_of(orbit%species)

  ! Calculate the field. 
  ! Important: Field is in frame of element. When ele%orientation = -1 => +z in -s direction.

  call  exact_bend_multipole_field (ele, param, orbit, .true., field, potential, make_matrix)

  ! Bend factor

  db_coef = mass_of(orbit%species)**2 * c_light / e_tot**3

  vel(1:2) = [orbit%vec(2), orbit%vec(4)] / rel_p
  v2 = vel(1)**2 + vel(2)**2
  if (v2 > 0.99999999_rp) return
  vel = beta * c_light * [vel(1), vel(2), sqrt(1 - v2) * direction]
  field%E = charge * field%E
  field%B = charge * field%B
  E_force = field%E
  B_force = cross_product(vel, field%B)

  f_bend = 1 + orbit%vec(1) * ele%value(g$) ! Longitudinal distance per unit s-distance. 
  dt_ds = orbit%direction * f_bend / abs(vel(3))
  dp_ds = dot_product(E_force, vel) * dt_ds / (beta * c_light)
  dbeta_ds = db_coef * dp_ds
  pz_p0 = rel_p * orbit%direction * abs(vel(3)) / (beta * c_light)  ! Pz / P0

  if (logic_option(.false., make_matrix)) then
    cc = ele%value(g$) * orbit%direction / c_light
    field%dE = charge * field%dE
    field%dB = charge * field%dB

    dbeta_dpz = (1 - beta**2) / e_tot
    dE_dpz = orbit%p0c * beta

    dvel_dpx = [1.0_rp, 0.0_rp, -vel(1) / vel(3)] * (beta * c_light / rel_p)
    dvel_dpy = [0.0_rp, 1.0_rp, -vel(2) / vel(3)] * (beta * c_light / rel_p)
    dvel_dpz = [-vel(1) * dE_dpz / E_tot, -vel(2) * dE_dpz / E_tot, (c_light * dbeta_dpz - ((vel(1)  + vel(2)) * dE_dpz / E_tot)) / vel(3)]

    dt_ds_dx  = orbit%direction * ele%value(g$) / abs(vel(3))
    dt_ds_dpx = -orbit%direction * f_bend * dvel_dpx(3) * sign_of(vel(3)) / vel(3)**2
    dt_ds_dpy = -orbit%direction * f_bend * dvel_dpy(3) * sign_of(vel(3)) / vel(3)**2
    dt_ds_dpz = -orbit%direction * f_bend * dvel_dpz(3) * sign_of(vel(3)) / vel(3)**2

    dp_ds_dx  = dot_product(field%dE(:,1), vel) * dt_ds / (beta * c_light) + dot_product(E_force, vel) * dt_ds_dx / (beta * c_light)
    dp_ds_dpx = (dot_product(E_force, vel) * dt_ds_dpx + dot_product(E_force, dvel_dpx) * dt_ds) / (beta * c_light)
    dp_ds_dy  = dot_product(field%dE(:,2), vel) * dt_ds / (beta * c_light)
    dp_ds_dpy = (dot_product(E_force, vel) * dt_ds_dpy + dot_product(E_force, dvel_dpy) * dt_ds) / (beta * c_light)
    dp_ds_dpz = (dot_product(E_force, vel) * dt_ds_dpz + dot_product(E_force, dvel_dpz) * dt_ds) / (beta * c_light) - dp_ds * dbeta_dpz / beta

    dmat6 = 0

    dmat6(2,1) = (field%dE(1,1) + cross(vel, field%dB(:,1), 1)) * dt_ds / p0 + (E_force(1) + B_force(1)) * dt_ds_dx / p0
    dmat6(2,2) = (E_force(1) + B_force(1)) * dt_ds_dpx / p0 + cross(dvel_dpx, field%B, 1) * dt_ds / p0 
    dmat6(2,3) = (field%dE(1,2) + cross(vel, field%dB(:,2), 1)) * dt_ds / p0
    dmat6(2,4) = (E_force(1) + B_force(1)) * dt_ds_dpy / p0 + cross(dvel_dpy, field%B, 1) * dt_ds / p0 
    dmat6(2,6) = (E_force(1) + B_force(1)) * dt_ds_dpz / p0 + cross(dvel_dpz, field%B, 1) * dt_ds / p0 

    dmat6(4,1) = (field%dE(2,1) + cross(vel, field%dB(:,1), 2)) * dt_ds / p0 + (E_force(2) + B_force(2)) * dt_ds_dx / p0
    dmat6(4,2) = (E_force(2) + B_force(2)) * dt_ds_dpx / p0 + cross(dvel_dpx, field%B, 2) * dt_ds / p0
    dmat6(4,3) = (field%dE(2,2) + cross(vel, field%dB(:,2), 2)) * dt_ds / p0
    dmat6(4,4) = (E_force(2) + B_force(2)) * dt_ds_dpy / p0 + cross(dvel_dpy, field%B, 2) * dt_ds / p0
    dmat6(4,6) = (E_force(2) + B_force(2)) * dt_ds_dpz / p0 + cross(dvel_dpz, field%B, 2) * dt_ds / p0

    dmat6(5,1) = db_coef * orbit%vec(5) * dp_ds_dx / beta
    dmat6(5,2) = db_coef * orbit%vec(5) * dp_ds_dpx / beta
    dmat6(5,3) = db_coef * orbit%vec(5) * dp_ds_dy / beta
    dmat6(5,4) = db_coef * orbit%vec(5) * dp_ds_dpy / beta
    dmat6(5,6) = db_coef * orbit%vec(5) * dp_ds_dpz / beta - dbeta_ds * orbit%vec(5) * (dbeta_dpz/beta + 3 * dE_dpz/e_tot) / beta

    dmat6(6,1) = dp_ds_dx / p0
    dmat6(6,2) = dp_ds_dpx / p0
    dmat6(6,3) = dp_ds_dy / p0
    dmat6(6,4) = dp_ds_dpy / p0
    dmat6(6,6) = dp_ds_dpz / p0

    dmat6 = dmat6 * (coef * step_len)
    forall (i = 1:6) dmat6(i,i) = 1 + dmat6(i,i)

    mat6 = matmul(dmat6, mat6)

  endif

  cc = coef * step_len
  orbit%vec(2) = orbit%vec(2) + cc * (E_force(1) + B_force(1)) * dt_ds / p0
  orbit%vec(4) = orbit%vec(4) + cc * (E_force(2) + B_force(2)) * dt_ds / p0
  orbit%vec(5) = orbit%vec(5) + cc * dbeta_ds * orbit%vec(5) / beta
  orbit%vec(6) = orbit%vec(6) + cc * dp_ds / p0

!

else
  if (has_nonzero_pole) call ab_multipole_kicks (an,      bn,      param%particle, orbit, magnetic$, coef * r_step,   mat6, make_matrix)
  if (has_nonzero_elec) call ab_multipole_kicks (an_elec, bn_elec, param%particle, orbit, electric$, coef * step_len, mat6, make_matrix)
endif

end subroutine apply_multipole_kicks

!---------------------------------------------------------------------------
! contains

function cross(A, B, ix) result (C_ix)

real(rp) A(3), B(3), C(3), C_ix
integer ix

C = cross_product(A, B)
C_ix = C(ix)

end function cross

end subroutine track_a_bend

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine sbend_body_with_k1_map (ele, param, n_step, orbit, mat6, make_matrix)
!
! Subroutine to calculate for a single step the transfer matrix and/or 
! ending coordinates for a sbend with a finite k1 but without a tilt.
!
! Modules Needed:
!   use bmad
!
! Input:
!   ele          -- Ele_struct: Sbend element.
!   param        -- Lat_param_struct: Branch parameters.
!   n_step       -- Integer: Number of steps to divide the bend into.
!                     Only one step is taken by this routine.
!   orbit        -- coord_struct: Orbit at beginning of the bend.
!   mat6(6,6)    -- Real(rp), optional: Transfer matrix before element.
!   make_matrix  -- logical, optional: Propagate the transfer matrix? Default is false.
!
! Output:
!   orbit      -- coord_struct: Ending coordinates.
!   mat6(6,6)  -- Real(rp), optional: Transfer matrix with body added in.
!-

subroutine sbend_body_with_k1_map (ele, param, n_step, orbit, mat6, make_matrix)

use bmad, except_dummy => sbend_body_with_k1_map

implicit none

type (ele_struct) ele
type (lat_param_struct) param
type (coord_struct) orbit

real(rp), optional :: mat6(6,6)
real(rp) mat6_i(6,6)
real(rp) g, g_err, length
real(rp) k_1, k_x, x_c, om_x, om_y, tau_x, tau_y, arg, s_x, c_x, s_y, c_y, r(6)
real(rp) z0, z1, z2, z11, z12, z22, z33, z34, z44
real(rp) dom_x, dom_xx, dx_c, dc_x, ds_x, dom_y, dom_yy, dc_y, ds_y, dcs_x, dcs_y
real(rp) g_tot, rel_p, rel_p2, charge_dir
real(rp) rel_pc, px, py, pxy2, pz

integer n_step
logical, optional :: make_matrix

! Degenerate case

charge_dir = relative_tracking_charge(orbit, param) * ele%orientation * orbit%direction

k_1 = ele%value(k1$) * charge_dir
g = ele%value(g$)
g_tot = (g + ele%value(g_err$)) * charge_dir
g_err = g_tot - g
length = ele%value(l$) / n_step

!

g_tot = g + g_err
rel_p = (1 + orbit%vec(6))
rel_p2 = rel_p**2


k_x = k_1 + g * g_tot
x_c = (g * rel_p - g_tot) / k_x

om_x = sqrt(abs(k_x) / rel_p)
om_y = sqrt(abs(k_1) / rel_p)

tau_x = -sign (1.0_rp, k_x)
tau_y =  sign (1.0_rp, k_1)

arg = om_x * length
if (arg < 1d-6) then
  s_x = (1 + tau_x * arg**2 / 6) * length
  c_x = 1 + tau_x * arg**2 / 2
  z2 = g * length**2 / (2 * rel_p)
elseif (k_x > 0) then
  s_x = sin(arg) / om_x
  c_x = cos(arg)
  z2 = tau_x * g * (1 - c_x) / (rel_p * om_x**2)
else
  s_x = sinh(arg) / om_x
  c_x = cosh(arg)
  z2 = tau_x * g * (1 - c_x) / (rel_p * om_x**2)
endif

arg = om_y * length
if (arg < 1d-6) then
  s_y = (1 + tau_y * arg**2 / 6) * length
  c_y = 1 + tau_y * arg**2 / 2
elseif (k_1 < 0) then
  s_y = sin(om_y * length) / om_y
  c_y = cos(om_y * length)
else
  s_y = sinh(om_y * length) / om_y
  c_y = cosh(om_y * length)
endif

r = orbit%vec
r(1) = r(1) - x_c

!

z0  = -g * x_c * Length
z1  = -g * s_x
z11 = tau_x * om_x**2 * (length - c_x * s_x) / 4
z12 = -tau_x * om_x**2 * s_x**2 / (2 * rel_p) 
z22 = -(length + c_x * s_x) / (4 * rel_p2) 
z33 = tau_y * om_y**2 * (length - c_y * s_y) / 4
z34 = -tau_y * om_y**2 * s_y**2 / (2 * rel_p) 
z44 = -(length + c_y * s_y) / (4 * rel_p2)

! Jacobian matrix

if (logic_option(.false., make_matrix)) then

  dom_x = -om_x / (2 * rel_p)
  dom_xx = -1 / (2 * rel_p)  ! [d(om_x) / d(p_z)] / om_x
  dx_c = g / k_x
  dc_x = tau_x * s_x * om_x * dom_x * length
  ds_x = (c_x * length - s_x) * dom_xx 

  dom_y = -om_y / (2 * rel_p)
  dom_yy = -1 / (2 * rel_p)  ! [d(om_y) / d(p_z)] / om_y
  dc_y = tau_y * s_y * om_y * dom_y * length
  ds_y = (c_y * length - s_y) * dom_yy

  dcs_x = c_x * ds_x + dc_x * s_x
  dcs_y = c_y * ds_y + dc_y * s_y

  mat6_i = 0

  mat6_i(1,1) = c_x
  mat6_i(1,2) = s_x / rel_p
  mat6_i(1,6) = dx_c * (1 - c_x) + dc_x * r(1) + &
              ds_x * r(2) / rel_p - s_x * r(2) / rel_p2
  mat6_i(2,1) = tau_x * om_x**2 * rel_p * s_x
  mat6_i(2,2) = c_x
  mat6_i(2,6) = tau_x * r(1) * 2 * om_x * dom_x * rel_p * s_x + &
                tau_x * r(1) * om_x**2 * s_x + &
                tau_x * r(1) * om_x**2 * rel_p * ds_x - &
                tau_x * dx_c * om_x**2 * rel_p * s_x + dc_x * r(2)

  mat6_i(3,3) = c_y
  mat6_i(3,4) = s_y / rel_p
  mat6_i(3,6) = dc_y * r(3) + ds_y * r(4) / rel_p - s_y * r(4) / rel_p2
  mat6_i(4,3) = tau_y * om_y**2 * rel_p * s_y
  mat6_i(4,4) = c_y
  mat6_i(4,6) = tau_y * r(3) * 2 * om_y * dom_y * rel_p * s_y + &
                tau_y * r(3) * om_y**2 * s_y + &
                tau_y * r(3) * om_y**2 * rel_p * ds_y + &
                dc_y * r(4)

  mat6_i(5,1) = z1 + 2 * z11 * r(1) +     z12 * r(2)  
  mat6_i(5,2) = z2 +     z12 * r(1) + 2 * z22 * r(2)
  mat6_i(5,3) =      2 * z33 * r(3) +     z34 * r(4)  
  mat6_i(5,4) =          z34 * r(3) + 2 * z44 * r(4)
  mat6_i(5,5) = 1
  mat6_i(5,6) = -dx_c * (z1 + 2 * z11 * r(1) + z12 * r(2)) - & 
                g * length * dx_c - g * ds_x * r(1) - &           ! dz0 & dz1
                tau_x * g * dc_x * r(2) / (om_x**2 * rel_p) - &   ! dz2
                (z11 / rel_p + tau_x * om_x**2 * dcs_x / 4) * r(1)**2 - &
                (2 * z12 / rel_p + tau_x * om_x**2 * s_x * ds_x / rel_p) * r(1) * r(2) - &
                (2 * z22 / rel_p + dcs_x / (4 * rel_p2)) * r(2)**2 - &
                (z33 / rel_p + tau_y * om_y**2 * dcs_y / 4) * r(3)**2 - &
                (2 * z34 / rel_p + tau_y * om_y**2 * s_y * ds_y / rel_p) * r(3) * r(4) - &
                (2 * z44 / rel_p + dcs_y / (4 * rel_p2)) * r(4)**2

  mat6_i(6,6) = 1

  mat6 = matmul(mat6_i, mat6)

endif

! Ending coords

orbit%vec(1) = c_x * r(1) + s_x * r(2) / rel_p + x_c
orbit%vec(2) = tau_x * om_x**2 * rel_p * s_x * r(1) + c_x * r(2)
orbit%vec(3) = c_y * r(3) + s_y * r(4) / rel_p
orbit%vec(4) = tau_y * om_y**2 * rel_p * s_y * r(3) + c_y * r(4)
orbit%vec(5) = r(5) + z0 + z1 * r(1) + z2 * r(2) + &
               z11 * r(1)**2 + z12 * r(1) * r(2) + z22 * r(2)**2 + &
               z33 * r(3)**2 + z34 * r(3) * r(4) + z44 * r(4)**2 

end subroutine sbend_body_with_k1_map
