!+
! subroutine track1_spin_integration (start_orb, ele, param, end_orb)
!
! Spin tracking through a single element by integrating the spin rotation vector omega.
!
! Note: spin tracking through a patch element is handled in track_a_patch since
! this is needed by runge_kutta tracking.
!
! Input :
!   start_orb  -- Coord_struct: Starting coords.
!   ele        -- Ele_struct: Element to track through.
!   param      -- lat_param_struct: Beam parameters.
!   end_orb    -- Coord_struct: Ending coords.
!
! Output:
!   end_orb    -- Coord_struct:
!     %spin(3)       -- Ending spin
!-

subroutine track1_spin_integration (start_orb, ele, param, end_orb)

use bmad_routine_interface, dummy => track1_spin_integration

implicit none

type (coord_struct) :: start_orb, end_orb
type (coord_struct) :: temp_start, temp_end, ave_orb, no_fringe_start
type (ele_struct) :: ele
type (lat_param_struct) :: param
type (fringe_field_info_struct) fringe_info, fringe_info2
type (em_field_struct) field

real(rp) spline_x(0:3), spline_y(0:3), omega(3), s_edge_track
real(rp) voltage, k_rf, phase
real(rp) fc, dxp, dyp, om(3), quat(0:3)

integer key, dir

logical magnus_tracking

character(*), parameter :: r_name = 'track1_spin_integration'

magnus_tracking = (ele%spin_tracking_method == magnus$)

! Spin tracking handled by track_a_patch for patch elements.

if (ele%key == patch$) return

! Beambeam
! This is used when *not* integrating.

if (ele%key == beambeam$) then
  fc = start_orb%p0c / ((1.0_rp + start_orb%beta**2) * charge_of(start_orb%species))
  dxp = (end_orb%vec(2) - start_orb%vec(2)) * fc
  dyp = (end_orb%vec(4) - start_orb%vec(4)) * fc
  field%E = [dxp,  dyp, 0.0_rp]
  field%B = [dyp, -dxp, 0.0_rp] / c_light
  om = spin_omega (field, end_orb, +1)
  quat = omega_to_quat(om)
  end_orb%spin = quat_rotate(quat, start_orb%spin)
  return
endif

! crab_cavity

if (ele%key == crab_cavity$) then
  ave_orb = start_orb
  ave_orb%vec  = (start_orb%vec  + end_orb%vec) / 2
  ave_orb%t    = (start_orb%t    + end_orb%t) / 2
  ave_orb%beta = (start_orb%beta + end_orb%beta) / 2

  voltage = e_accel_field(ele, voltage$) * rel_tracking_charge_to_mass(ave_orb, param%particle)
  k_rf = twopi * ele%value(rf_frequency$) / c_light
  phase = twopi * (ele%value(phi0$) + ele%value(phi0_multipass$) - &
          (particle_rf_time (ave_orb, ele, .false.) - rf_ref_time_offset(ele)) * ele%value(rf_frequency$))

  field%B(2) = -voltage * sin(phase) / c_light
  field%E(3) = voltage * k_rf * ave_orb%beta * ave_orb%vec(1) * cos(phase)

  omega = spin_omega(field, ave_orb, start_orb%direction * ele%orientation) * start_orb%time_dir
  call rotate_spin(omega, end_orb%spin)

  return
endif

! A slice_slave may or may not span a fringe. calc_next_fringe_edge will figure this out.

temp_start = start_orb
call calc_next_fringe_edge(ele, s_edge_track, fringe_info, temp_start, .true.)

call offset_particle(ele, set$, temp_start, set_hvkicks = magnus_tracking, set_spin = .true.)

no_fringe_start = temp_start

if (fringe_info%has_fringe .and. fringe_info%particle_at == first_track_edge$) then
  if (fringe_info%ds_edge /= 0) call track_a_drift(temp_start, ele%orientation*fringe_info%ds_edge)
  call apply_element_edge_kick(temp_start, fringe_info, ele, param, .true.)
  call calc_next_fringe_edge (ele, s_edge_track, fringe_info, end_orb, .false.)
endif

!

temp_end  = end_orb
temp_end%time_dir = -temp_end%time_dir
call calc_next_fringe_edge(ele, s_edge_track, fringe_info2, temp_end, .true.)
call offset_particle (ele, set$, temp_end, & ! drift_to_edge = temp_end%location, &
                                                  set_hvkicks = magnus_tracking, s_pos = temp_end%s-ele%s_start)

if (fringe_info2%has_fringe .and. fringe_info2%particle_at == first_track_edge$) then
  if (fringe_info2%ds_edge /= 0) call track_a_drift (temp_end, ele%orientation*fringe_info2%ds_edge)
  call apply_element_edge_kick (temp_end, fringe_info2, ele, param, .true.)
endif

temp_end%time_dir = -temp_end%time_dir

temp_end%spin = temp_start%spin

! temp_start%s = temp_end%s when, for example, in an RF cavity when ele is a super_slave covering
! a region outside the active region.

if (ele%value(l$) == 0 .or. ele%key == multipole$ .or. ele%key == ab_multipole$ .or. ele%key == taylor$) then
  temp_end%vec = (temp_end%vec + temp_start%vec) / 2
  call multipole_spin_tracking (ele, param, temp_end)
elseif (temp_start%s /= temp_end%s) then
  if (magnus_tracking) then
    omega = magnus_omega (ele, param, no_fringe_start, temp_start, temp_end)
  else
    call spline_fit_orbit (temp_start, temp_end, spline_x, spline_y)
    omega = trapzd_omega (ele, spline_x, spline_y, temp_start, temp_end, param)
    dir = start_orb%direction * ele%orientation * start_orb%time_dir
    if (ele%key == sbend$ .or. ele%key == rf_bend$) omega = omega + [0.0_rp, ele%value(g$)*ele%value(l$)*dir, 0.0_rp]
  endif
  call rotate_spin (omega, temp_end%spin)
endif
  

!----------

if (fringe_info%has_fringe .and. fringe_info%particle_at == second_track_edge$) then
  call apply_element_edge_kick (temp_end, fringe_info, ele, param, .true.)
endif

call offset_particle (ele, unset$, temp_end, set_hvkicks = magnus_tracking, set_spin = .true.)

end_orb%spin = temp_end%spin

!---------------------------------------------------------------------------------------------
contains

function trapzd_omega (ele, spline_x, spline_y, start_orb, end_orb, param) result (omega)

use super_recipes_mod, only: super_polint

implicit none

type q_array_struct
  real(rp) h
  real(rp) omega(3)
end type

integer, parameter ::  j_max = 10

type (q_array_struct) q_array(j_max), z(0:512)
type (ele_struct) ele
type (coord_struct) start_orb, end_orb, orb
type (lat_param_struct) param

real(rp) s0, s1, del_s, s, spline_x(0:3), spline_y(0:3), omega(3), deps
real(rp) dint, eps, quat(0:3)
real(rp), parameter :: eps_rel = 1d-5, eps_abs = 1d-8

integer j, k, n, n_pts

! Only integrate over where the field is finite.
! This will be the whole element except for RF cavities.

deps = 0.1_rp * sign_of(end_orb%s - start_orb%s) * bmad_com%significant_length

s0 = deps
s1 = end_orb%s - start_orb%s - deps

q_array(1)%h = 1
z(0)%omega = omega_func(s0, spline_x, spline_y, start_orb, end_orb, ele, param)
z(1)%omega = omega_func(s1, spline_x, spline_y, start_orb, end_orb, ele, param)

del_s = s1 - s0
q_array(1)%omega = quat_to_omega(quat_mul(omega_to_quat(z(1)%omega * abs(del_s) / 2), omega_to_quat(z(0)%omega * abs(del_s) / 2)))

do j = 2, j_max
  ! This is trapzd from NR
  n_pts = 2**(j-2)
  del_s = (s1 - s0) / (2 * n_pts)
  quat = omega_to_quat(z(0)%omega * abs(del_s) / 2)

  z(2:2*n_pts:2) = z(1:n_pts)

  do n = 1, n_pts
    s = s0 + del_s * (2*n - 1)
    z(2*n-1)%omega = omega_func(s, spline_x, spline_y, start_orb, end_orb, ele, param)
    quat = quat_mul(omega_to_quat(z(2*n-1)%omega * abs(del_s)), quat)
    if (n == n_pts) del_s = del_s / 2
    quat = quat_mul(omega_to_quat(z(2*n)%omega * abs(del_s)), quat)
  enddo
  
  q_array(j)%omega = quat_to_omega(quat)
  q_array(j)%h = q_array(j-1)%h / 4

  eps = eps_abs + eps_rel * sum(abs(q_array(j)%omega))

  select case (ele%key)
  case (wiggler$, rfcavity$, lcavity$);  if (j < 5) cycle  ! Cannot trust until have enough points
  end select

  do k = 1, 3
    call super_polint (q_array(1:j)%h, q_array(1:j)%omega(k), 0.0_rp, omega(k), dint)
    if (abs(dint) > eps .and. j < j_max) exit ! Failed test. Note: Last loop with j = j_max -> no test.
    if (k == 3) return                        ! Passed all tests or last loop
  enddo

enddo

end function trapzd_omega

!-----------------------------------------------------------------------------------
! contains

function omega_func (ds, spline_x, spline_y, start_orb, end_orb, ele, param) result (omega)

implicit none

type (coord_struct) start_orb, end_orb, orb
type (ele_struct) ele
type (em_field_struct) field
type (lat_param_struct) param

real(rp) spline_x(0:3), spline_y(0:3), omega(3), B(3)
real(rp) ds, dss, s_tot, s2, df

!

s_tot = end_orb%s - start_orb%s
orb = end_orb

orb%vec(5) = start_orb%vec(5) * (s_tot - ds) / s_tot + end_orb%vec(5) * ds / s_tot
orb%vec(6) = start_orb%vec(6) * (s_tot - ds) / s_tot + end_orb%vec(6) * ds / s_tot

df = orb%direction * (1 + orb%vec(6)) 
orb%vec(1) =       spline_x(0) + spline_x(1) * ds + spline_x(2) * ds**2 + spline_x(3) * ds**3
orb%vec(2) = df * (spline_x(1) + 2 * spline_x(2) * ds + 3 * spline_x(3) * ds**2)
orb%vec(3) =       spline_y(0) + spline_y(1) * ds + spline_y(2) * ds**2 + spline_y(3) * ds**3
orb%vec(4) = df * (spline_y(1) + 2 * spline_y(2) * ds + 3 * spline_y(3) * ds**2)

orb%t      = start_orb%t      * (s_tot - ds) / s_tot + end_orb%t      * ds / s_tot
orb%beta   = start_orb%beta   * (s_tot - ds) / s_tot + end_orb%beta   * ds / s_tot


if (ele%orientation == 1) then
  s2 = ds + start_orb%s - ele%s_start
else
  s2 = ele%s - (ds + start_orb%s)
endif

call em_field_calc (ele, param, s2, orb, .true., field)

! 1 + g*x term comes from the curved coordinates.

omega = spin_omega (field, orb, start_orb%direction * ele%orientation) * start_orb%time_dir
if (ele%key == sbend$ .or. ele%key == rf_bend$) omega = (1 + ele%value(g$) * orb%vec(1)) * omega

end function omega_func

!-----------------------------------------------------------------------------------
! contains

function magnus_omega (ele, param, no_fringe_start, start_orb, end_orb) result (omega_vec)

implicit none

type (coord_struct) no_fringe_start, start_orb, end_orb
type (ele_struct) ele
type (lat_param_struct) param

real(rp) entrance_orb(6), body_orb(6), exit_orb(6), quat(4), omega_vec(3)
real(rp) m, e, gma, anom, q, l, chi, xi, rel_p, x0, px0, y0, py0, k2, alpha, taux, tauy, xc, omegay, pr
real(rp) ks, c, s, a, b, cc, zeta, sc, omega, e1, e2, k0, kx, g, k1, t, omegax, sx, cx, sy, cy, sh, ch
integer key

! Constants

m = mass_of(start_orb%species)
rel_p = 1.0_rp + start_orb%vec(6)
e = start_orb%p0c * rel_p / start_orb%beta
gma = e/m
anom = anomalous_moment_of(start_orb%species)
q = charge_of(start_orb%species)
l = ele%value(l$)

chi = 1.0_rp + anom*gma
xi = anom*(gma - 1.0_rp)

entrance_orb = no_fringe_start%vec
body_orb = start_orb%vec
exit_orb = end_orb%vec

x0 = body_orb(1)
px0 = body_orb(2)
y0 = body_orb(3)
py0 = body_orb(4)

! Spin tracking

quat(1)= 1.0_rp
quat(2) = 0.0_rp
quat(3) = 0.0_rp
quat(4) = 0.0_rp

key = ele%key
if (key == sbend$ .and. abs(ele%value(g$)) < 1d-20 .and. abs(ele%value(dg$)) < 1d-20) key = quadrupole$
if (key == quadrupole$ .and. abs(ele%value(k1$)) < 1d-20) key = pipe$

select case (key)

case (pipe$, drift$, monitor$, instrument$, kicker$, hkicker$, vkicker$, rfcavity$)
  quat = quat
  
! Solenoid

case (solenoid$)
  ks = ele%value(ks$)
  
  pr = sqrt(rel_p**2-(px0+y0*ks/2.0_rp)**2-(py0-x0*ks/2.0_rp)**2)

  c = cos(ks*L/pr)
  s = sin(ks*L/pr)
  
  a = 2.0_rp*pr*(c-1.0_rp)*py0-2.0_rp*pr*s*px0+ks**2*L*y0
  a = a - ks*(2.0_rp*L*px0+(c-1.0_rp)*pr*x0+pr*s*y0)
  a = a * xi/(8.0_rp*rel_p**2)
  
  b = ks*L*(2.0_rp*py0+ks*x0)
  b = b + pr*(2.0_rp*(c-1.0_rp)*x0+2.0_rp*s*py0+c*ks*y0-ks*(s*x0+y0))
  b = -b * xi/(8.0_rp*rel_p**2)
  
  cc = 0.5_rp*(1.0_rp+anom)*ks/rel_p*L
  
  zeta = sqrt(a**2 + b**2 + cc**2)
  sc = sinc(zeta)
  
  quat(1) = -cos(zeta)
  quat(2) = a*sc
  quat(3) = b*sc
  quat(4) = cc*sc
  
  
! SBend

case (sbend$)
  e1 = ele%value(e1$)
  e2 = ele%value(e2$)
  k1 = ele%value(k1$)
  g = ele%value(g$)
  k0 = g + ele%value(dg$)
  kx = k1+g*k0
  omegay = sqrt(abs(k1)/rel_p)
  omegax = sqrt(abs(kx)/rel_p)
  
  if (abs(kx) > 0) then
    xc = (g*rel_p-k0)/kx
  else
    xc = 0.0_rp
  endif
  
  if (kx > 0) then
    cx = cos(omegax*L)
    sx = sin(omegax*L)/omegax
    taux = -1.0_rp
  elseif (kx < 0) then
    cx = cosh(omegax*L)
    sx = sinh(omegax*L)/omegax
    taux = 1.0_rp
  else
    cx = 1.0_rp
    sx = 0.0_rp
    taux = 0.0_rp
  endif

  if (k1 == 0) then
    if (g == 0) then
      a = 0.5_rp*chi*L*k0/rel_p
      b = -0.5_rp*xi*L*k0/rel_p*py0/rel_p
    else
      alpha = taux*(cx-1.0_rp)*g*px0
      alpha = alpha + rel_p*omegax**2*(L*(1.0_rp+g*xc)+g*sx*(x0-xc))
      alpha = alpha * k0/(2.0_rp*rel_p**2*omegax**2)
    
      a = chi*alpha - 0.5_rp*g*L
      b = -xi*py0*alpha/rel_p
    endif
    
    zeta = sqrt(a**2 + b**2)
    sc = sinc(zeta)
    
    quat(1) = -cos(zeta)
    quat(3) = a*sc
    quat(4) = b*sc
    
  else
    if (k1 > 0) then
      cy = cosh(omegay*L)
      sy = sinh(omegay*L)/omegay
      tauy = 1.0_rp
    else
      cy = cos(omegay*L)
      sy = sin(omegay*L)/omegay
      tauy = -1.0_rp
    endif
    
    a = 0.5_rp*chi*k1/rel_p*(sy*y0+tauy*(cy-1.0_rp)/omegay**2*py0/rel_p)
    b = 0.5_rp*chi*L*k0/rel_p-0.5_rp*g*L+0.5_rp*chi*kx/rel_p*(sx*(x0-xc)+L*xc+taux*(cx-1.0_rp)/omegax**2*px0/rel_p)
    c = 0.5_rp*(1.0_rp+anom)*k0/rel_p*((cy-1.0_rp)*y0+sy*py0/rel_p)
    
    zeta = sqrt(a**2 + b**2 + c**2)
    sc = sinc(zeta)
    
    quat(1) = -cos(zeta)
    quat(2) = a*sc
    quat(3) = b*sc
    quat(4) = c*sc
    
  endif 
  
! Quadrupole

case (quadrupole$)
  k1 = ele%value(k1$)/rel_p

  omega = sqrt(abs(k1))
  
  s = sin(omega*L)
  sh = sinh(omega*L)
  c = cos(omega*L)
  ch = cosh(omega*L)
 
  if (k1 > 0) then
   cx = (1.0_rp-c)/omega**2
   sx = s/omega
   cy = (-1.0_rp+ch)/omega**2
   sy = sh/omega
  else
   cy = (1.0_rp-c)/omega**2
   sy = s/omega
   cx = (-1.0_rp+ch)/omega**2
   sx = sh/omega
  endif
 
  a = 0.5_rp*chi*k1*(sy*y0+cy*py0/rel_p)
  b = 0.5_rp*chi*k1*(sx*x0+cx*px0/rel_p)
  
  zeta = sqrt(a**2 + b**2)
  sc = sinc(zeta)
  
  quat(1) = -cos(zeta)
  quat(2) = a*sc
  quat(3) = b*sc

! Sextupole

case (sextupole$)
  k2 = ele%value(k2$)
  
  a = L*py0*(2.0_rp*L*px0+3.0_rp*rel_p*x0)
  a = a + 3.0_rp*rel_p*y0*(L*px0+2.0_rp*rel_p*x0)
  a = a * k2*L*chi/(12.0_rp*rel_p**3)
  
  b = L**2*(px0**2-py0**2)
  b = b + 3.0_rp*rel_p**2*(x0**2-y0**2)
  b = b + 3.0_rp*L*rel_p*(px0*x0-py0*y0)
  b = b * k2*L*chi/(12.0_rp*rel_p**3)
  
  zeta = sqrt(a**2 + b**2)
  sc = sinc(zeta)
  
  quat(1) = -cos(zeta)
  quat(2) = a*sc
  quat(3) = b*sc


case default
  print *, 'HELP! I SHOULD NOT BE HERE!!!!', key_name(key)
  return

end select

omega_vec = quat_to_omega(quat/norm2(quat))

end function magnus_omega

end subroutine track1_spin_integration

