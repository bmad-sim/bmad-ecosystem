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
real(rp) m, e, gma, anom, q, l, chi, xi, rel_p, x0, px0, y0, py0, k2, alpha, taux, tauy, xc, omegay, pr, upsilon, sigma, a2, b2, cc2, nu
real(rp) ks, c, s, a, b, cc, zeta, sc, omega, k0, kx, g, k1, t, omegax, sx, cx, sy, cy, sh, ch, cx2, sx2, pry, eta, mu, iota, pl, xd
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
iota = 2.0_rp*gma - 1.0_rp
nu = 2.0_rp*(1.0_rp + anom) - chi

entrance_orb = no_fringe_start%vec
body_orb = start_orb%vec
exit_orb = end_orb%vec

x0 = body_orb(1)
px0 = body_orb(2)
y0 = body_orb(3)
py0 = body_orb(4)

! Spin tracking

quat(1) = 1.0_rp
quat(2) = 0.0_rp
quat(3) = 0.0_rp
quat(4) = 0.0_rp

key = ele%key
if (key == solenoid$ .and. abs(ele%value(ks$)) < 1d-20) key = pipe$
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
  b = b + pr*(2.0_rp*(c-1.0_rp)*px0+2.0_rp*s*py0+c*ks*y0-ks*(s*x0+y0))
  b = -b * xi/(8.0_rp*rel_p**2)
  
  cc = 4.0_rp*pr*(px0**2+py0**2)*s*(1.0_rp+anom*iota)+ks**3*L*(x0**2+y0**2)*(1.0_rp+anom*iota)- &
  ks**2*pr*s*(x0**2+y0**2)*(1.0_rp+anom*iota)+4.0_rp*ks*((-1.0_rp+c)*pr*(px0*x0+py0*y0)* &
  (1.0_rp+anom*iota)+L*(px0**2+py0**2+4.0_rp*rel_p**2+anom*(4.0_rp*rel_p**2+(px0**2+py0**2)*iota)))
  cc = cc/(32.0_rp*rel_p**3)
  
  cc2 = 2.0_rp*ks**2*L*(x0*(-2.0_rp*py0+2.0_rp*px0*s+ks*x0)+2.0_rp*(px0+py0*s)*y0+ks*y0**2)+pr* &
  (4.0_rp*s*(px0**2+py0**2)+4.0_rp*ks*py0*(s*x0-2.0_rp*y0)-4.0_rp*ks*px0*(2.0_rp*x0+s*y0)- &
  3.0_rp*ks**2*s*(x0**2+y0**2))+c*ks*(8.0_rp*pr*(px0*x0+py0*y0)+L*(-4.0_rp*(px0**2+py0**2)+ks**2*(x0**2+y0**2)))
  cc2 = -cc2*pr*xi**2/(64.0_rp*rel_p**4)
   
  cc = cc + cc2
  
  zeta = sqrt(a**2 + b**2 + cc**2)
  sc = sinc(zeta)
  
  quat(1) = -cos(zeta)
  quat(2) = a*sc
  quat(3) = b*sc
  quat(4) = cc*sc
  
! SBend

case (sbend$)
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
  
  xd = x0-xc
  
  if (kx > 0) then
    cx = cos(omegax*L)
    sx = sin(omegax*L)/omegax
    cx2 = cos(omegax*L/2.0_rp)
    sx2 = sin(omegax*L/2.0_rp)
    taux = -1.0_rp
  elseif (kx < 0) then
    cx = cosh(omegax*L)
    sx = sinh(omegax*L)/omegax
    cx2 = cosh(omegax*L/2.0_rp)
    sx2 = sinh(omegax*L/2.0_rp)
    taux = 1.0_rp
  else
    cx = 1.0_rp
    sx = 0.0_rp
    taux = 0.0_rp
  endif

  if (k1 == 0) then
    upsilon = py0**2*xi-rel_p**2*chi
    pry = sqrt(rel_p**2-py0**2)
    if (g == 0) then
      a = -0.5_rp*k0*py0*xi*(L*px0-0.5_rp*k0*L**2)/(rel_p**2*pry)
      b = -0.5_rp*k0*(L*px0**2-k0*L**2*px0+k0**2*L**3/3.0_rp)*upsilon/(2.0_rp*rel_p**2*pry**3)-0.5_rp*k0*L*upsilon/(rel_p**2*pry)
      cc = -0.5_rp*k0*py0*xi*L/(rel_p**2)
    elseif (k0 == 0) then
      a = 0.0_rp
      b = -0.5_rp*g*L
      cc = 0.0_rp
    else
      sigma = rel_p*sx*(x0-3.0_rp*xc)*xd+L*(-(rel_p*xd*(x0-(2.0_rp+cx)*xc))+px0*(sx*xc+taux*px0/(rel_p*omegax**2)))-px0*(px0*sx+2.0_rp*(cx-1.0_rp)*rel_p*xc)*taux/(rel_p*omegax**2)
      
      a = -k0*py0*sx2*xi*(g*px0*sx*omegax+rel_p*(2+g*(x0*(1.0_rp+cx)+xc*(1.0_rp-cx)))*omegax)*(cx2*px0+rel_p*sx2*xd*taux*omegax)/(2.0_rp*rel_p**3*omegax**2*pry)
      
      b = -k0*L*upsilon/(2.0_rp*rel_p**2*pry)
      b = b + (-g*k0*upsilon*(sx*xd+L*xc+(cx-1.0_rp)*px0*taux/(rel_p*omegax**2))/(2.0_rp*rel_p**2*pry))
      b = b + (-k0*upsilon*(cx**2*px0*rel_p*xd+px0*rel_p*xd*taux*(-taux+sx**2*omegax**2)+L*(px0**2-rel_p**2*xd**2*taux*omegax**2)+cx*sx*(px0**2+rel_p**2*xd**2*taux*omegax**2))/(8.0_rp*rel_p**2*pry**3))
      b = b + (sigma*g*k0**2*py0**2*xi**2/(4.0_rp*rel_p**4*pry))
      b = b - 0.5_rp*g*L
      
      cc = -k0*py0*xi*(L+g*sx*xd+g*L*xc+(taux*(cx-1.0_rp)*g*px0/(rel_p*omegax**2)))/(2.0_rp*rel_p**2)
      cc = cc + (-sigma*g*k0**2*py0*xi*upsilon/(4.0_rp*rel_p**4*pry**2))
    endif
    
    zeta = sqrt(a**2 + b**2 + cc**2)
    sc = sinc(zeta)
    
    quat(1) = -cos(zeta)
    quat(2) = a*sc
    quat(3) = b*sc
    quat(4) = cc*sc
    
  else
    if (k1 > 0) then
      cy = cosh(omegay*L)
      sy = sinh(omegay*L)/omegay
    else
      cy = cos(omegay*L)
      sy = sin(omegay*L)/omegay
    endif
    
    eta = omegax**2 + omegay**2
    mu = omegax**2 - omegay**2
    c = cx*cy - 1.0_rp

    if (kx > 0 .and. k1 > 0) then
      a = k0*xi*omegay**2*(py0*omegax**2*(-cy*px0*sx-c*rel_p*xd*omegax)-(cx*px0*py0*sy+c*px0*rel_p*y0+px0*rel_p*sx*sy*y0*omegax**2-rel_p*xd*(py0*sx*sy+rel_p*(cy*sx-cx*sy)*y0)*omegax**3)*omegay**2)
      a = a + k1*chi*((cy-1.0_rp)*py0*rel_p*(1.0_rp+g*xc)*eta+g*py0*(px0*(cy*sx-cx*sy)+cx*cy*rel_p*xd+rel_p*xd*(sx*sy*omegax**2-1.0_rp))*omegay**2+rel_p*y0*omegay**2*(rel_p*sy*eta+g*px0*(sx*sy*omegay**2-c)+ &
      g*rel_p*(sy*xc*eta+cy*sx*xd*omegax**2+cx*sy*xd*omegay**2)))
      a = a * 0.5_rp/(rel_p**3*eta*omegay**2)
      
      b = 2.0_rp*g*px0*(k1*L*px0-2.0_rp*(cx-1.0_rp)*k0*rel_p+k1*(rel_p*(x0+3.0_rp*xc)-cx*(px0*sx+4.0_rp*rel_p*xc)-cx**2*rel_p*xd))*chi
      b = b + k0*py0*(L*py0+cy*py0*sy+(cy**2-1.0_rp)*y0)*nu*omegax**2-2.0_rp*g*rel_p*(-sx*xd*(2.0_rp*k0*rel_p+k1*px0*sx+k1*rel_p*(4.0_rp*xc+cx*xd))*chi+L*rel_p*(2.0_rp*rel_p-(2.0_rp*k0*xc+k1* &
      (x0**2-2.0_rp*x0*xc+3.0_rp*xc**2))*chi))*omegax**2
      b = b + 4.0_rp*k1*rel_p*chi*(px0*(1.0_rp-cx)+rel_p*(L*xc+sx*xd)*omegax**2)+k0*omegax**2*(L*(4.0_rp+px0**2)*chi+chi*(8.0_rp*L*(rel_p-1.0_rp)+4.0_rp*L*(rel_p-1.0_rp)**2+cx*px0**2*sx+ &
      (cx**2-1.0_rp)*px0*rel_p*xd*omegax-px0*rel_p*sx**2*xd*omegax**3+rel_p**3*(L-cx*sx)*xd**2*omegax**4)-rel_p**2*(L-cy*sy)*y0**2*nu*omegay**2+py0*y0*nu*((cy**2-1.0_rp)*(rel_p-1.0_rp)+rel_p*sy**2*omegay**2))
      b = b/(8.0_rp*rel_p**3*omegax**2)
      
      cc = k0*rel_p*(1.0_rp+g*xc)*(py0*sy+(cy-1.0_rp)*rel_p*y0)*eta+kx*py0*(-c*px0+cy*rel_p*sx*xd*omegax**2)
      cc = cc + kx*(px0*py0*sx*sy+px0*rel_p*(cy*sx-cx*sy)*y0+rel_p*xd*(cx*py0*sy+cx*cy*rel_p*y0+rel_p*y0*(sx*sy*omegax**2-1.0_rp)))*omegay**2
      cc = cc + k1*(px0*(py0*(c+sx*sy*omegax**2)+rel_p*y0*(cy*sx*omegax**2+cx*sy*omegay**2))+rel_p*(py0*(sy*xc*eta-cy*sx*xd*omegax**3+cx*sy*xd*omegax**3)+rel_p*y0*((cy-1.0_rp)*xc*eta+xd*omegax**3*(c-sx*sy*omegay**2))))
      cc = -0.5_rp*cc*xi/(rel_p**3*eta)
      
      a2 = rel_p*y0*omegay**2*((px0-rel_p*(L*xc-2.0_rp*sy*xc+sx*xd)*eta)*omegax**2+cx*px0*(eta-cy*omegax**2)+px0*(-1.0_rp+2.0_rp*sx*sy*omegax**2)*omegay**2+cx*(cy*px0+2.0_rp*rel_p*sy*xd*omegax**2)*omegay**2+cy* &
      (-px0*eta+rel_p*omegax**2*(-L*xc*eta+sx*xd*mu)))
      a2 = a2 + py0*(-(rel_p*xc*omegax**2*(eta-2.0_rp*cy*eta+omegax**2+L*sy*eta*omegay**2))+omegay**2*(-(rel_p*omegax**2*(x0+xd-2.0_rp*cx*cy*xd-sx*sy*xd*mu))+px0*(2.0_rp*cy*sx*omegax**2-sy*(eta+cx*mu))))
      a2 = 0.25_rp*a2*k0*kx*xi*chi/(rel_p**4*eta*omegax**2*omegay**2)
      
      b2 = k0*k1*(L-sy)*xi*chi*(py0**2-rel_p**2*y0**2*omegay**2)/(4.0_rp*rel_p**4*omegay**2)
      
      cc2 = px0*(py0*((cx-cy-1.0_rp)*eta+cx*cy*omegay**2+omegax**2*(2.0_rp-cx*cy+2.0_rp*sx*sy*omegay**2))+rel_p*y0*omegay**2*(2.0_rp*cy*sx*omegax**2-sy*(eta+cx*mu)))
      cc2 = cc2 - rel_p*omegax**2*(L*xc*eta*((cy+1.0_rp)*py0+rel_p*sy*y0*omegay**2)+py0*(-2.0_rp*sy*(xc*eta+cx*xd*omegay**2)+sx*xd*(eta-cy*mu))+rel_p*y0*(2.0_rp*xc*omegax**2-2.0_rp*cy*(xc*eta+cx*xd*omegay**2)+ &
      omegay**2*(2.0_rp*x0-sx*sy*xd*mu)))
      cc2 = 0.25_rp*cc2*k1*kx*chi**2/(rel_p**4*eta*omegax**2*omegay**2)
    elseif (kx > 0 .and. k1 < 0) then
      if (abs(mu) > 1d-20) then
        a = k0*xi*omegay**2*(py0*omegax**2*(-cy*px0*sx-c*rel_p*xd*omegax)+(cx*px0*py0*sy+c*px0*rel_p*y0+px0*rel_p*sx*sy*y0*omegax**2-rel_p*xd*(py0*sx*sy+rel_p*(cy*sx-cx*sy)*y0)*omegax**3)*omegay**2)
        a = a + k1*chi*(-((-1.0_rp+cy)*py0*rel_p*(1.0_rp+g*xc)*omegax**2)+(rel_p*y0*(-c*g*px0+rel_p*(sy+g*sy*xc+cy*g*sx*xd)*omegax**2)+py0*(-rel_p-cx*g*px0*sy+cy*(rel_p+g*(rel_p-1.0_rp)*(xc+cx*xd)+g*(px0*sx+xc+cx*xd))- &
        g*rel_p*(xc+xd-sx*sy*xd*omegax**2)))*omegay**2-rel_p*sy*(rel_p+g*(rel_p-1.0_rp)*(xc+cx*xd)+g*(px0*sx+xc+cx*xd))*y0*omegay**4)
        a = 0.5_rp*a/(rel_p**3*mu*omegay**2)
      
        cc = -k0*rel_p*(1.0_rp+g*xc)*(py0*sy+(-1.0_rp+cy)*rel_p*y0)*mu+kx*(c*px0*py0+px0*(py0*sx*sy+rel_p*(cy*sx-cx*sy)*y0)*omegay**2+rel_p*xd*(-cy*py0*sx*omegax**2+(cx*py0*sy+cx*cy*rel_p*y0+rel_p*y0*(-1.0_rp+sx*sy*omegax**2))*omegay**2))
        cc = cc - k1*(px0*(py0*(c+sx*sy*omegax**2)+rel_p*y0*(cy*sx*omegax**2-cx*sy*omegay**2))+rel_p*(py0*omegax**2*(-cy*sx*xd*omegax+sy*(xc+cx*xd*omegax))-py0*sy*xc*omegay**2+rel_p*y0*((-1.0_rp+cy)*xc*mu+xd*omegax**3*(-c+sx*sy*omegay**2))))
        cc = 0.5_rp*cc*xi/(rel_p**3*mu)
      
        a2 = -(rel_p*y0*omegay**2*((1.0_rp+cx)*(-1.0_rp+cy)*px0*omegax**2+px0*((-1.0_rp+cx)*(1.0_rp+cy)+2.0_rp*sx*sy*omegax**2)*omegay**2+rel_p*omegax**2*(((1.0_rp+cy)*L*xc-2.0_rp*sy*xc+sx*xd-cy*sx*xd)*omegax**2-((1.0_rp+cy)*L*xc+(1.0_rp+cy)* &
        sx*xd-2.0_rp*sy*(xc+cx*xd))*omegay**2)))
        a2 = a2 + py0*(-(rel_p*xc*omegax**2*mu*(-2.0_rp+2.0_rp*cy+L*sy*omegay**2))+omegay**2*(2.0_rp*cy*(px0*sx+cx*rel_p*xd)*omegax**2-px0*sy*((1.0_rp+cx)*omegax**2+(-1.0_rp+cx)*omegay**2)+rel_p*xd*omegax**2*(-2.0_rp+sx*sy*eta)))
        a2 = 0.25_rp*a2*k0*kx*xi*chi/(rel_p**4*omegax**2*omegay**2*mu)
      
        cc2 = px0*((1.0_rp+cx)*(-1.0_rp+cy)*py0*omegax**2+py0*((-1.0_rp+cx)*(1.0_rp+cy)+2.0_rp*sx*sy*omegax**2)*omegay**2+rel_p*y0*omegay**2*((2.0_rp*cy*sx-(1.0_rp+cx)*sy)*omegax**2-(-1.0_rp+cx)*sy*omegay**2))
        cc2 = cc2 - rel_p*omegax**2*(py0*(2.0_rp*sy*xc+(-1.0_rp+cy)*sx*xd)*omegax**2+py0*((1.0_rp+cy)*sx*xd-2.0_rp*sy*(xc+cx*xd))*omegay**2+L*xc*mu*(-((1.0_rp+cy)*py0)+rel_p*sy*y0*omegay**2)+rel_p*y0*(2.0_rp*(-1.0_rp+cy)*xc*omegax**2+ &
        (x0+xc-2.0_rp*cy*xc+xd-2.0_rp*cx*cy*xd-sx*sy*xd*omegax**2)*omegay**2-sx*sy*xd*omegay**4))
        cc2 = 0.25_rp*cc2*k1*kx*chi**2/(rel_p**4*omegax**2*omegay**2*mu)
      else
        a = k0*xi*omegax**2*(-cx**2*rel_p*(px0*y0+py0*xd*omegax)+rel_p*(px0*y0+py0*xd*omegax)*(1.0_rp+sx**2*omegax**2)+2.0_rp*cx*sx*(-px0*py0+rel_p**2*xd*y0*omegax**3)-2.0_rp*L*(px0*py0+rel_p**2*xd*y0*omegax**3))
        a = a + k1*chi*(rel_p*y0*(-((-1.0_rp+cx**2)*g*px0)+(g*px0*sx**2+2.0_rp*g*L*rel_p*xd+2.0_rp*rel_p*sx*(2.0_rp+2.0_rp*g*xc+cx*g*xd))*omegax**2)+py0*(4.0_rp+2.0_rp*g*L*px0+4.0_rp*(rel_p-1.0_rp)-2.0_rp*cx* &
        (2.0_rp*rel_p+g*px0*sx+2.0_rp*g*rel_p*xc)-cx**2*g*rel_p*xd+g*rel_p*(4.0_rp*xc+xd+sx**2*xd*omegax**2)))
        a = a/(8.0_rp*rel_p**3*omegax**2)
        
        cc = c*kx*px0*py0-(kx*py0*(px0*sx**2+2.0_rp*rel_p*(L+cx*sx)*xd)+kx*rel_p*(-2.0_rp*L*px0+2.0_rp*cx*px0*sx+c*rel_p*xd)*y0+4.0_rp*k0*rel_p*(1.0_rp+g*xc)*(py0*sx+(-1.0_rp+cx)*rel_p*y0))*omegax**2
        cc = cc + kx*rel_p**2*sx**2*xd*y0*omegax**4+k1*(c*px0*py0-px0*(py0*sx**2+2.0_rp*rel_p*(L+cx*sx)*y0)*omegax**2+rel_p*omegax**2*(-4.0_rp*py0*sx*xc+2.0_rp*py0*(L-cx*sx)*xd*omegax-rel_p*y0*(4.0_rp*(-1.0_rp+cx)*xc+xd*omegax*(c-sx**2*omegax**2))))
        cc = cc*xi/(8.0_rp*rel_p**3*omegax**2)
        
        a2 = L*px0*py0-px0*py0*sx-L*rel_p*(py0*sx*xc+rel_p*(-x0+(2.0_rp+cx)*xc)*y0)*omegax**2+rel_p*(-2.0_rp*(-1.0_rp+cx)*py0*xc-rel_p*sx*(x0-3.0_rp*xc)*y0*omegax**2)
        a2 = 0.25_rp*a2*k0*kx*xi*chi/(rel_p**4*omegax**2)
        
        cc2 = -py0*sx*(x0-3.0_rp*xc)+px0*sx*y0+2.0_rp*(-1.0_rp+cx)*rel_p*xc*y0+L*(py0*(x0-(2.0_rp+cx)*xc)-px0*y0+rel_p*sx*xc*y0*omegax**2)
        cc2 = -0.25_rp*cc2*k1*kx*chi**2/(rel_p**3*omegax**2)
      endif
      b = 2.0_rp*g*px0*(k1*L*px0-2.0_rp*(-1.0_rp+cx)*k0*rel_p+k1*(rel_p*(x0+3.0_rp*xc)-cx*(px0*sx+4.0_rp*rel_p*xc)-cx**2*rel_p*xd))*chi
      b = b + k0*py0*(L*py0+cy*py0*sy+(cy**2-1.0_rp)*y0)*nu*omegax**2-2.0_rp*g*rel_p*(-sx*xd*(2.0_rp*k0*rel_p+k1*px0*sx+k1*rel_p*(4.0_rp*xc+cx*xd))*chi+L*rel_p*(2.0_rp*rel_p-(2.0_rp*k0*xc+k1*(x0**2-2.0_rp*x0*xc+3.0_rp*xc**2))*chi))*omegax**2
      b = b + 4.0_rp*k1*rel_p*chi*(px0*(1.0_rp-cx)+rel_p*(L*xc+sx*xd)*omegax**2)+k0*omegax**2*(L*(4.0_rp+px0**2)*chi+chi*(8.0_rp*L*(rel_p-1.0_rp)+4.0_rp*L*(rel_p-1.0_rp)**2+cx*px0**2*sx+(-1.0_rp+cx**2)*px0*rel_p*xd*omegax- &
      px0*rel_p*sx**2*xd*omegax**3+rel_p**2*(L-cx*sx)*xd**2*omegax**4)+rel_p**2*(L-cy*sy)*y0**2*nu*omegay**2+py0*y0*nu*((-1.0_rp+cy**2)*(rel_p-1.0_rp)-rel_p*sy**2*omegay**2))
      b = b/(8.0_rp*rel_p**3*omegax**2)  
      
      b2 = -k0*k1*(L-sy)*xi*chi*(py0**2+rel_p**2*y0**2*omegay**2)/(4.0_rp*rel_p**4*omegay**2)
    elseif (kx < 0 .and. k1 > 0) then
      if (abs(mu) > 1d-20) then
        a = k0*xi*omegay**2*(py0*omegax**2*(-cy*px0*sx-c*rel_p*xd*omegax)+(cx*px0*py0*sy+c*px0*rel_p*y0-px0*rel_p*sx*sy*y0*omegax**2+rel_p*xd*(py0*sx*sy+rel_p*(cy*sx-cx*sy)*y0)*omegax**3)*omegay**2)
        a = a + k1*chi*((-1.0_rp+cy)*py0*rel_p*(1.0_rp+g*xc)*omegax**2+(rel_p*y0*(c*g*px0+rel_p*(sy+g*sy*xc+cy*g*sx*xd)*omegax**2)+py0*(rel_p-cy*(rel_p+g*px0*sx+g*rel_p*xc+cx*g*rel_p*xd)+g*(cx*px0*sy+rel_p*(xc+xd+sx*sy*xd*omegax**2))))*omegay**2-rel_p*sy*(rel_p+g*px0*sx+g*rel_p*xc+cx*g*rel_p*xd)*y0*omegay**4)
        a = 0.5_rp*a/(rel_p**3*mu*omegay**2)
      
        cc = kx*py0*(-px0*c-cy*rel_p*sx*xd*omegax**2)+kx*(px0*py0*sx*sy+px0*rel_p*(cy*sx-cx*sy)*y0+rel_p*xd*(cx*py0*sy+cx*cy*rel_p*y0-rel_p*y0*(1.0_rp+sx*sy*omegax**2)))*omegay**2
        cc = cc - k0*rel_p*(1.0_rp+g*xc)*(py0*sy+(-1.0_rp+cy)*rel_p*y0)*mu+k1*(px0*(py0*(c-sx*sy*omegax**2)-rel_p*y0*(cy*sx*omegax**2-cx*sy*omegay**2))+rel_p*(py0*omegax**2*(cy*sx*xd*omegax-sy*(xc+cx*xd*omegax))+py0*sy*xc*omegay**2+rel_p*y0*(-((-1.0_rp+cy)*xc*mu)+xd*omegax**3*(-c+sx*sy*omegay**2))))
        cc = 0.5_rp*cc*xi/(rel_p**3*mu)
      
        a2 = rel_p*y0*omegay**2*((1.0_rp+cx)*(-1.0_rp+cy)*px0*omegax**2+px0*((-1.0_rp+cx)*(1.0_rp+cy)-2.0_rp*sx*sy*omegax**2)*omegay**2+rel_p*omegax**2*((2.0_rp*sy*xc+(-1.0_rp+cy)*sx*xd)*omegax**2+((1.0_rp+cy)*sx*xd-2.0_rp*sy*(xc+cx*xd))*omegay**2-(1.0_rp+cy)*L*xc*mu))
        a2 = a2 + py0*(-((rel_p*xc*omegax**2*mu*(2.0_rp*(1.0_rp-cy)+L*sy*omegay**2))+omegay**2*(-2.0_rp*cy*(px0*sx+cx*rel_p*xd)*omegax**2+px0*sy*((1.0_rp+cx)*omegax**2+(-1.0_rp+cx)*omegay**2)+rel_p*xd*omegax**2*(2.0_rp+sx*sy*eta))))
        a2 = 0.25_rp*a2*k0*kx*xi*chi/(rel_p**4*omegax**2*omegay**2*mu)
      
        cc2 = px0*((1.0_rp+cx)*(-1.0_rp+cy)*py0*omegax**2+((-1.0_rp+cx)*(1.0_rp+cy)*py0+(-2.0_rp*py0*sx*sy+rel_p*(-2.0_rp*cy*sx+sy+cx*sy)*y0)*omegax**2)*omegay**2+(-1.0_rp+cx)*rel_p*sy*y0*omegay**4)
        cc2 = cc2 + rel_p*omegax**2*(py0*(2.0_rp*sy*xc+(-1.0_rp+cy)*sx*xd)*omegax**2+py0*((1.0_rp+cy)*sx*xd-2.0_rp*sy*(xc+cx*xd))*omegay**2-L*xc*mu*((1.0_rp+cy)*py0+rel_p*sy*y0*omegay**2)+ &
        rel_p*y0*(2.0_rp*(-1.0_rp+cy)*xc*omegax**2+(x0+xc-2.0_rp*cy*xc+xd-2.0_rp*cx*cy*xd+sx*sy*xd*omegax**2)*omegay**2+sx*sy*xd*omegay**4))
        cc2 = 0.25_rp*cc2*k1*kx*chi**2/(rel_p**4*omegax**2*omegay**2*mu)
      else
        a = k0*xi*omegax**2*(-cx**2*rel_p*(px0*y0+py0*xd*omegax)-rel_p*(px0*y0+py0*xd*omegax)*(-1.0_rp+sx**2*omegax**2)-2.0_rp*cx*sx*(px0*py0+rel_p**2*xd*y0*omegax**3)+L*(-2.0_rp*px0*py0+2.0_rp*rel_p**2*xd*y0*omegax**3))
        a = a + k1*chi*(rel_p*y0*(c*g*px0+(g*px0*sx**2+2.0_rp*g*L*rel_p*xd+2.0_rp*rel_p*sx*(2.0_rp+2.0_rp*g*xc+cx*g*xd))*omegax**2)+py0*(-2.0_rp*g*L*px0-4.0_rp*rel_p+2.0_rp*cx*(2.0_rp*rel_p+g*px0*sx+2.0_rp*g*rel_p*xc)+cx**2*g*rel_p*xd-g*rel_p*(x0+3.0_rp*xc-sx**2*xd*omegax**2)))
        a = a/(8.0_rp*rel_p**3*omegax**2)
        
        cc = c*kx*px0*py0+(kx*py0*(px0*sx**2+2.0_rp*rel_p*(L+cx*sx)*xd)+kx*rel_p*(-2.0_rp*L*px0+2.0_rp*cx*px0*sx+c*rel_p*xd)*y0+4.0_rp*k0*rel_p*(1.0_rp+g*xc)*(py0*sx+(-1.0_rp+cx)*rel_p*y0))*omegax**2
        cc = cc + kx*rel_p**2*sx**2*xd*y0*omegax**4+k1*(c*px0*py0+px0*(py0*sx**2+2.0_rp*rel_p*(L+cx*sx)*y0)*omegax**2+rel_p*omegax**2*(2.0_rp*py0*(2.0_rp*sx*xc-L*xd*omegax+cx*sx*xd*omegax)+rel_p*y0*(4.0_rp*(-1.0_rp+cx)*xc+xd*omegax*(c+sx**2*omegax**2))))
        cc = -cc*xi/(8.0_rp*rel_p**3*omegax**2)
        
        a2 = L*px0*py0-px0*py0*sx+L*rel_p*(py0*sx*xc+rel_p*(-x0+(2.0_rp+cx)*xc)*y0)*omegax**2+rel_p*(-2.0_rp*(-1.0_rp+cx)*py0*xc+rel_p*sx*(x0-3.0_rp*xc)*y0*omegax**2)
        a2 = -a2*k0*kx*xi*chi/(4.0_rp*rel_p**4*omegax**2)
        
        cc2 = -py0*sx*(x0-3.0_rp*xc)+L*py0*(x0-(2.0_rp+cx)*xc)+px0*sx*y0+2.0_rp*(-1.0_rp+cx)*rel_p*xc*y0-L*y0*(px0+rel_p*sx*xc*omegax**2)
        cc2 = 0.25_rp*cc2*k1*kx*chi**2/(rel_p**3*omegax**2)
      endif
      b = 2.0_rp*g*px0*(2.0_rp*(-1.0_rp+cx)*k0*rel_p-k1*(L*px0-cx*px0*sx+rel_p*x0)+(-3.0_rp+4.0_rp*cx)*k1*rel_p*xc+cx**2*k1*rel_p*xd)*chi
      b = b + k0*py0*(L*py0+cy*py0*sy+(cy**2-1.0_rp)*y0)*nu*omegax**2-2.0_rp*g*rel_p*(-sx*xd*(2.0_rp*k0*rel_p+k1*px0*sx+k1*rel_p*(4.0_rp*xc+cx*xd))*chi+L*rel_p*(2.0_rp*rel_p-(2.0_rp*k0*xc+k1*(x0**2-2.0_rp*x0*xc+3.0_rp*xc**2))*chi))*omegax**2
      b = b + 4.0_rp*k1*rel_p*chi*((-1.0_rp+cx)*px0+rel_p*(L*xc+sx*xd)*omegax**2)+k0*omegax**2*(L*(4.0_rp+px0**2)*chi+chi*(8.0_rp*L*(rel_p-1.0_rp)+4.0_rp*L*(rel_p-1.0_rp)**2+cx*px0**2*sx+(-1.0_rp+cx**2)*px0*rel_p*xd*omegax+px0*rel_p*sx**2*xd &
      *omegax**3-rel_p**2*(L-cx*sx)*xd**2*omegax**4)-rel_p**2*(L-cy*sy)*y0**2*nu*omegay**2+py0*y0*nu*((-1.0_rp+cy**2)*(rel_p-1.0_rp)+rel_p*sy**2*omegay**2))
      b = b/(8.0_rp*rel_p**3*omegax**2)
      
      b2 = k0*k1*(L-sy)*xi*chi*(py0**2-rel_p**2*y0**2*omegay**2)/(4.0_rp*rel_p**4*omegay**2)     
    else
      a = k0*xi*omegay**2*(py0*omegax**2*(-cy*px0*sx-c*rel_p*xd*omegax)-(cx*px0*py0*sy+c*px0*rel_p*y0-px0*rel_p*sx*sy*y0*omegax**2+rel_p*xd*(py0*sx*sy+rel_p*(cy*sx-cx*sy)*y0)*omegax**3)*omegay**2)
      a = a + k1*chi*(-((-1.0_rp+cy)*py0*rel_p*(1.0_rp+g*xc)*eta)+g*py0*(cx*px0*sy-cy*(px0*sx+cx*rel_p*xd)+rel_p*xd*(1.0_rp+sx*sy*omegax**2))*omegay**2+rel_p*y0*omegay**2*(rel_p*sy*eta+g*px0*(c+sx*sy*omegay**2)+g*rel_p* &
      (sy*xc*eta+cy*sx*xd*omegax**2+cx*sy*xd*omegay**2)))
      a = 0.5_rp*a/(rel_p**3*eta*omegay**2)
      
      b = 2.0_rp*g*px0*(2.0_rp*(-1.0_rp+cx)*k0*rel_p-k1*(L*px0-cx*px0*sx+x0+(rel_p-1.0_rp)*x0)+(-3.0_rp+4.0_rp*cx)*k1*rel_p*xc+cx**2*k1*rel_p*xd)*chi
      b = b + k0*py0*(L*py0+cy*py0*sy+(-1.0_rp+cy**2)*y0)*nu*omegax**2-2.0_rp*g*rel_p*(-sx*xd*(2.0_rp*k0*rel_p+k1*px0*sx+k1*rel_p*(4.0_rp*xc+cx*xd))*chi+L*rel_p*(2.0_rp*rel_p-(2.0_rp*k0*xc+k1*(x0**2-2.0_rp*x0*xc*3.0_rp*xc**2))*chi))* &
      omegax**2+4.0_rp*k1*rel_p*chi*((-1.0_rp+cx)*px0+rel_p*(L*xc+sx*xd)*omegax**2)
      b = b + k0*omegax**2*(L*(4.0_rp+px0**2)*chi+chi*(8.0_rp*L*(rel_p-1.0_rp)+4.0_rp*L*(rel_p-1.0_rp)**2+cx*px0**2*sx+(-1.0_rp+cx**2)*px0*rel_p*xd*omegax+px0*rel_p*sx**2*xd*omegax**3-rel_p**2*(L-cx*sx)*xd**2*omegax**4)+rel_p**2*(L-cy*sy)* &
      y0**2*nu*omegay**2+py0*y0*nu*((-1.0_rp+cy**2)*(rel_p-1.0_rp)-rel_p*sy**2*omegay**2))
      b = b/(8.0_rp*rel_p**3*omegax**2)
      
      cc = -k0*rel_p*(1.0_rp+g*xc)*(py0*sy+(-1.0_rp+cy)*rel_p*y0)*eta-kx*(c*px0*py0+px0*(py0*sx*sy+rel_p*(cy*sx-cx*sy)*y0)*omegay**2+rel_p*xd*(cy*py0*sx*omegax**2+(cx*py0*sy+cx*cy*rel_p*y0-rel_p*y0*(1.0_rp+sx*sy*omegax**2))*omegay**2))
      cc = cc + k1*(px0*(py0*(c-sx*sy*omegax**2)-rel_p*y0*(cy*sx*omegax**2+cx*sy*omegay**2))+rel_p*(-py0*sy*xc*eta+py0*(cy*sx-cx*sy)*xd*omegax**3-rel_p*y0*((-1.0_rp+cy)*xc*eta+xd*omegax**3*(c+sx*sy*omegay**2))))
      cc = 0.5_rp*cc*xi/(rel_p**3*eta)
      
      a2 = rel_p*y0*omegay**2*(px0*((1.0_rp+cx-cy)*eta-cx*cy*omegax**2+(-2.0_rp+cx*cy-2.0_rp*sx*sy*omegax**2)*omegay**2)+rel_p*omegax**2*((1.0_rp+cy)*L*xc*eta+sx*xd*(eta-cy*mu)-2.0_rp*sy*(xc*eta+cx*xd*omegay**2)))
      a2 = a2 + py0*(rel_p*xc*eta*omegax**2*(-2.0_rp+2.0_rp*cy+L*sy*omegay**2)+omegay**2*(2.0_rp*cy*px0*sx*omegax**2-px0*sy*(eta+cx*mu)+rel_p*xd*omegax**2*(2.0_rp*c-sx*sy*mu)))
      a2 = -0.25_rp*a2*k0*kx*xi*chi/(rel_p**4*eta*omegax**2*omegay**2)
      
      b2 = -k0*k1*(L-sy)*xi*chi*(py0**2+rel_p**2*y0**2*omegay**2)/(4.0_rp*rel_p**4*omegay**2)
      
      cc2 = px0*(py0*((1.0_rp+cx-cy)*eta-cx*cy*omegax**2+(-2.0_rp+cx*cy-2.0_rp*sx*sy*omegax**2)*omegay**2)+rel_p*y0*omegay**2*(-2.0_rp*cy*sx*omegax**2+sy*(eta+cx*mu)))
      cc2 = cc2 + rel_p*omegax**2*(L*xc*eta*((1.0_rp+cy)*py0-rel_p*sy*y0*omegay**2)+py0*(sx*xd*(eta-cy*omegax**2+cy*omegay**2)-2.0_rp*sy*(xc*eta+cx*xd*omegay**2))-rel_p*y0*(2.0_rp*(-1.0_rp+cy)*xc*eta+xd*omegay**2*(2.0_rp*c-sx*sy*mu)))
      cc2 = 0.25_rp*cc2*k1*kx*chi**2/(rel_p**4*eta*omegax**2*omegay**2)
    endif
    
    a = a + a2
    b = b + b2
    cc = cc + cc2
    
    zeta = sqrt(a**2 + b**2 + cc**2)
    sc = sinc(zeta)
    
    quat(1) = -cos(zeta)
    quat(2) = a*sc
    quat(3) = b*sc
    quat(4) = cc*sc
    
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
   cx = c
   sx = s/omega
   cy = ch
   sy = sh/omega
   tauy = 1.0_rp
  else
   cy = c
   sy = s/omega
   cx = ch
   sx = sh/omega
   tauy = -1.0_rp
  endif
 
  a = 0.5_rp*chi*k1*(sy*y0+tauy*(cy-1.0_rp)/omega**2*py0/rel_p)
  b = 0.5_rp*chi*k1*(sx*x0+tauy*(1.0_rp-cx)/omega**2*px0/rel_p)
  cc = px0*sx*(py0*sy+cy*rel_p*y0)+rel_p*x0*(cx*py0*sy+rel_p*(cx*cy-1.0_rp)*y0)
  cc = cc * (-0.5_rp*xi*k1/rel_p**2)
  cc = cc + 0.25_rp*chi**2/rel_p**2*((cx-cy)*px0*py0 + tauy*omega**2* &
  (px0*py0*sx*sy+px0*rel_p*(cy*sx-sy)*y0+rel_p*x0*(py0*(cx*sy-sx)+(cx*cy-1.0_rp)*rel_p*y0)))
  
  zeta = sqrt(a**2 + b**2 + cc**2)
  sc = sinc(zeta)
  
  quat(1) = -cos(zeta)
  quat(2) = a*sc
  quat(3) = b*sc
  quat(4) = cc*sc

! Sextupole

case (sextupole$)
  k2 = ele%value(k2$)
  
  ! First kick of a kick-drift-kick split
  px0 = px0 + 0.25_rp*k2*L*(y0**2-x0**2)
  py0 = py0 + 0.5_rp*k2*L*x0*y0
  
  pl = sqrt(rel_p**2-px0**2-py0**2)
  
  a = px0*xi*(L**2*(-3.0_rp*px0**2*py0+py0**3)-3.0_rp*L*pl*(2.0_rp*px0*py0*x0+y0*(px0**2-py0**2))- &
  3.0_rp*pl**2*(2.0_rp*px0*x0*y0+py0*(x0**2-y0**2)))
  a = a + rel_p**2*chi*(L*py0*(2.0_rp*L*px0+3.0_rp*pl*x0)+3.0_rp*pl*(L*px0+2.0_rp*pl*x0)*y0)
  a = a * k2*L/(12.0_rp*pl**3*rel_p**2)
  
  b = py0*xi*(L**2*(-3.0_rp*px0**2*py0+py0**3)-3.0_rp*L*pl*(2.0_rp*px0*py0*x0+y0*(px0**2-py0**2))- &
  3.0_rp*pl**2*(2.0_rp*px0*x0*y0+py0*(x0**2-y0**2)))
  b = b + rel_p**2*chi*(L**2*(px0**2-py0**2)+3.0_rp*pl**2*(x0**2-y0**2)+3.0_rp*L*pl*(px0*x0-py0*y0))
  b = b * k2*L/(12.0_rp*pl**3*rel_p**2)
  
  cc = L**2*(-3.0_rp*px0**2*py0+py0**3)-3.0_rp*L*pl*(2.0_rp*px0*py0*x0+y0*(px0**2-py0**2))- &
  3.0_rp*pl**2*(2.0_rp*px0*x0*y0+py0*(x0**2-y0**2))
  cc = cc * k2*L*xi/(12.0_rp*pl**2*rel_p**2)
  
  zeta = sqrt(a**2 + b**2 + cc**2)
  sc = sinc(zeta)
  
  quat(1) = -cos(zeta)
  quat(2) = a*sc
  quat(3) = b*sc
  quat(4) = cc*sc


case default
  print *, 'HELP! I SHOULD NOT BE HERE!!!!', key_name(key)
  return

end select

omega_vec = quat_to_omega(quat)

end function magnus_omega

end subroutine track1_spin_integration

