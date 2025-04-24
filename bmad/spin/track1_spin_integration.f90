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
!print *, "spin:", temp_end%spin

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

real(rp) entrance_orb(6), body_orb(6), exit_orb(6), quat(4), quat_entrance(4), quat_exit(4), omega_vec(3)
real(rp) m, e, gma, anom, q, l, chi, xi
real(rp) ks, c, s, a, b, cc, zeta, sc, omega, ch, sh, e1, e2, k0, kx, g, k1, t, omegax, sx, cx, shx, chx, theta
logical spin_fringe
integer key

! Constants

m = mass_of(ele%ref_species)
e = ele%value(e_tot$)
gma = e/m
anom = anomalous_moment_of(ele%ref_species)
q = charge_of(ele%ref_species)
l = ele%value(l$)

chi = 1.0_rp + anom*gma
xi = anom*(gma - 1.0_rp)

spin_fringe = is_true(ele%value(spin_fringe_on$))

entrance_orb = no_fringe_start%vec
body_orb = start_orb%vec
exit_orb = end_orb%vec

! Spin tracking

quat_entrance(1) = 1.0_rp
quat_entrance(2) = 0.0_rp
quat_entrance(3) = 0.0_rp
quat_entrance(4) = 0.0_rp

quat(1)= 1.0_rp
quat(2) = 0.0_rp
quat(3) = 0.0_rp
quat(4) = 0.0_rp

quat_exit(1) = 1.0_rp
quat_exit(2) = 0.0_rp
quat_exit(3) = 0.0_rp
quat_exit(4) = 0.0_rp

key = ele%key
if (key == quadrupole$ .and. abs(ele%value(k1$)) < 1d-20) key = pipe$

select case (key)

case (pipe$)
  quat = quat
  
! Solenoid

case (solenoid$)
  ks = ele%value(ks$)

  c = cos(ks*l)
  s = sin(ks*l)

  if (spin_fringe) then
    a = -0.25_rp*chi*ks*entrance_orb(1)
    b = -0.25_rp*chi*ks*entrance_orb(3)
    zeta = sqrt(a**2 + b**2)
    sc = sinc(zeta)
   
    quat_entrance(1) = -cos(zeta)
    quat_entrance(2) = a*sc
    quat_entrance(3) = b*sc
    
    a = 0.25_rp*chi*ks*exit_orb(1)
    b = 0.25_rp*chi*ks*exit_orb(3)
    zeta = sqrt(a**2 + b**2)
    sc = sinc(zeta)
    
    quat_exit(1) = -cos(zeta)
    quat_exit(2) = a*sc
    quat_exit(3) = b*sc
  endif
  
  a = -0.5_rp*xi*(s*body_orb(2)+(1.0_rp-c)*body_orb(4))
  b = -0.5_rp*xi*((c-1.0_rp)*body_orb(2)+s*body_orb(4))
  cc = -0.5_rp*(1.0_rp+anom)*ks*l*(body_orb(6)-1.0_rp)
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


  if (spin_fringe) then
    t = tan(e1)
    
    a = -0.5_rp*chi*k0*t*entrance_orb(3)
    b = -0.5_rp*chi*k0*t*entrance_orb(1)
    zeta = sqrt(a**2 + b**2)
    sc = sinc(zeta)
    
    quat_entrance(1) = -cos(zeta)
    quat_entrance(2) = a*sc
    quat_entrance(3) = b*sc
    
    t = tan(e2)
    a = -0.5_rp*chi*k0*t*exit_orb(3)
    b = -0.5_rp*chi*k0*t*exit_orb(1)
    zeta = sqrt(a**2 + b**2)
    sc = sinc(zeta)
    
    quat_exit(1) = -cos(zeta)
    quat_exit(2) = a*sc
    quat_exit(3) = b*sc
    
  endif

  if (k1 == 0) then
    theta = k0*l
    c = cos(theta)
    s = sin(theta)
    a = -theta/(2.0_rp*gma)*(gma+(anom+gma)*body_orb(6))
    a = a + 0.5_rp*chi*(k0*s*body_orb(1)+(1.0_rp-c)*body_orb(2)-s*body_orb(6)+(1.0_rp+body_orb(6))*theta)
    b = -0.5_rp*xi*theta*body_orb(4)
    zeta = sqrt(a**2 + b**2)
    sc = sinc(zeta)
    
    quat(1) = -cos(zeta)
    quat(3) = a*sc
    quat(4) = b*sc
    
  else
    omega = sqrt(abs(k1))
    c = cos(omega*l)
    ch = cosh(omega*l)
    s = sin(omega*l)
    sh = sinh(omega*l)
    
    omegax = sqrt(abs(kx))
    cx = cos(omegax*l)
    chx = cosh(omegax*l)
    sx = sin(omegax*l)
    shx = sinh(omegax*l)
    
    if (k1 > 0 .and. kx > 0) then
      a = 0.5_rp*chi*(sh*omega*body_orb(3)+(ch-1.0_rp)*body_orb(4))
      b = -chi*omega**2/(2.0_rp*omegax**3)*(omegax*cx*body_orb(2)+sx*(k0*body_orb(6)-omegax**2*body_orb(1)))
      b = b + 0.5_rp*(chi*omega**2/(omegax**2)*(body_orb(2)+k0*l*body_orb(6))-k0*l/gma*(anom*body_orb(6)+gma*(1.0_rp+body_orb(6)-chi)))
      cc = -0.5_rp*k0*xi*((ch-1.0_rp)*body_orb(3)+sh/omega*body_orb(4))
    else if (k1 < 0 .and. kx > 0) then
      a = 0.5_rp*chi*(-s*omega*body_orb(3)+(c-1.0_rp)*body_orb(4))
      b = chi*omega**2/(2.0_rp*omegax**3)*(omegax*cx*body_orb(2)+sx*(k0*body_orb(6)-omegax**2*body_orb(1)))
      b = b + 0.5_rp*(-chi*omega**2/(omegax**2)*(body_orb(2)+k0*l*body_orb(6))-k0*l/gma*(anom*body_orb(6)+gma*(1.0_rp+body_orb(6)-chi)))
      cc = -0.5_rp*k0*xi*((c-1.0_rp)*body_orb(3)+s/omega*body_orb(4))
    else ! (k1 < 0 .and. kx > 0)
      a = 0.5_rp*chi*(-s*omega*body_orb(3)+(c-1.0_rp)*body_orb(4))
      b = -chi*omega**2/(2.0_rp*omegax**3)*(-omegax*chx*body_orb(2)-sx*(k0*body_orb(6)+omegax**2*body_orb(1)))
      b = b + 0.5_rp*(-chi*omega**2/(omegax**2)*(body_orb(2)+k0*l*body_orb(6))-k0*l/gma*(anom*body_orb(6)+gma*(1.0_rp+body_orb(6)-chi)))
      cc = -0.5_rp*k0*xi*((c-1.0_rp)*body_orb(3)+s/omega*body_orb(4))
    endif
    
    zeta = sqrt(a**2 + b**2 + cc**2)
    sc = sinc(zeta)
    
    quat(1) = -cos(zeta)
    quat(2) = a*sc
    quat(3) = b*sc
    quat(4) = cc*sc
    
  endif 
  
! Quadrupole

case (quadrupole$)
  k1 = ele%value(k1$)

  omega = sqrt(abs(k1))
  
  c = cos(omega*l)
  ch = cosh(omega*l)
  s = sin(omega*l)
  sh = sinh(omega*l)
  
  if (k1 > 0) then
    a = 0.5_rp*chi*(body_orb(4)*(ch-1.0_rp)+body_orb(3)*omega*sh)
    b = 0.5_rp*chi*(body_orb(2)*(1.0_rp-c)+body_orb(1)*omega*s)
  else
    a = 0.5_rp*chi*(body_orb(4)*(c-1.0_rp)-body_orb(3)*omega*s)
    b = 0.5_rp*chi*(body_orb(2)*(1.0_rp-ch)-body_orb(1)*omega*sh)
  endif
  
  zeta = sqrt(a**2 + b**2)
  sc = sinc(zeta)
  
  
  quat(1) = -cos(zeta)
  quat(2) = a*sc
  quat(3) = b*sc

  


case default
  print *, 'HELP! I SHOULD NOT BE HERE!!!!', key_name(key)
  return

end select

if (spin_fringe) then
  quat = quat_mul(quat_exit,quat,quat_entrance)
endif

omega_vec = quat_to_omega(quat/norm2(quat))

end function magnus_omega

end subroutine track1_spin_integration

