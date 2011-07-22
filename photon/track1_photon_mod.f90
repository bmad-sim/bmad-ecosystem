module track1_photon_mod

use bmad_struct
use bmad_interface

! This is for passing info into the field_calc_routine

type, private :: crystal_param_struct
  real(rp) cap_gamma, dtheta_sin_2theta, b_eff
  complex(rp) f0, fh
end type

! This is for passing info to the photon_depth_in_crystal routine used by zbrent.

type (ele_struct), private, pointer, save :: ele_com
type (coord_struct), private, pointer, save :: end_orb_com
real(rp), private, pointer, save :: vec0_com(:)
real(rp), private, pointer, save :: k0_outside_norm_com(:)

private reflection_point, photon_depth_in_crystal, e_field_calc

contains

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine track1_multilayer_mirror (ele, param, end_orb)
!
! Routine to track reflection from a multilayer_mirror.
!
! Input:
!   ele    -- ele_struct: Element tracking through.
!   param  -- lat_param_struct: lattice parameters.
!   end_orb    -- Coord_struct: phase-space coords to be transformed
!
! Output:
!   end_orb    -- Coord_struct: final phase-space coords
!-

subroutine track1_multilayer_mirror (ele, param, end_orb)

implicit none

type (ele_struct), target:: ele
type (coord_struct), target:: end_orb
type (lat_param_struct) :: param

real(rp) k0_outside_norm(3), hit_point(3), wavelength, sin_g, cos_g, f, m_in(3,3)
real(rp) vec0(3), s_len
real(rp), pointer :: v(:)

complex(rp) one, eps_1, eps_2, kz1, kz2, c1, c2

logical curved_surface

!

wavelength = ele%value(ref_wavelength$) / (1 + end_orb%vec(6))
v => ele%value

! (px, py, sqrt(1-px^2+py^2)) coords are with respect to laboratory reference trajectory.
! Convert this vector to k0_outside_norm which are coords with respect to crystal surface.
! k0_outside_norm is normalized to 1.

sin_g = sin(ele%value(graze_angle$))
cos_g = cos(ele%value(graze_angle$))
f = sqrt (1 - end_orb%vec(2)**2 - end_orb%vec(4)**2)

k0_outside_norm(1) =  cos_g * end_orb%vec(2) + f * sin_g
k0_outside_norm(2) = end_orb%vec(4)
k0_outside_norm(3) = -sin_g * end_orb%vec(2) + f * cos_g

! m_in 

m_in(1, 1:3) = [ cos_g, 0.0_rp, sin_g]
m_in(2, 1:3) = [0.0_rp, 1.0_rp, 0.0_rp]
m_in(3, 1:3) = [-sin_g, 0.0_rp, cos_g]

! If there is curvature, compute the reflection point which is where 
! the photon intersects the surface.

curved_surface = has_curved_surface(ele)
if (curved_surface) then
  ! vec0 is the body element coordinates of the photon.
  ! Note: The photon position in element entrance coords is [vec(1), vec(3), 0].
  vec0 = matmul (m_in, [end_orb_com%vec(1), end_orb_com%vec(3), 0.0_rp])
  call reflection_point (ele, end_orb, k0_outside_norm, vec0, hit_point, s_len)
  call to_curved_body_coords (ele, hit_point, k0_outside_norm, set$)
endif

!

eps_1 = 1 + r_e * wavelength**2 / (pi * ele%value(v1_unitcell$)) * cmplx(v(f0_re1$), v(f0_im1$))
eps_2 = 1 + r_e * wavelength**2 / (pi * ele%value(v2_unitcell$)) * cmplx(v(f0_re2$), v(f0_im2$))

kz1 = twopi * sqrt(k0_outside_norm(3)**2 + (eps_1)**2) / wavelength
kz2 = twopi * sqrt(k0_outside_norm(3)**2 + (eps_2)**2) / wavelength

c1 = exp(I_imaginary * kz1 * ele%value(d1_thickness$) / 2)
c2 = exp(I_imaginary * kz1 * ele%value(d2_thickness$) / 2)

one = cmplx(1.0_rp, 0.0_rp)

call multilayer_track (eps_1, eps_2, end_orb%e_field_x, end_orb%phase_x) ! pi polarization
call multilayer_track (one, one, end_orb%e_field_y, end_orb%phase_y) ! sigma polarization

!! s_len not being used !!!!!!!!!!!!!!!!!!!

!-----------------------------------------------------------------------------------------------
contains

subroutine multilayer_track (e1_coef, e2_coef, e_field, e_phase)

real(rp) e_field, e_phase
complex(rp) e1_coef, e2_coef, r11, r22, t12, t21, denom, k1, k2
complex(rp) a, v, f_minus, f_plus, r_ratio

!

k1 = e2_coef * kz1
k2 = e1_coef * kz2
denom = k1 + k2

r11 = c1**2 * (kz1 - kz2) / denom
r22 = c2**2 * (kz2 - kz1) / denom

t12 = 2 * c1 * c2 * kz2 / denom
t21 = 2 * c1 * c2 * kz1 / denom

end subroutine multilayer_track 

end subroutine track1_multilayer_mirror

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine track1_crystal (ele, param, end_orb)
!
! Routine to track reflection from a crystal.
!
! Input:
!   ele    -- ele_struct: Element tracking through.
!   param  -- lat_param_struct: lattice parameters.
!   end_orb    -- Coord_struct: phase-space coords to be transformed
!
! Output:
!   end_orb    -- Coord_struct: final phase-space coords
!-

subroutine track1_crystal (ele, param, end_orb)

implicit none

type (ele_struct), target:: ele
type (coord_struct), target:: end_orb
type (lat_param_struct) :: param
type (crystal_param_struct) c_param

real(rp), target :: m_in(3, 3)
real(rp), target :: k0_outside_norm(3)
real(rp) f, sin_alpha, cos_alpha, sin_psi, cos_psi, wavelength
real(rp) cos_g, sin_g, cos_tc, sin_tc, tilt_cor_mat(3,3)
real(rp) h_norm(3), kh_outside_norm(3), e_tot, pc, p_factor
real(rp) m_out(3, 3), y_out(3), x_out(3), k_out(3)
real(rp) temp_vec(3), direction(3), s_len
real(rp) gamma_0, gamma_h, b_err, vec0(3)
real(rp) hit_point(3)

logical curved_surface

!

wavelength = ele%value(ref_wavelength$) / (1 + end_orb%vec(6))

! (px, py, sqrt(1-px^2+py^2)) coords are with respect to laboratory reference trajectory.
! Convert this vector to k0_outside_norm which are coords with respect to crystal surface.
! k0_outside_norm is normalized to 1.

sin_g = sin(ele%value(graze_angle_in$))
cos_g = cos(ele%value(graze_angle_in$))
f = sqrt (1 - end_orb%vec(2)**2 - end_orb%vec(4)**2)

k0_outside_norm(1) =  cos_g * end_orb%vec(2) + f * sin_g
k0_outside_norm(2) = end_orb%vec(4)
k0_outside_norm(3) = -sin_g * end_orb%vec(2) + f * cos_g

! m_in 

m_in(1, 1:3) = [ cos_g, 0.0_rp, sin_g]
m_in(2, 1:3) = [0.0_rp, 1.0_rp, 0.0_rp]
m_in(3, 1:3) = [-sin_g, 0.0_rp, cos_g]

! If there is curvature, compute the reflection point which is where 
! the photon intersects the surface.

curved_surface = has_curved_surface(ele)
if (curved_surface) then
  ! vec0 is the body element coordinates of the photon.
  ! Note: The photon position in element entrance coords is [vec(1), vec(3), 0].
  vec0 = matmul (m_in, [end_orb_com%vec(1), end_orb_com%vec(3), 0.0_rp])
  call reflection_point (ele, end_orb, k0_outside_norm, vec0, hit_point, s_len)
  call to_curved_body_coords (ele, hit_point, k0_outside_norm, set$)
endif

! Construct h_norm = H * wavelength.

sin_alpha = sin(ele%value(alpha_angle$))
cos_alpha = cos(ele%value(alpha_angle$))
sin_psi = sin(ele%value(psi_angle$))
cos_psi = cos(ele%value(psi_angle$))
h_norm = [-cos_alpha, sin_alpha * sin_psi, sin_alpha * cos_psi] * wavelength / ele%value(d_spacing$)

! kh_outside_norm is the normalized outgoing wavevector outside the crystal

c_param%cap_gamma = r_e * wavelength**2 / (pi * ele%value(v_unitcell$)) 

kh_outside_norm = k0_outside_norm + h_norm

if (ele%value(b_param$) < 0) then ! Bragg
  kh_outside_norm(1) = -sqrt(1 - kh_outside_norm(2)**2 - kh_outside_norm(3)**2)
else
  kh_outside_norm(3) = sqrt(1 - kh_outside_norm(1)**2 - kh_outside_norm(2)**2)
endif

!-------------------------------------
! Calculate phase and intensity

if (ele%value(b_param$) < 0) then ! Bragg
  gamma_0 = k0_outside_norm(1)
  gamma_h = kh_outside_norm(1)
else
  gamma_0 = k0_outside_norm(3)
  gamma_h = kh_outside_norm(3)
endif

c_param%b_eff             = gamma_0 / gamma_h
c_param%dtheta_sin_2theta = dot_product(h_norm + 2 * k0_outside_norm, h_norm) / 2
c_param%f0                = cmplx(ele%value(f0_re$), ele%value(f0_im$)) 
c_param%fh                = cmplx(ele%value(fh_re$), ele%value(fh_im$))

p_factor = cos(2.0_rp*ele%value(graze_angle_in$))
call e_field_calc (c_param, ele, p_factor, end_orb%e_field_x, end_orb%phase_x)
call e_field_calc (c_param, ele, 1.0_rp,   end_orb%e_field_y, end_orb%phase_y)   ! Sigma polarization



!--------------------------------- 
! (x, px, y, py, x, pz) - Phase Space Calculations
! Translate to outgoing basis

! y_out = inverse(m_in) . m_tiltcorr . m_in . (0, 1, 0)

if (ele%value(tilt_corr$) == 0) then
  y_out = [0, 1, 0]
else
  sin_tc = sin(ele%value(tilt_corr$))
  cos_tc = cos(ele%value(tilt_corr$))
  tilt_cor_mat(1,1:3) = [cos_tc, -sin_tc, 0.0_rp]  
  tilt_cor_mat(2,1:3) = [sin_tc, cos_tc, 0.0_rp]
  tilt_cor_mat(3,1:3) = [0.0_rp, 0.0_rp, 1.0_rp]
  y_out = matmul(m_in, [0.0_rp, 1.0_rp, 0.0_rp])
  y_out = matmul(tilt_cor_mat, y_out)
  y_out = matmul(transpose(m_in), y_out)
endif

! x_out = vector orthogonal to y and z

x_out(1) =  y_out(2) * ele%value(kh_z_norm$) - y_out(3) * ele%value(kh_y_norm$)
x_out(2) = -y_out(1) * ele%value(kh_z_norm$) + y_out(3) * ele%value(kh_x_norm$)
x_out(3) =  y_out(1) * ele%value(kh_y_norm$) - y_out(2) * ele%value(kh_x_norm$)
  
k_out(1) = ele%value(kh_x_norm$)
k_out(2) = ele%value(kh_y_norm$)
k_out(3) = ele%value(kh_z_norm$)

m_out = reshape([x_out, y_out, k_out], [3, 3])

direction = matmul(transpose(m_out), kh_outside_norm)
end_orb%vec(2) = direction(1)
end_orb%vec(4) = direction(2)

! Compute position in phase space, backpropagating the ray

temp_vec = matmul(transpose(m_out), hit_point)
temp_vec = temp_vec - direction * temp_vec(3)

end_orb%vec(1) = temp_vec(1)
end_orb%vec(3) = temp_vec(2)
! %vec(5) doesn't include phase change due to wave nature of radiation
end_orb%vec(5) = temp_vec(3) 
if (curved_surface) end_orb%vec(5) = end_orb%vec(5) + s_len  ! IS this correct !!!!!!!!!!!!!1


end subroutine track1_crystal

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine e_field_calc (cp, ele, p_factor, e_field, e_phase)
!
! Routine to compute position where crystal reflects in crystal coordinates.
!
! Input:
!   cp       -- crystal_param_struct: Crystal parameters.
!   ele      -- ele_struct: Crystal element.
!   p_factor -- Real(rp)
!
! Output:
!   e_field -- Real(rp)
!   e_phase -- Real(rp)
!-

subroutine e_field_calc (cp, ele, p_factor, e_field, e_phase)

type (crystal_param_struct) cp
type (ele_struct) ele

real(rp) p_factor, e_field, e_phase

complex(rp) e_rel, e_rel_a, e_rel_b, eta, eta1, f_cmp, xi_0k_a, xi_hk_a, xi_0k_b, xi_hk_b

! Construct xi_0k = xi_0 / k and xi_hk = xi_h / k

eta = (-cp%b_eff * cp%dtheta_sin_2theta + cp%f0 * cp%cap_gamma * (1.0_rp - cp%b_eff)/2) / &
          (cp%cap_gamma * abs(p_factor) * sqrt(abs(cp%b_eff)) * cp%fh) 
eta1 = sqrt(eta**2 + sign(1.0_rp, cp%b_eff))
f_cmp = abs(p_factor) * sqrt(abs(cp%b_eff)) * cp%cap_gamma * cp%fh / 2

xi_0k_b = f_cmp * (eta - eta1)
xi_hk_b = f_cmp / (abs(cp%b_eff) * (eta - eta1))

xi_0k_a = f_cmp * (eta + eta1)
xi_hk_a = f_cmp / (abs(cp%b_eff) * (eta + eta1))

! Bragg

if (ele%value(b_param$) < 0) then 
  if (abs(eta+eta1) > abs(eta-eta1)) then
    e_rel = -2.0_rp * xi_0k_b / (p_factor * cp%cap_gamma * cp%fh)
  else
    e_rel = -2.0_rp * xi_0k_a / (p_factor * cp%cap_gamma * cp%fh)
  endif

  e_field = e_field * abs(e_rel)
  e_phase = atan2(aimag(e_rel), real(e_rel)) + e_phase

!---------------
! Laue calc

else 

  e_rel_a = -2.0_rp * xi_0k_a / (p_factor * cp%cap_gamma * cp%fh)
  e_rel_b = -2.0_rp * xi_0k_b / (p_factor * cp%cap_gamma * cp%fh)

  


endif

end subroutine e_field_calc

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine reflection_point (ele, end_orb, k0_outside_norm, vec0, hit_point, s_len)
!
! Private routine to compute position where crystal reflects in crystal coordinates.
!
! Output:
!   hit_point(3) -- Real(rp): point of reflection in crystal coordinates
!   s_len       -- Real(rp): distance for photon to propagate to hit_point
!-

subroutine reflection_point (ele, end_orb, k0_outside_norm, vec0, hit_point, s_len)

use nr, only: zbrent

implicit none

type (ele_struct), target :: ele
type (coord_struct), target :: end_orb

real(rp), target :: vec0(3), k0_outside_norm(3)
real(rp) :: s_len, hit_point(3)
real(rp) :: x1, x2, x_center, fa, fb

! init

ele_com => ele
end_orb_com => end_orb
vec0_com => vec0
k0_outside_norm_com => k0_outside_norm

! Assume flat crystal, compute s required to hit the intersection
! Choose a Bracket of 1m around this point.

if (ele_com%value(b_param$) < 0) then ! Bragg
  x_center = vec0_com(1) / k0_outside_norm_com(1)
else
  x_center = vec0_com(3) / k0_outside_norm_com(3)
endif

x1 = x_center
x2 = x_center
if (photon_depth_in_crystal(x_center) > 0) then
  do
    x1 = x_center - 0.1
    if (photon_depth_in_crystal(x1) < 0) exit
  enddo
else
  do
    x2 = x_center + 0.1
    if (photon_depth_in_crystal(x2) > 0) exit
  enddo
endif

s_len = zbrent (photon_depth_in_crystal, x1, x2, 1d-10)

! Compute the intersection point
hit_point = s_len * k0_outside_norm_com + vec0_com

end subroutine reflection_point

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Function photon_depth_in_crystal (s_len) result (delta_h)
! 
! Private routine to be used as an argument in zbrent. Propagates
! photon forward by a distance s_len. Returns delta_h = x-x0
! where x0 is the height of the crystal surface. 
! Since positive x points inward, positive delta_h => inside crystal.
!
! Input:
!   s_len   -- Real(rp): Place to position the photon.
!
! Output:
!   delta_h -- Real(rp): Depth of photon below surface in crystal coordinates.
!-

function photon_depth_in_crystal (s_len) result (delta_h)

implicit none

real(rp), intent(in) :: s_len
real(rp) :: delta_h
real(rp) :: c2, c3, c4
real(rp) :: point(3), r

c2 = ele_com%value(c2_curve_tot$)
c3 = ele_com%value(c3_curve_tot$)
c4 = ele_com%value(c4_curve_tot$)

point = s_len * k0_outside_norm_com + vec0_com

if (ele_com%value(b_param$) < 0) then ! Bragg
  r = point(3)
  delta_h = point(1) + c2*r*r + c3*r*r*r + c4*r*r*r*r
else
  r = point(1)
  delta_h = point(3) + c2*r*r + c3*r*r*r + c4*r*r*r*r
endif

end function photon_depth_in_crystal

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine to_curved_body_coords (ele, hit_point, vector, set)
!
! Routine to rotate between element body coords and effective body coords ("curved body coords") with 
! respect to the surface at the point of photon impact.
!
! Input:
!   ele         -- ele_struct: reflecting element
!   hit_point(3) -- real(rp): Point of photon impact.
!   vector(3)   -- real(rp): Vector to rotate.
!   set         -- Logical: True -> Transform body to curved body. False -> Transform curved body to body.
!
! Output:
!   vector(3)   -- Real(rp): Rotated vector.
!-

! Compute the slope of the crystal at that the point of impact.
! curverot transforms from standard body element coords to body element coords at point of impact.

subroutine to_curved_body_coords (ele, hit_point, vector, set)

implicit none

type (ele_struct) ele

real(rp) hit_point(3), vector(3)
real(rp), save :: c2_old = real_garbage$, c3_old = 0, c4_old = 0
real(rp), save :: curverot(3,3)
real(rp) c2, c3, c4, slope, curveangle, cos_c, sin_c

logical set

!

c2 = ele%value(c2_curve_tot$)
c3 = ele%value(c3_curve_tot$)
c4 = ele%value(c4_curve_tot$)

if (c2 /= c2_old .or. c3 /= c3_old .or. c4 /= c4_old) then
  slope = -2.0_rp * c2 * hit_point(3) - 3.0_rp * c3 * hit_point(3)**2 - 4.0_rp * c4 * hit_point(3)**4
  curveangle = atan(slope)
  cos_c = cos(curveangle)
  sin_c = sin(curveangle)
  curverot(1, 1:3) = [ cos_c, 0.0_rp, sin_c]
  curverot(2, 1:3) = [0.0_rp, 1.0_rp, 0.0_rp]
  curverot(3, 1:3) = [-sin_c, 0.0_rp, cos_c]
endif

vector = matmul(curverot, vector) ! Goto body element coords at point of photon impact

end subroutine to_curved_body_coords

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Function has_curved_surface (ele) result (has_curve)
!
! Routine to determine if the reflecting surface of an element is curved or not.
!
! Input:
!   ele -- ele_struct: Reflecting element.
!
! Output:
!   has_curve -- Logical: True if curved. False if not.
!-

function has_curved_surface (ele) result (has_curve)

implicit none

type (ele_struct) ele
logical has_curve

!

if (ele%value(c2_curve_tot$) /= 0 .or. ele%value(c3_curve_tot$) /= 0 .or. ele%value(c4_curve_tot$) /= 0) then
  has_curve = .true.
else
  has_curve = .false.
endif

end function has_curved_surface

end module
