module track1_photon_mod

use bmad_struct
use bmad_interface
use track1_mod

! This is for passing info into the field_calc_routine

type, private :: crystal_param_struct
  real(rp) cap_gamma, dtheta_sin_2theta, b_eff
  complex(rp) f0, fh
end type

private reflection_point, e_field_calc

contains

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine track1_multilayer_mirror (ele, param, end_orb)
!
! Routine to track reflection from a multilayer_mirror.
! Basic equations are from Kohn, "On the Theory of Reflectivity of an X-Ray Multilayer Mirror".
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

real(rp) k0_outside_norm(3), hit_point(3), wavelength, sin_g, cos_g
real(rp) s_len, kz_air, temp_vec(3)
real(rp), pointer :: val(:)

complex(rp) zero, xi_1, xi_2, kz1, kz2, c1, c2

logical curved_surface

character(32), parameter :: r_name = 'track1_multilayer_mirror'

!

val => ele%value
wavelength = c_light * h_planck / end_orb%p0c

! (px, py, sqrt(1-px^2+py^2)) coords are with respect to laboratory reference trajectory.
! Convert this vector to k0_outside_norm which are coords with respect to crystal surface.
! k0_outside_norm is normalized to 1.

sin_g = sin(val(graze_angle$))
cos_g = cos(val(graze_angle$))

k0_outside_norm(1) =  cos_g * end_orb%vec(2) + sin_g * end_orb%vec(6)
k0_outside_norm(2) =          end_orb%vec(4)
k0_outside_norm(3) = -sin_g * end_orb%vec(2) + cos_g * end_orb%vec(6)

! If there is curvature, compute the reflection point which is where 
! the photon intersects the surface.

curved_surface = has_curved_surface(ele)
if (curved_surface) then
  call reflection_point (ele, cos_g, sin_g, end_orb, k0_outside_norm, hit_point, s_len)
  call to_curved_body_coords (ele, hit_point, k0_outside_norm, set$)
else
  hit_point = [end_orb%vec(1) * cos_g, end_orb%vec(3), -end_orb%vec(1) * sin_g]
endif

! Check aperture

if (ele%aperture_at == surface$) then
  end_orb%vec(1:5:2) = hit_point ! This is temporary
  call check_aperture_limit (end_orb, ele, surface$, param)
  if (end_orb%state /= alive$) return
endif

! Note: Koln z-axis = Bmad x-axis.
! Note: f0_re and f0_im are both positive.

xi_1 = cmplx(-val(f0_re1$), val(f0_im1$)) * r_e * wavelength**2 / (pi * val(v1_unitcell$)) 
xi_2 = cmplx(-val(f0_re2$), val(f0_im2$)) * r_e * wavelength**2 / (pi * val(v2_unitcell$)) 

kz1 = twopi * sqrt(k0_outside_norm(1)**2 + xi_1) / wavelength
kz2 = twopi * sqrt(k0_outside_norm(1)**2 + xi_2) / wavelength
kz_air = twopi * k0_outside_norm(1) / wavelength

c1 = exp(I_imaginary * kz1 * val(d1_thickness$) / 2)
c2 = exp(I_imaginary * kz2 * val(d2_thickness$) / 2)

zero = cmplx(0.0_rp, 0.0_rp)

call multilayer_track (xi_1, xi_2, end_orb%e_field_x, end_orb%phase_x)     ! pi polarization
call multilayer_track (zero, zero, end_orb%e_field_y, end_orb%phase_y)     ! sigma polarization

! Reflect momentum vector

k0_outside_norm(1) = -k0_outside_norm(1)

! Rotate back to uncurved element coords

if (curved_surface) then
  call to_curved_body_coords (ele, hit_point, k0_outside_norm, unset$)
endif

! Translate momentum to laboratory exit coords

end_orb%vec(2) = k0_outside_norm(1) * cos_g + k0_outside_norm(3) * sin_g
end_orb%vec(4) = k0_outside_norm(2)
end_orb%vec(6) = sqrt(1 - end_orb%vec(2)**2 - end_orb%vec(4)**2)

! Compute position, backpropagating the ray
!! end_orb%vec(5) not computed properly

temp_vec(1) = cos_g * hit_point(1) - sin_g * hit_point(3)
temp_vec(2) = hit_point(2)
temp_vec(3) = sin_g * hit_point(1) + cos_g * hit_point(3)

temp_vec = temp_vec - end_orb%vec(2:6:2) * temp_vec(3)

end_orb%vec(1) = temp_vec(1)
end_orb%vec(3) = temp_vec(2)
end_orb%vec(5) = temp_vec(3) 


!-----------------------------------------------------------------------------------------------
contains

subroutine multilayer_track (xi1_coef, xi2_coef, e_field, e_phase)

real(rp) e_field, e_phase

complex(rp) xi1_coef, xi2_coef, r_11, r_22, tt, denom, k1, k2
complex(rp) a, v, f_minus, f_plus, r_ratio, f, r, ttbar, nu, exp_half, exp_n
complex(rp) Rc_n_top, R_tot, k_a, r_aa, Rc_1_bot, Rc_1_top

integer i, n1

! Rc_n_top is the field ratio for the top layer of cell n.
! Rc_n_bot is the field ratio for the bottom layer of cell n.
! Rc_1_bot for the bottom layer just above the substrate is assumed to be zero.
! Upgrade: If we knew the substrate material we would not have to make this assumption.

Rc_1_bot = 0

! Compute Rc_1_top.
! The top layer of a cell is labeled "2" and the bottom "1". See Kohn Eq 6.

k1 = (1 + xi2_coef) * kz1
k2 = (1 + xi1_coef) * kz2
denom = k1 + k2

r_11 = c1**2 * (k1 - k2) / denom
r_22 = c2**2 * (k2 - k1) / denom

tt = 4 * k1 * k2 * (c1 * c2 / denom)**2    ! = t_12 * t_21

Rc_1_top = r_22 + tt * Rc_1_bot / (1 - r_11 * Rc_1_bot)

! Now compute the single cell factors. See Kohn Eq. 12.
! Note: If you go through the math you will find r = r_bar.

f = tt / (1 - r_11**2)
r = r_22 + r_11 * f
ttbar = f**2

! Calc Rc_n_top. See Kohn Eq. 21. Note that there are n-1 cells in between. 

a = (1 - ttbar + r**2) / 2
nu = (1 + ttbar - r**2) / (2 * sqrt(ttbar))
n1 = nint(val(n_cell$)) - 1

exp_half = nu + I_imaginary * sqrt(1 - nu**2)
exp_n = exp_half ** (2 * n1)
f_plus  = a - I_imaginary * sqrt(ttbar) * sqrt(1 - nu**2)
f_minus = a + I_imaginary * sqrt(ttbar) * sqrt(1 - nu**2)
Rc_n_top = r * (r - Rc_1_top * f_minus - (r - Rc_1_top * f_plus) * exp_n) / &
               (f_plus * (r - Rc_1_top * f_minus) - f_minus * (r - Rc_1_top * f_plus) * exp_n)

! Now factor in the air interface

k_a = kz_air
denom = k_a + k2

tt = 4 * k_a * k2 * (c2 / denom)**2
r_aa = (k_a - k2) / denom
r_22 = c2**2 * (k2 - k_a) / denom

R_tot = r_aa + tt * Rc_n_top / (1 - r_22 * Rc_n_top)

e_field = e_field * abs(R_tot)
e_phase = e_phase + atan2(aimag(R_tot), real(R_tot))

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
real(rp) sin_alpha, cos_alpha, sin_psi, cos_psi, wavelength
real(rp) cos_g, sin_g, cos_tc, sin_tc, tilt_cor_mat(3,3)
real(rp) h_norm(3), kh_outside_norm(3), e_tot, pc, p_factor
real(rp) m_out(3, 3), y_out(3)
real(rp) temp_vec(3), direction(3), s_len
real(rp) gamma_0, gamma_h, b_err
real(rp) hit_point(3)

logical curved_surface

!

wavelength = c_light * h_planck / end_orb%p0c

! (px, py, sqrt(1-px^2+py^2)) coords are with respect to laboratory reference trajectory.
! Convert this vector to k0_outside_norm which are coords with respect to crystal surface.
! k0_outside_norm is normalized to 1.

sin_g = sin(ele%value(graze_angle_in$))
cos_g = cos(ele%value(graze_angle_in$))

k0_outside_norm(1) =  cos_g * end_orb%vec(2) + sin_g * end_orb%vec(6)
k0_outside_norm(2) = end_orb%vec(4)
k0_outside_norm(3) = -sin_g * end_orb%vec(2) + cos_g * end_orb%vec(6)

! m_in 

m_in(1, 1:3) = [ cos_g, 0.0_rp, sin_g]
m_in(2, 1:3) = [0.0_rp, 1.0_rp, 0.0_rp]
m_in(3, 1:3) = [-sin_g, 0.0_rp, cos_g]

! If there is curvature, compute the reflection point which is where 
! the photon intersects the surface.

curved_surface = has_curved_surface(ele)
if (curved_surface) then
  call reflection_point (ele, cos_g, sin_g, end_orb, k0_outside_norm, hit_point, s_len)
  call to_curved_body_coords (ele, hit_point, k0_outside_norm, set$)
else
  hit_point = [end_orb%vec(1) * cos_g, end_orb%vec(3), -end_orb%vec(1) * sin_g]
endif

! Check aperture

if (ele%aperture_at == surface$) then
  end_orb%vec(1:5:2) = hit_point ! This is temporary
  call check_aperture_limit (end_orb, ele, surface$, param)
  if (end_orb%state /= alive$) return
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
c_param%dtheta_sin_2theta = -dot_product(h_norm + 2 * k0_outside_norm, h_norm) / 2
c_param%f0                = cmplx(ele%value(f0_re$), ele%value(f0_im$)) 
c_param%fh                = cmplx(ele%value(fh_re$), ele%value(fh_im$))

p_factor = cos(2.0_rp*ele%value(graze_angle_in$))
call e_field_calc (c_param, ele, p_factor, end_orb%e_field_x, end_orb%phase_x)
call e_field_calc (c_param, ele, 1.0_rp,   end_orb%e_field_y, end_orb%phase_y)   ! Sigma polarization

! Rotate back from curved body coords to element coords

if (curved_surface) then
  call to_curved_body_coords (ele, hit_point, kh_outside_norm, unset$)
  end_orb%vec(5) = end_orb%vec(5) + s_len  ! IS this correct ?!!!!!!!!!!!!!
endif

!--------------------------------- 
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

! m_out(:,1) = x_out = vector orthogonal to y and z

m_out(1,:) = [y_out(2) * ele%value(kh_z_norm$) - y_out(3) * ele%value(kh_y_norm$), &
             -y_out(1) * ele%value(kh_z_norm$) + y_out(3) * ele%value(kh_x_norm$), &
              y_out(1) * ele%value(kh_y_norm$) - y_out(2) * ele%value(kh_x_norm$)]
m_out(2,:) = y_out
m_out(3,:) = [ele%value(kh_x_norm$), ele%value(kh_y_norm$), ele%value(kh_z_norm$)]   ! k_out

!

direction = matmul(m_out, kh_outside_norm)
end_orb%vec(2) = direction(1)
end_orb%vec(4) = direction(2)
end_orb%vec(6) = sqrt(1 - end_orb%vec(2)**2 - end_orb%vec(4)**2)

! Compute position, backpropagating the ray

temp_vec = matmul(m_out, hit_point)
temp_vec = temp_vec - direction * temp_vec(3)

end_orb%vec(1) = temp_vec(1)
end_orb%vec(3) = temp_vec(2)
! %vec(5) doesn't include phase change due to wave nature of radiation
end_orb%vec(5) = temp_vec(3) 

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

real(rp) p_factor, e_field, e_phase, sqrt_b

complex(rp) e_rel, e_rel_a, e_rel_b, eta, eta1, f_cmp, xi_0k_a, xi_hk_a, xi_0k_b, xi_hk_b

! Construct xi_0k = xi_0 / k and xi_hk = xi_h / k

sqrt_b = sqrt(abs(cp%b_eff))

eta = (cp%b_eff * cp%dtheta_sin_2theta + cp%f0 * cp%cap_gamma * (1.0_rp - cp%b_eff)/2) / &
                                              (cp%cap_gamma * abs(p_factor) * sqrt_b * cp%fh) 
eta1 = sqrt(eta**2 + sign(1.0_rp, cp%b_eff))
f_cmp = abs(p_factor) * sqrt_b * cp%cap_gamma * cp%fh / 2

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

  ! Factor of sqrt_b comes from geometrical change in the transverse width of the photon beam

  e_field = e_field * abs(e_rel) / sqrt_b
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
! Subroutine reflection_point (ele, sin_g, cos_g, end_orb, k0_outside_norm, hit_point, s_len)
!
! Private routine to compute position where crystal reflects in crystal coordinates.
!
! Input:
!   ele                -- ele_struct: Element
!   cos_g              -- Real(rp): cosine of grazing angle.
!   sin_g              -- Real(rp): sine of grazing angle.
!   end_orb            -- coord_struct: Input coordinates
!   k0_outside_norm(3) -- Real(rp): Direction vector. 
!
! Output:
!   hit_point(3) -- Real(rp): point of reflection in crystal coordinates
!   s_len        -- Real(rp): distance for photon to propagate to hit_point
!-

subroutine reflection_point (ele, cos_g, sin_g, end_orb, k0_outside_norm, hit_point, s_len)

use nr, only: zbrent

implicit none

type (ele_struct), target :: ele
type (coord_struct), target :: end_orb

real(rp), target :: vec0(3), k0_outside_norm(3), cos_g, sin_g
real(rp) :: s_len, hit_point(3)
real(rp) :: s1, s2, s_center, fa, fb

! vec0 is the body element coordinates of the photon.
! Note: The photon position in element entrance coords is [vec(1), vec(3), 0].

vec0 = [end_orb%vec(1) * cos_g, end_orb%vec(3), -end_orb%vec(1) * sin_g]

! Assume flat crystal, compute s required to hit the intersection
! Choose a Bracket of 1m around this point.

if (ele%value(b_param$) < 0) then ! Bragg
  s_center = vec0(1) / k0_outside_norm(1)
else
  s_center = vec0(3) / k0_outside_norm(3)
endif

s1 = s_center
s2 = s_center
if (photon_depth_in_crystal(s_center) > 0) then
  do
    s1 = s_center - 0.1
    if (photon_depth_in_crystal(s1) < 0) exit
  enddo
else
  do
    s2 = s_center + 0.1
    if (photon_depth_in_crystal(s2) > 0) exit
  enddo
endif

s_len = zbrent (photon_depth_in_crystal, s1, s2, 1d-10)

! Compute the intersection point
hit_point = s_len * k0_outside_norm + vec0

contains

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
real(rp) :: point(3), r, y

point = s_len * k0_outside_norm + vec0

if (ele%value(b_param$) < 0) then ! Bragg
  r = point(3)
  delta_h = point(1) 
else
  r = point(1)
  delta_h = point(3) 
endif

y = point(2)

delta_h = delta_h + ele%value(a2_trans_curve$) * y**2 + ele%value(a3_trans_curve$) * y**3 + &
                    ele%value(a4_trans_curve$) * y**4 + ele%value(c2_curve_tot$) * r**2 + &
                    ele%value(c3_curve_tot$) * r**3 + ele%value(c4_curve_tot$) * r**4

end function photon_depth_in_crystal

end subroutine reflection_point

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

type (ele_struct), target :: ele

real(rp) hit_point(3), vector(3)
real(rp) curverot(3,3)
real(rp) slope_t, slope_c, cos_c, sin_c, cos_t, sin_t
real(rp), pointer :: v(:)
logical set

! The transformation curverot is
!     curve_rot = W_c W_t
! where W_c is the rotation in the reflection plane and W_t is the rotation
! in the transverse plane. [Note: There is no good reason for choosing W_c W_t over W_t W_c.]
! Want to have the transformation curverot match the slopes slope_c and slope_t.
! Since rotations do not commute, the cos_t and sin_t terms, which are used
! in W_t (the second rotation), are modified by cos_c.

v => ele%value

slope_c = 2.0_rp * v(c2_curve_tot$) * hit_point(3) + 3.0_rp * v(c3_curve_tot$) * hit_point(3)**2 + &
          4.0_rp * v(c4_curve_tot$) * hit_point(3)**3 
slope_t = 2.0_rp * v(a2_trans_curve$) * hit_point(2) + 3.0_rp * v(a3_trans_curve$) * hit_point(2)**2 + &
          4.0_rp * v(a4_trans_curve$) * hit_point(2)**3

if (.not. set) then
  slope_c = -slope_c
  slope_t = -slope_t
endif

cos_c =       1 / sqrt(1 + slope_c**2)
sin_c = slope_c / sqrt(1 + slope_c**2)

cos_t =               1 / sqrt(1 + (cos_c * slope_t)**2)
sin_t = cos_c * slope_t / sqrt(1 + (cos_c * slope_t)**2)

curverot(1, 1:3) = [ cos_c*cos_t,  cos_c*sin_t,  sin_c]
curverot(2, 1:3) = [      -sin_t,        cos_t, 0.0_rp]
curverot(3, 1:3) = [-sin_c*cos_t, -sin_c*sin_t,  cos_c]

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

if (ele%value(c2_curve_tot$) /= 0 .or. ele%value(c3_curve_tot$) /= 0 .or. ele%value(c4_curve_tot$) /= 0 .or. &
    ele%value(a2_trans_curve$) /= 0 .or. ele%value(a3_trans_curve$) /= 0 .or. ele%value(a4_trans_curve$) /= 0) then
  has_curve = .true.
else
  has_curve = .false.
endif

end function has_curved_surface

end module
