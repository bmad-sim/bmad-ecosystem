module track1_photon_mod

use bmad_struct
use bmad_interface
use track1_mod

! This is for passing info into the field_calc_routine

type, private :: crystal_param_struct
  real(rp) cap_gamma, dtheta_sin_2theta, b_eff
  complex(rp) f0, fh
end type

private e_field_calc

contains

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine track1_mirror (ele, param, orbit)
!
! Routine to track reflection from a mirror.
!
! Input:
!   ele    -- ele_struct: Element tracking through.
!   param  -- lat_param_struct: lattice parameters.
!   orbit    -- Coord_struct: phase-space coords to be transformed
!
! Output:
!   orbit    -- Coord_struct: final phase-space coords
!-

subroutine track1_mirror (ele, param, orbit)

implicit none

type (ele_struct), target:: ele
type (coord_struct), target:: orbit
type (lat_param_struct) :: param

real(rp) wavelength
real(rp), pointer :: val(:)

character(32), parameter :: r_name = 'track1_mirror'

!

val => ele%value
wavelength = c_light * h_planck / orbit%p0c

call to_surface_coords (ele, orbit)

! Check aperture

if (ele%aperture_at == surface$) then
  call check_aperture_limit (orbit, ele, surface$, param)
  if (orbit%state /= alive$) return
endif

! Reflect momentum vector

orbit%vec(2) = -orbit%vec(2)

! Rotate back to uncurved element coords

if (ele%surface%has_curvature) then
  call to_curved_body_coords (ele, orbit, unset$)
endif

end subroutine track1_mirror

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine track1_multilayer_mirror (ele, param, orbit)
!
! Routine to track reflection from a multilayer_mirror.
! Basic equations are from Kohn, "On the Theory of Reflectivity of an X-Ray Multilayer Mirror".
!
! Input:
!   ele    -- ele_struct: Element tracking through.
!   param  -- lat_param_struct: lattice parameters.
!   orbit    -- Coord_struct: phase-space coords to be transformed
!
! Output:
!   orbit    -- Coord_struct: final phase-space coords
!-

subroutine track1_multilayer_mirror (ele, param, orbit)

implicit none

type (ele_struct), target:: ele
type (coord_struct), target:: orbit
type (lat_param_struct) :: param

real(rp) wavelength, kz_air
real(rp), pointer :: val(:)

complex(rp) zero, xi_1, xi_2, kz1, kz2, c1, c2

character(32), parameter :: r_name = 'track1_multilayer_mirror'

!

val => ele%value
wavelength = c_light * h_planck / orbit%p0c

call to_surface_coords (ele, orbit)

! Check aperture

if (ele%aperture_at == surface$) then
  call check_aperture_limit (orbit, ele, surface$, param)
  if (orbit%state /= alive$) return
endif

! Note: Koln z-axis = Bmad x-axis.
! Note: f0_re and f0_im are both positive.

xi_1 = cmplx(-val(f0_re1$), val(f0_im1$)) * r_e * wavelength**2 / (pi * val(v1_unitcell$)) 
xi_2 = cmplx(-val(f0_re2$), val(f0_im2$)) * r_e * wavelength**2 / (pi * val(v2_unitcell$)) 

kz1 = twopi * sqrt(orbit%vec(2)**2 + xi_1) / wavelength
kz2 = twopi * sqrt(orbit%vec(2)**2 + xi_2) / wavelength
kz_air = twopi * orbit%vec(2) / wavelength

c1 = exp(I_imaginary * kz1 * val(d1_thickness$) / 2)
c2 = exp(I_imaginary * kz2 * val(d2_thickness$) / 2)

zero = cmplx(0.0_rp, 0.0_rp)

call multilayer_track (xi_1, xi_2, orbit%field(1), orbit%phase(1))     ! pi polarization
call multilayer_track (zero, zero, orbit%field(2), orbit%phase(2))     ! sigma polarization

! Reflect momentum vector

orbit%vec(2) = -orbit%vec(2)

! Rotate back to uncurved element coords

if (ele%surface%has_curvature) then
  call to_curved_body_coords (ele, orbit, unset$)
endif

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
! Subroutine track1_crystal (ele, param, orbit)
!
! Routine to track reflection from a crystal.
!
! Input:
!   ele      -- ele_struct: Element tracking through.
!   param    -- lat_param_struct: lattice parameters.
!   orbit    -- Coord_struct: phase-space coords to be transformed
!
! Output:
!   orbit    -- Coord_struct: final phase-space coords
!-

subroutine track1_crystal (ele, param, orbit)

implicit none

type (ele_struct), target:: ele
type (coord_struct), target:: orbit
type (lat_param_struct) :: param
type (crystal_param_struct) c_param

real(rp) h_norm(3), e_tot, pc, p_factor, wavelength
real(rp) gamma_0, gamma_h, old_vec(6)

character(*), parameter :: r_name = 'track1_cyrstal'

! A graze angle of zero means the wavelength of the reference photon was too large
! for the bragg condition. 

if (ele%value(bragg_angle_in$) == 0) then
  call out_io (s_fatal$, r_name, 'REFERENCE ENERGY TOO SMALL TO SATISFY BRAGG CONDITION!')
  orbit%state = lost_z_aperture$
  if (global_com%exit_on_error) call err_exit
  return
endif

!

wavelength = c_light * h_planck / orbit%p0c

! (px, py, pz) coords are with respect to laboratory reference trajectory.
! Convert this vector to k0_outside_norm which are coords with respect to crystal surface.
! k0_outside_norm is normalized to 1.

call to_surface_coords (ele, orbit)
old_vec = orbit%vec

! Check aperture

if (ele%aperture_at == surface$) then
  call check_aperture_limit (orbit, ele, surface$, param)
  if (orbit%state /= alive$) return
endif

! Construct h_norm = H * wavelength.

h_norm = [-ele%value(h_x_norm$), ele%value(h_y_norm$), ele%value(h_z_norm$)] * wavelength / ele%value(d_spacing$)

! kh_outside_norm is the normalized outgoing wavevector outside the crystal

c_param%cap_gamma = r_e * wavelength**2 / (pi * ele%value(v_unitcell$)) 

orbit%vec(2:6:2) = orbit%vec(2:6:2) + h_norm

if (ele%value(b_param$) < 0) then ! Bragg
  orbit%vec(2) = -sqrt(1 - orbit%vec(4)**2 - orbit%vec(6)**2)
  gamma_0 = old_vec(2)
  gamma_h = orbit%vec(2)
else
  orbit%vec(6) = sqrt(1 - orbit%vec(2)**2 - orbit%vec(4)**2)
  gamma_0 = old_vec(6)
  gamma_h = orbit%vec(6)
endif

!-------------------------------------
! Calculate phase and intensity

c_param%b_eff             = gamma_0 / gamma_h
c_param%dtheta_sin_2theta = -dot_product(h_norm + 2 * old_vec(2:6:2), h_norm) / 2
c_param%f0                = cmplx(ele%value(f0_re$), ele%value(f0_im$)) 
c_param%fh                = cmplx(ele%value(fh_re$), ele%value(fh_im$))

p_factor = cos(ele%value(bragg_angle_in$) + ele%value(bragg_angle_out$))
call e_field_calc (c_param, ele, p_factor, orbit%field(1), orbit%phase(1))
call e_field_calc (c_param, ele, 1.0_rp,   orbit%field(2), orbit%phase(2))   ! Sigma polarization

! Rotate back from curved body coords to element coords

if (ele%surface%has_curvature) then
  call to_curved_body_coords (ele, orbit, unset$)
endif

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
! Subroutine to_surface_coords (ele, orbit)
!
! Routine to adjust the photon position to be at the surface
!
! Input:
!   ele                -- ele_struct: Element
!   orbit              -- coord_struct: Input coordinates
!
! Output:
!   orbit              -- coord_struct: Input coordinates
!-

subroutine to_surface_coords (ele, orbit)

use nr, only: zbrent

implicit none

type (ele_struct) ele
type (coord_struct) orbit

real(rp) :: s_len, s1, s2, s_center


! If there is curvature, compute the reflection point which is where 
! the photon intersects the surface.

if (ele%surface%has_curvature) then

  ! Assume flat crystal, compute s required to hit the intersection
  ! Choose a Bracket of 1m around this point.

  if (ele%value(b_param$) < 0) then ! Bragg
    s_center = orbit%vec(1) / orbit%vec(2)
  else
    s_center = orbit%vec(5) / orbit%vec(6)
  endif

  s1 = s_center
  s2 = s_center
  if (photon_depth_in_crystal(s_center) > 0) then
    do
      s1 = s1 - 0.1
      if (photon_depth_in_crystal(s1) < 0) exit
    enddo
  else
    do
      s2 = s2 + 0.1
      if (photon_depth_in_crystal(s2) > 0) exit
    enddo
  endif

  s_len = zbrent (photon_depth_in_crystal, s1, s2, 1d-10)

  ! Compute the intersection point

  orbit%vec(1:5:2) = s_len * orbit%vec(2:6:2) + orbit%vec(1:5:2)
  orbit%t = orbit%t + s_len / c_light

  call to_curved_body_coords (ele, orbit, set$)

else
  s_len = -orbit%vec(1) / orbit%vec(2)
  orbit%vec(1:5:2) = orbit%vec(1:5:2) + s_len * orbit%vec(2:6:2) ! Surface is at x = 0
  orbit%t = orbit%t + s_len / c_light
endif

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

type (photon_surface_struct), pointer :: surface

real(rp), intent(in) :: s_len
real(rp) :: delta_h
real(rp) :: point(3), z, y
integer iz, iy

!

point = s_len * orbit%vec(2:6:2) + orbit%vec(1:5:2)

if (ele%value(b_param$) < 0) then ! Bragg
  z = point(3)
  delta_h = point(1)
else
  z = point(1)
  delta_h = point(3)
endif

y = point(2)
surface => ele%surface

do iz = 0, ubound(surface%curvature_zy, 1)
do iy = 0, ubound(surface%curvature_zy, 2) - iz
  if (ele%surface%curvature_zy(iz, iy) == 0) cycle
  delta_h = delta_h + surface%curvature_zy(iz, iy) * y**iy * z**iz
enddo
enddo

end function photon_depth_in_crystal

end subroutine to_surface_coords

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine to_curved_body_coords (ele, orbit, set)
!
! Routine to rotate between element body coords and effective body coords ("curved body coords") with 
! respect to the surface at the point of photon impact.
!
! Input:
!   ele      -- ele_struct: reflecting element
!   orbit    -- coord_struct: Photon position.
!   set      -- Logical: True -> Transform body to curved body. 
!                        False -> Transform curved body to body.
!
! Output:
!   orbit    -- coord_struct: Photon position.
!-

! Compute the slope of the crystal at that the point of impact.
! curverot transforms from standard body element coords to body element coords at point of impact.

subroutine to_curved_body_coords (ele, orbit, set)

implicit none

type (ele_struct), target :: ele
type (coord_struct) orbit
type (photon_surface_struct), pointer :: s

real(rp) curverot(3,3), angle
real(rp) slope_t, slope_c, cos_c, sin_c, cos_t, sin_t, y, z
integer iz, iy

logical set

! Compute slopes

s => ele%surface

slope_c = 0
slope_t = 0
y = orbit%vec(3)
z = orbit%vec(5)

do iz = 0, ubound(s%curvature_zy, 1)
do iy = 0, ubound(s%curvature_zy, 2) - iz
  if (s%curvature_zy(iz, iy) == 0) cycle
  if (iz > 0) slope_c = slope_c - iz * s%curvature_zy(iz, iy) * y**iy * z**(iz-1)
  if (iy > 0) slope_t = slope_t - iy * s%curvature_zy(iz, iy) * y**(iy-1) * z**iz
enddo
enddo

if (slope_c == 0 .and. slope_t == 0) return

! Compute rotation matrix and goto body element coords at point of photon impact

angle = atan2(sqrt(slope_c**2 + slope_t**2), 1.0_rp)
if (set) angle = -angle
call axis_angle_to_w_mat ([0.0_rp, slope_c, -slope_t], angle, curverot)

orbit%vec(2:6:2) = matmul(curverot, orbit%vec(2:6:2))

end subroutine to_curved_body_coords

end module
