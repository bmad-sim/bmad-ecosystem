!+
! Subroutine crystal_attribute_bookkeeper (ele)
! 
! Routine to orient the crystal element.
!
! Input:
!   ele   -- Ele_struct: Crystal element.
!     %value(bragg_angle$)     -- Bragg angle.
!     %value(e_tot$)           -- Photon reference energy.
!     %value(alpha_angle$)     -- H-vector orientation
!     %value(d_spacing$)       -- Crystal plane spacing.
!     ... etc ...
!
! Output
!   ele       -- Ele_struct: Crystal element.
!     %value(graze_angle_in$)
!     %value(graze_angle_out$)
!     %value(tilt_corr$)
!     %value(kh_x_norm$
!     %value(kh_y_norm$)
!     %value(kh_z_norm$)
!     ... etc.
!-

subroutine crystal_attribute_bookkeeper (ele)

use bmad_struct

implicit none

type (ele_struct) ele

real(rp) lambda, gamma, delta1, lambda_in, d, alpha, psi, theta0
real(rp) cos_theta0, sin_theta0, graze_angle_in, ang_tot
real(rp) h_x, h_y, h_z, kh_x_norm, kh_y_norm, kh_z_norm, ent_kh_x_norm, ent_kh_y_norm, ent_kh_z_norm
real(rp) cos_graze_in, sin_graze_in, s_vec(3)
real(rp) source_r,detect_r,r_c

! If the photon energy or the bragg angle has not been set then cannot do the calc yet.

if (ele%value(e_tot$) == 0) return
if (ele%value(bragg_angle$) == 0) return

lambda = ele%value(ref_wavelength$)
gamma = lambda**2 * r_e / (pi * ele%value(v_unitcell$))
delta1 = 1 / sqrt( 1 - gamma * ele%value(f0_re$) )
lambda_in = lambda * delta1
d = ele%value(d_spacing$)

alpha = ele%value(alpha_angle$)
psi   = ele%value(psi_angle$)

h_x = -cos(alpha) / d
h_y = sin(alpha) * sin(psi) / d
h_z = sin(alpha) * cos(psi) / d

theta0 = asin(lambda_in / (2 * d * sqrt(1 - (d * h_y)**2))) - atan(h_z/h_x)
cos_theta0 = cos(theta0)
sin_theta0 = sin(theta0)

! Compute graze_angle_in

if (ele%value(b_param$) < 0) then ! Bragg
  cos_graze_in   = cos_theta0/delta1
  graze_angle_in = acos(cos_graze_in)
  sin_graze_in   = sin(graze_angle_in)
else                              ! Laue
  sin_graze_in   = sin_theta0/delta1
  graze_angle_in = asin(sin_graze_in)
  cos_graze_in   = cos(graze_angle_in)
endif

ele%value(graze_angle_in$) = graze_angle_in

! kh_norm is the normalized k_H vector.
! Obtained from: k_H = H + K_0 + q*n_surface

if (ele%value(b_param$) < 0) then ! Bragg
  kh_y_norm = lambda * h_y
  kh_z_norm = lambda * h_z + cos_theta0 / delta1
  kh_x_norm = -sqrt(1 - kh_y_norm**2 - kh_z_norm**2)
else                              ! Laue
  kh_x_norm = lambda * h_x + sin_theta0 / delta1  
  kh_y_norm = lambda * h_y
  kh_z_norm = sqrt(1 - kh_x_norm**2 - kh_y_norm**2)
endif

ele%value(kh_x_norm$) = kh_x_norm
ele%value(kh_y_norm$) = kh_y_norm
ele%value(kh_z_norm$) = kh_z_norm

! ent_kh_norm is kn_norm in element entrance coordinates

ent_kh_x_norm = kh_z_norm * sin_graze_in - kh_x_norm * cos_graze_in
ent_kh_y_norm = kh_y_norm
ent_kh_z_norm = kh_z_norm * cos_graze_in + kh_x_norm * sin_graze_in

! There is no tilt correction if we are following the non-diffracted beam.

if (nint(ele%value(ref_orbit_follows$)) == undiffracted$) then
  ele%value(tilt_corr$) = 0
else
  ele%value(tilt_corr$) = atan2(ent_kh_y_norm, ent_kh_x_norm)
endif

! total graze angle

if (nint(ele%value(ref_orbit_follows$)) == undiffracted$) then
  ele%value(graze_angle_out$) = -graze_angle_in
else
  ang_tot = atan2(sqrt(ent_kh_x_norm**2 + ent_kh_y_norm**2), ent_kh_z_norm)
  ele%value(graze_angle_out$) = ang_tot - graze_angle_in
endif

! displacement L due to finite crystal thickness for Laue diffraction.

if (ele%value(b_param$) > 0) then
  ! Energy flow direction is K_0 + K_H = 2*K_0 + H
  s_vec = 2 * [sin_theta0, 0.0_rp, cos_theta0] / (delta1 * lambda) + [h_x, h_y, h_z] 
  ele%value(l_x$:l_z$) = dot_product(s_vec, [0.0_rp, 0.0_rp, ele%value(thickness$)]) * &
                              s_vec / dot_product(s_vec, s_vec)
endif

!

ele%value(ref_cap_gamma$) = gamma
ele%value(darwin_width_sigma$) = 2 * gamma * ele%value(fh_re$) / &
                      (abs(sin(2 * ele%value(graze_angle$))) * sqrt(abs(ele%value(b_param$))))
ele%value(darwin_width_pi$) = ele%value(darwin_width_sigma$) * abs(cos(2 * ele%value(graze_angle$)))

end subroutine
