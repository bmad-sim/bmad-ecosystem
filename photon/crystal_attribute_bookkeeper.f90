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
!     %value(bragg_angle_in$)
!     %value(bragg_angle_out$)
!     %value(tilt_corr$)
!     ... etc.
!-

subroutine crystal_attribute_bookkeeper (ele)

use bmad_struct

implicit none

type (ele_struct) ele

real(rp) lambda, gamma, delta1, lambda_in, d, alpha, psi, theta0
real(rp) bragg_angle_in, ang_tot, k0_x_norm, k0_y_norm, k0_z_norm
real(rp) kh_x_norm, kh_y_norm, kh_z_norm, h_x_norm, h_y_norm, h_z_norm
real(rp) ent_kh_x_norm, ent_kh_y_norm, ent_kh_z_norm
real(rp) cos_graze_in, sin_graze_in, s_vec(3), beta
real(rp) source_r,detect_r,r_c, total_angle

! If the photon energy or the bragg angle has not been set then cannot do the calc yet.

if (ele%value(e_tot$) == 0) return
if (ele%value(bragg_angle$) == 0) return

lambda = ele%value(ref_wavelength$)
gamma = lambda**2 * r_e / (pi * ele%value(v_unitcell$))
delta1 = 1 / sqrt(1 - gamma * real(ele%photon%material%f_0))
lambda_in = lambda * delta1
d = ele%value(d_spacing$)

alpha = ele%value(alpha_angle$)
psi   = ele%value(psi_angle$)

h_x_norm = -sin(alpha) * cos(psi)
h_y_norm =  sin(alpha) * sin(psi) 
h_z_norm = -cos(alpha) 

ele%photon%material%h_norm = [h_x_norm, h_y_norm, h_z_norm]

! Compute bragg_angle_in

beta = lambda_in / (2 * d)

if (ele%value(b_param$) < 0) then ! Bragg
  theta0 = asin((-beta * h_z_norm - h_x_norm * sqrt(h_x_norm**2 + h_z_norm**2 - beta**2)) / (h_x_norm**2 + h_z_norm**2))
  bragg_angle_in = acos(cos(theta0)/delta1)
  k0_x_norm = -cos(bragg_angle_in)
  k0_z_norm =  sin(bragg_angle_in)

else                              ! Laue
  theta0 = asin((-beta * h_x_norm + h_z_norm * sqrt(h_x_norm**2 + h_z_norm**2 - beta**2)) / (h_x_norm**2 + h_z_norm**2))
  bragg_angle_in = asin(sin(theta0)/delta1)
  k0_x_norm =  sin(bragg_angle_in)
  k0_z_norm =  cos(bragg_angle_in)
endif

ele%value(bragg_angle_in$) = bragg_angle_in

! kh_norm is the normalized k_H vector.
! Obtained from: k_H = H + K_0 + q*n_surface

kh_x_norm = lambda * h_x_norm/d + k0_x_norm
kh_y_norm = lambda * h_y_norm/d
kh_z_norm = sqrt(1 - kh_x_norm**2 - kh_y_norm**2)
if (ele%value(b_param$) < 0) kh_z_norm = -kh_z_norm   ! Bragg

! ent_kh_norm is kh_norm in entrance coordinates 

ent_kh_x_norm = kh_z_norm * k0_x_norm - kh_x_norm * k0_z_norm
ent_kh_y_norm = kh_y_norm
ent_kh_z_norm = kh_z_norm * k0_z_norm + kh_x_norm * k0_x_norm

! There is no tilt correction if we are following the non-diffracted beam.

if (nint(ele%value(ref_orbit_follows$)) == bragg_diffracted$) then
  ele%value(tilt_corr$) = atan2(ent_kh_y_norm, ent_kh_x_norm)
else
  ele%value(tilt_corr$) = 0
endif

ang_tot = atan2(sqrt(ent_kh_x_norm**2 + ent_kh_y_norm**2), ent_kh_z_norm)
ele%value(bragg_angle_out$) = ang_tot - bragg_angle_in

! Displacement due to finite crystal thickness for Laue diffraction.

if (ele%value(b_param$) > 0) then
  select case (nint(ele%value(ref_orbit_follows$)))
  case (forward_diffracted$, undiffracted$)
    ! reference orbit direction is same as k_0 (outside)
    ele%photon%material%l_ref = [tan(bragg_angle_in), 0.0_rp, 1.0_rp] * ele%value(thickness$) 
  case (bragg_diffracted$)
    ! Energy flow direction is K_0 + K_H = 2*K_0 + H
    s_vec = 2 * [sin(theta0), 0.0_rp, cos(theta0)] / (delta1 * lambda) + [h_x_norm, h_y_norm, h_z_norm] / d 
    ele%photon%material%l_ref = s_vec * ele%value(thickness$) / s_vec(3) 
  end select
else
  ele%photon%material%l_ref = 0
endif

ele%value(l$) = norm2(ele%photon%material%l_ref)

!

ele%value(ref_cap_gamma$) = gamma
total_angle = ele%value(bragg_angle_in$) + ele%value(bragg_angle_out$)
ele%value(darwin_width_sigma$) = 2 * gamma * real(ele%photon%material%f_hkl) / &
                      (abs(sin(total_angle)) * sqrt(abs(ele%value(b_param$))))
ele%value(darwin_width_pi$) = ele%value(darwin_width_sigma$) * abs(cos(total_angle))

end subroutine
