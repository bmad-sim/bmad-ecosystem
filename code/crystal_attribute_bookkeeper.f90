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
!   ele   -- Ele_struct: Crystal element.
!     %value(graze_angle_in$)
!     %value(graze_angle_out$)
!     %value(tilt_corr$)
!     %value(nx_out$
!     %value(ny_out$)
!     %value(nz_out$)
!-

subroutine crystal_attribute_bookkeeper (ele)

use bmad_struct

implicit none

type (ele_struct) ele

real(rp) lambda, gamma, delta1, lambda_in, d, alpha, psi, theta0
real(rp) cos_theta0, sin_theta0, graze_angle_in, ang_tot
real(rp) h_x, h_y, h_z, nx_out, ny_out, nz_out, nxx_out, nyy_out, nzz_out
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
  cos_graze_in = cos_theta0/delta1
  graze_angle_in = acos(cos_graze_in)
  sin_graze_in = sin(graze_angle_in)
else                              ! Laue
  sin_graze_in = sin_theta0/delta1
  graze_angle_in = asin(sin_graze_in)
  cos_graze_in = cos(graze_angle_in)
endif

ele%value(graze_angle_in$) = graze_angle_in

! (nx, ny, nz)_out is the normalized k_H vector.
! obtained from: k_H = H + K_0 + q*n_surface

if (ele%value(b_param$) < 0) then ! Bragg
  ny_out = lambda * h_y
  nz_out = lambda * h_z + cos_theta0 / delta1
  nx_out = -sqrt(1 - ny_out**2 - nz_out**2)
else                              ! Laue
  nx_out = lambda * h_x + sin_theta0 / delta1  
  ny_out = lambda * h_y
  nz_out = sqrt(1 - nx_out**2 - ny_out**2)
endif

ele%value(nx_out$) = nx_out
ele%value(ny_out$) = ny_out
ele%value(nz_out$) = nz_out

! (nxx, nyy, nzz)_out is (nx, ny, nz)_out in element entrance coordinates

nxx_out = nz_out * sin_graze_in - nx_out * cos_graze_in
nyy_out = ny_out
nzz_out = nz_out * cos_graze_in + nx_out * sin_graze_in

! There is no tilt correction if we are following the non-diffracted beam.

if (ele%value(follow_diffracted_beam$) == 0) then
  ele%value(tilt_corr$) = 0
else
  ele%value(tilt_corr$) = atan2(nyy_out, nxx_out)
endif

! total graze angle

ang_tot = atan2(sqrt(nxx_out**2 + nyy_out**2), nzz_out)
ele%value(graze_angle_out$) = ang_tot - graze_angle_in

! displacement L due to finite crystal thickness for Laue diffraction.

if (ele%value(b_param$) > 0) then
  ! Energy flow direction is K_0 + K_H = 2*K_0 + H
  s_vec = 2 * [sin_theta0, 0.0_rp, cos_theta0] / (delta1 * lambda) + [h_x, h_y, h_z] 
  ele%value(l_x$:l_x$+2) = dot_product(s_vec, [0.0_rp, 0.0_rp, ele%value(thickness$)]) * &
                              s_vec / dot_product(s_vec, s_vec)
endif

! Attributes for to curved crystal
! source_r and detect_r are reciprocrals of distance multipled by angle

if ( ele%value(d_source$) == 0) then
  source_r = 0.0_rp
else
  source_r=ele%value(graze_angle_in$)/ele%value(d_source$)
endif

if ( ele%value(d_detec$) == 0) then
  detect_r = 0.0_rp
else
  detect_r = ele%value(graze_angle_out$)/ele%value(d_detec$)
endif

ele%value(c2_curve_tot$)=ele%value(c2_curve$)+(detect_r+source_r)*0.25_rp
ele%value(c4_curve_tot$)=ele%value(c4_curve$)+0.03125_rp*(detect_r+source_r)**3


end subroutine
