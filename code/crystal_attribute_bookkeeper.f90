subroutine crystal_attribute_bookkeeper (ele)

use bmad_struct

implicit none

type (ele_struct) ele

real(rp) lambda, gamma, delta1, lambda_in, d, sin_alpha, cos_alpha, sin_psi, cos_psi, theta0
real(rp) cos_theta0, sin_theta0, graze_angle_in
real(rp) kx_in

! If the photon energy is not set then cannot do the calc yet.

if (ele%value(e_tot$) == 0) return

!

lambda = ele%value(ref_wave_length$)
gamma = lambda**2 * r_e / (pi * ele%value(v_unitcell$))
delta1 = 1 + gamma * ele%value(f0_re$) / 2
lambda_in = lambda * delta1
d = ele%value(d_spacing$)

sin_alpha = sin(ele%value(alpha_angle$))
cos_alpha = cos(ele%value(alpha_angle$))

sin_psi = sin(ele%value(psi_angle$))
cos_psi = cos(ele%value(psi_angle$))

theta0 = asin(lambda_in / (2 * d * sqrt(1 - (sin_alpha * sin_psi)**2))) - atan(sin_alpha * cos_psi / cos_alpha)
cos_theta0 = cos(theta0)
sin_theta0 = sin(theta0)

graze_angle_in = acos(cos_theta0/delta1)
ele%value(graze_angle_in$) = graze_angle_in

kx_in = 0

end subroutine
