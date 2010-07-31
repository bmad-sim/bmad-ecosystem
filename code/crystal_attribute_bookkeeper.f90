subroutine crystal_attribute_bookkeeper (ele)

use bmad_struct

implicit none

type (ele_struct) ele

real(rp) lambda, gamma, delta1, lambda_in, d, alpha, psi, theta0
real(rp) cos_theta0, sin_theta0, graze_angle_in
real(rp) h_x, h_y, h_z, kx_in

! If the photon energy is not set then cannot do the calc yet.

if (ele%value(e_tot$) == 0) return

!

lambda = ele%value(ref_wave_length$)
gamma = lambda**2 * r_e / (pi * ele%value(v_unitcell$))
delta1 = 1 + gamma * ele%value(f0_re$) / 2
lambda_in = lambda * delta1
d = ele%value(d_spacing$)

alpha = ele%value(alpha_angle$)
psi   = ele%value(psi_angle$)

h_x = -cos(alpha) / d
h_y = sin(alpha) * sin(psi) / d
h_z = sin(alpha) * cos(psi) / d

theta0 = asin(lambda_in / (2 * sqrt(d**2 - h_y**2))) - atan(h_z/h_x)
cos_theta0 = cos(theta0)
sin_theta0 = sin(theta0)

graze_angle_in = acos(cos_theta0/delta1)
ele%value(graze_angle_in$) = graze_angle_in

ny_out =  lambda * h_y
nz_out = -lambda * h_z + cos_theta0 / delta1
nx_out = -sqrt(1 - ny_out**2 - nz_out**2)




end subroutine
