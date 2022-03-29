!+
! Subroutine strong_beam_sigma_calc (ele, s_pos, z_strong, sig_x, sig_y, bbi_const, x_center, y_center)
!
! Routine to calculate the strong beam sigmas, strong beam centroid offsets (due to crabbing), and
! the BBI force constant for a beambeam element.
!
! Input:
!   ele           -- ele_struct: Beambeam element.
!   s_pos         -- real(rp): Longitudinal position.
!   z_strong      -- real(rp): Position within strong beam. Positive z_strong is at the tail of the bunch
!
! Output:
!   sig_x, sig_y  -- real(rp): Strong beam sigmas.
!   bbi_const
!   x_center, y_center -- real(rp): Strong beam centroid

subroutine strong_beam_sigma_calc (ele, s_pos, z_strong, sig_x, sig_y, bbi_const, x_center, y_center)

use bmad_struct

implicit none

type (ele_struct) ele

real(rp) s_pos, z_strong, bbi_const, x_center, y_center, sig_x, sig_y
real(rp) beta, sig_x0, sig_y0, beta_a0, beta_b0, alpha_a0, alpha_b0, r

!

sig_x0 = ele%value(sig_x$)
sig_y0 = ele%value(sig_y$)

if (ele%value(beta_a_strong$) == 0) then
  beta_a0 = ele%a%beta
  alpha_a0 = ele%a%alpha
else
  beta_a0 = ele%value(beta_a_strong$)
  alpha_a0 = ele%value(alpha_a_strong$)
endif

if (ele%value(beta_b_strong$) == 0) then
  beta_b0 = ele%b%beta
  alpha_b0 = ele%b%alpha
else
  beta_b0 = ele%value(beta_b_strong$)
  alpha_b0 = ele%value(alpha_b_strong$)
endif

if (beta_a0 == 0) then
  sig_x = sig_x0
  sig_y = sig_y0
else
  beta = beta_a0 - 2 * alpha_a0 * s_pos + (1 + alpha_a0**2) * s_pos**2 / beta_a0
  sig_x = sig_x0 * sqrt(beta / beta_a0)
  beta = beta_b0 - 2 * alpha_b0 * s_pos + (1 + alpha_b0**2) * s_pos**2 / beta_b0
  sig_y = sig_y0 * sqrt(beta / beta_b0)
endif

bbi_const = -ele%branch%param%n_part * ele%value(charge$) * classical_radius_factor /  &
                                      (2 * pi * ele%value(p0c$) * (sig_x + sig_y))

r = ((((ele%value(crab_x5$) * z_strong + ele%value(crab_x4$)) * z_strong + ele%value(crab_x3$)) * z_strong + & 
                                         ele%value(crab_x2$)) * z_strong + ele%value(crab_x1$)) * z_strong

x_center = r * cos(ele%value(crab_tilt$))
y_center = r * sin(ele%value(crab_tilt$))

end subroutine
