!+
! Subroutine strong_beam_sigma_calc (ele, s_pos, sig_x, sig_y, bbi_const)
!
! Routine to calculate the strong beam sigmas, strong beam centroid offsets (due to crabbing), and
! the BBI force constant for a beambeam element.
!
! Input:
!   ele           -- ele_struct: Beambeam element.
!   s_pos         -- real(rp): Longitudinal position in lab coords of slice (used with hourglass effect correction).
!
! Output:
!   sig_x, sig_y       -- real(rp): Strong beam sigmas.
!   bbi_const          -- real(rp): BBI kick scale factor.
!-

subroutine strong_beam_sigma_calc (ele, s_pos, sig_x, sig_y, bbi_const)

use bmad_routine_interface, dummy => strong_beam_sigma_calc

implicit none

type (ele_struct) ele
type (ele_struct), pointer :: ele0

real(rp) s_pos, bbi_const, x_center, y_center, sig_x, sig_y
real(rp) beta, sig_x0, sig_y0, beta_a0, beta_b0, alpha_a0, alpha_b0

!

sig_x0 = ele%value(sig_x$)
sig_y0 = ele%value(sig_y$)

if (ele%value(beta_a_strong$) == 0) then
  ele0 => pointer_to_next_ele(ele, -1)
  beta_a0 = ele%a%beta
  alpha_a0 = 0.5_rp * (ele0%a%alpha + ele%a%alpha)
else
  beta_a0 = ele%value(beta_a_strong$)
  alpha_a0 = ele%value(alpha_a_strong$)
endif

if (ele%value(beta_b_strong$) == 0) then
  ele0 => pointer_to_next_ele(ele, -1)
  beta_b0 = ele%b%beta
  alpha_b0 = 0.5_rp * (ele0%b%alpha + ele%b%alpha)
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

bbi_const = -strong_beam_strength(ele) * classical_radius_factor /  &
                                      (2 * pi * ele%value(p0c$) * (sig_x + sig_y))

end subroutine
