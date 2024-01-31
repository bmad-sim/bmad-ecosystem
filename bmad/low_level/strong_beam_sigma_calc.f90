!+
! Subroutine strong_beam_sigma_calc (ele, s_pos, sigma, bbi_const, dsigma_ds)
!
! Routine to calculate the strong beam sigmas, strong beam centroid offsets (due to crabbing), and
! the BBI force constant for a beambeam element.
!
! Input:
!   ele             -- ele_struct: Beambeam element.
!   s_pos           -- real(rp): Longitudinal position in lab coords of slice (used with hourglass effect correction).
!
! Output:
!   sigma(2)        -- real(rp): Strong beam x,y sigmas.
!   bbi_const       -- real(rp): BBI kick scale factor.
!   dsigma_ds(2)    -- real(rp): sig_x and sig_y longitudinal derivatives.
!-

subroutine strong_beam_sigma_calc (ele, s_pos, sigma, bbi_const, dsigma_ds)

use bmad_routine_interface, dummy => strong_beam_sigma_calc

implicit none

type (ele_struct) ele
type (ele_struct), pointer :: ele0

real(rp) s_pos, bbi_const, x_center, y_center, sigma(2), dsigma_ds(2)
real(rp) beta, sig_x0, sig_y0, beta_a0, beta_b0, alpha_a0, alpha_b0, gamma0

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
  sigma = [sig_x0, sig_y0]
  dsigma_ds = 0
else
  gamma0 = (1 + alpha_a0**2) / beta_a0
  beta = beta_a0 - 2 * alpha_a0 * s_pos + gamma0 * s_pos**2
  sigma(1) = sig_x0 * sqrt(beta / beta_a0)
  dsigma_ds(1) = -(alpha_a0 - s_pos * gamma0) * sigma(1) / beta

  gamma0 = (1 + alpha_b0**2) / beta_b0
  beta = beta_b0 - 2 * alpha_b0 * s_pos + gamma0 * s_pos**2
  sigma(2) = sig_y0 * sqrt(beta / beta_b0)
  dsigma_ds(2) = -(alpha_b0 - s_pos * gamma0) * sigma(2) / beta
endif

bbi_const = -strong_beam_strength(ele) * classical_radius_factor /  &
                                      (2 * pi * ele%value(p0c$) * (sigma(1) + sigma(2)))

end subroutine
