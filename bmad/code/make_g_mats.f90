!+
! Subroutine make_g_mats (ele, g_mat, g_inv_mat)
!
! Subroutine make the matrices needed to go from normal mode coords 
! to coordinates with the beta function removed.
!
! Input:
!   ele        -- Ele_struct: Element
!
! Output:
!   g_mat(4,4)     -- Real(rp): Normal mode to betaless coords
!   g_inv_mat(4,4) -- Real(rp): The inverse of G_MAT
!-

subroutine make_g_mats (ele, g_mat, g_inv_mat)

  use bmad_interface, except_dummy => make_g_mats

  implicit none

  type (ele_struct) ele

  real(rp) g_mat(4,4), g_inv_mat(4,4)
  real(rp) sqrt_beta_a, sqrt_beta_b, alpha_a, alpha_b
!

  sqrt_beta_a = sqrt(ele%a%beta)
  alpha_a     = ele%a%alpha
  sqrt_beta_b = sqrt(ele%b%beta)
  alpha_b     = ele%b%alpha

  g_mat = 0
  g_mat(1,1) = 1 / sqrt_beta_a
  g_mat(2,1) = alpha_a / sqrt_beta_a
  g_mat(2,2) = sqrt_beta_a
  g_mat(3,3) = 1 / sqrt_beta_b
  g_mat(4,3) = alpha_b / sqrt_beta_b
  g_mat(4,4) = sqrt_beta_b
                                  
  g_inv_mat = 0
  g_inv_mat(1,1) = sqrt_beta_a
  g_inv_mat(2,1) = -alpha_a / sqrt_beta_a
  g_inv_mat(2,2) = 1 / sqrt_beta_a
  g_inv_mat(3,3) = sqrt_beta_b
  g_inv_mat(4,3) = -alpha_b / sqrt_beta_b
  g_inv_mat(4,4) = 1 / sqrt_beta_b

end subroutine
