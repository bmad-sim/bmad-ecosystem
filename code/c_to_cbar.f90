!+
! subroutine c_to_cbar (ele, cbar_mat)
!
! Subroutine to compute Cbar from the C matrix and the Twiss parameters.
!
! Modules Needed:
!   use bmad
!
! Input:
!   ele -- Ele_struct: Element with C matrix and Twiss parameters.
!
! Output:
!   cbar_mat(2,2) -- Real(rp): Cbar matrix.
!-

#include "CESR_platform.inc"

subroutine c_to_cbar (ele, cbar_mat)

  use bmad_struct
  
  implicit none

  type (ele_struct)  ele

  real(rp) cbar_mat(2,2), g_a(2,2), g_b_inv(2,2)
  real(rp) sqrt_beta_a, sqrt_beta_b, alpha_a, alpha_b

!

  sqrt_beta_a  = sqrt(ele%x%beta)
  sqrt_beta_b  = sqrt(ele%y%beta)
  alpha_a = ele%x%alpha
  alpha_b = ele%y%alpha

  g_a(1,1) = 1 / sqrt_beta_a
  g_a(1,2) = 0
  g_a(2,1) = alpha_a / sqrt_beta_a
  g_a(2,2) = sqrt_beta_a

  g_b_inv(1,1) = sqrt_beta_b
  g_b_inv(1,2) = 0
  g_b_inv(2,1) = -alpha_b / sqrt_beta_b
  g_b_inv(2,2) = 1 / sqrt_beta_b

  cbar_mat = matmul (matmul (g_a, ele%c_mat), g_b_inv)

end subroutine

