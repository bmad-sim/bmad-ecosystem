!+
! subroutine cbar_to_c (cbar_mat, ele)
!
! Subroutine to compute the C coupling matrix from the Cbar matrix and 
! the Twiss parameters.
!
! Modules Needed:
!   use bmad
!
! Input:
!   cbar_mat(2,2) -- Real(rp): Cbar matrix.
!   ele -- Ele_struct: Element with the Twiss parameters.
!     %x          -- a-mode Twiss parameters
!     %y          -- b-mode Twiss parameters
!
! Output:
!   ele -- Ele_struct: Element with the C matrix.
!     %c_mat(2,2) -- C matrix.
!-

#include "CESR_platform.inc"

subroutine cbar_to_c (cbar_mat, ele)

  use bmad_struct
  
  implicit none

  type (ele_struct)  ele

  real(rp) cbar_mat(2,2), g_a_inv(2,2), g_b(2,2)
  real(rp) sqrt_beta_a, sqrt_beta_b, alpha_a, alpha_b

!

  sqrt_beta_a  = sqrt(ele%x%beta)
  sqrt_beta_b  = sqrt(ele%y%beta)
  alpha_a = ele%x%alpha
  alpha_b = ele%y%alpha

  g_b(1,1) = 1 / sqrt_beta_b
  g_b(1,2) = 0
  g_b(2,1) = alpha_b / sqrt_beta_b
  g_b(2,2) = sqrt_beta_b

  g_a_inv(1,1) = sqrt_beta_a
  g_a_inv(1,2) = 0
  g_a_inv(2,1) = -alpha_a / sqrt_beta_a
  g_a_inv(2,2) = 1 / sqrt_beta_a

  ele%c_mat = matmul (matmul (g_a_inv, cbar_mat), g_b)

end subroutine

