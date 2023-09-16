!+
! Subroutine cbar_to_c (cbar_mat, a, b, c_mat)
!
! Subroutine to compute the C coupling matrix from the Cbar matrix and 
! the Twiss parameters.
!
! Input:
!   a             -- Twiss_struct: a-mode Twiss parameters
!   b             -- Twiss_struct: b-mode Twiss parameters
!   cbar_mat(2,2) -- Real(rp): Cbar matrix.
!
! Output:
!   c_mat(2,2) -- Real(rp): C matrix.
!-

subroutine cbar_to_c (cbar_mat, a, b, c_mat)

use bmad_interface, except_dummy => cbar_to_c

implicit none

type (twiss_struct) a, b

real(rp) cbar_mat(2,2), g_a_inv(2,2), g_b(2,2), c_mat(2,2)
real(rp) sqrt_beta_a, sqrt_beta_b, alpha_a, alpha_b

!

sqrt_beta_a  = sqrt(a%beta)
sqrt_beta_b  = sqrt(b%beta)
alpha_a = a%alpha
alpha_b = b%alpha

g_b(1,1) = 1 / sqrt_beta_b
g_b(1,2) = 0
g_b(2,1) = alpha_b / sqrt_beta_b
g_b(2,2) = sqrt_beta_b

g_a_inv(1,1) = sqrt_beta_a
g_a_inv(1,2) = 0
g_a_inv(2,1) = -alpha_a / sqrt_beta_a
g_a_inv(2,2) = 1 / sqrt_beta_a

c_mat = matmul (matmul (g_a_inv, cbar_mat), g_b)

end subroutine

