!+
! Subroutine make_g2_mats (twiss, g_mat, g_inv_mat)
!
! Subroutine make the matrices needed to go from normal mode coords
! to coordinates with the beta function removed.
!
! Input:
!   twiss        -- Twiss_struct: Twiss parameters.
!
! Output:
!   g_mat(2,2)     -- Real(rp): Normal mode to betaless coords.
!   g_inv_mat(2,2) -- Real(rp): The inverse of g_mat.
!-

subroutine make_g2_mats (twiss, g2_mat, g2_inv_mat)

use bmad_struct

implicit none

type (twiss_struct) twiss

real(rp) g2_mat(2,2), g2_inv_mat(2,2)
real(rp) sqrt_beta, alpha
!

sqrt_beta = sqrt(twiss%beta)
alpha     = twiss%alpha

g2_mat(1,1) = 1 / sqrt_beta
g2_mat(1,2) = 0
g2_mat(2,1) = alpha / sqrt_beta
g2_mat(2,2) = sqrt_beta

g2_inv_mat = g2_mat
g2_inv_mat(2,1) = -g2_mat(2,1)

end subroutine make_g2_mats

