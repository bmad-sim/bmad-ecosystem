!+
! Subroutine MAKE_G2_MATS (TWISS, G_MAT, G_INV_MAT)
!
! Subroutine make the matrices needed to go from normal mode coords 
! to coordinates with the beta function removed.
!
! Modules Needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!     TWISS        -- Twiss_struct: Twiss parameters
!
! Output:
!     G_MAT(2,2)     -- Real: Normal mode to betaless coords
!     G_INV_MAT(2,2) -- Real: The inverse of G_MAT
!-

subroutine make_g2_mats (twiss, g2_mat, g2_inv_mat)

  use bmad_struct
  use bmad_interface

  implicit none

  type (twiss_struct) twiss

  real g2_mat(2,2), g2_inv_mat(2,2)
  real sqrt_beta, alpha
!

  sqrt_beta = sqrt(twiss%beta)
  alpha     = twiss%alpha

  g2_mat(1,1) = 1 / sqrt_beta
  g2_mat(1,2) = 0
  g2_mat(2,1) = alpha / sqrt_beta
  g2_mat(2,2) = sqrt_beta
                                  
  g2_inv_mat = g2_mat
  g2_inv_mat(2,1) = -g2_mat(2,1)

end subroutine
