!+
! Subroutine MAKE_G2_MATS (TWISS, G_MAT, G_INV_MAT)
!
! Subroutine make the matrices needed to go from normal mode coords 
! to coordinates with the beta function removed.
!
! Modules Needed:
!   use bmad
!
! Input:
!     TWISS        -- Twiss_struct: Twiss parameters
!
! Output:
!     G_MAT(2,2)     -- Real(rdef): Normal mode to betaless coords
!     G_INV_MAT(2,2) -- Real(rdef): The inverse of G_MAT
!-

!$Id$
!$Log$
!Revision 1.3  2002/02/23 20:32:17  dcs
!Double/Single Real toggle added
!
!Revision 1.2  2001/09/27 18:31:52  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine make_g2_mats (twiss, g2_mat, g2_inv_mat)

  use bmad

  implicit none

  type (twiss_struct) twiss

  real(rdef) g2_mat(2,2), g2_inv_mat(2,2)
  real(rdef) sqrt_beta, alpha
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
