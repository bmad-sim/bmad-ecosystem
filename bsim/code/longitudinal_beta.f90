!........................................................................
!+
! Subroutine longitudinal_beta(ring, mode)
!
! Description:
!
! Arguments  :
!
! Mod/Commons:
!
! Calls      :
!
! Author     :
!
! Modified   :
!-
!........................................................................
!
! $Id$
!
! $Log$
! Revision 1.2  2007/01/30 16:14:31  dcs
! merged with branch_bmad_1.
!
! Revision 1.1.1.1.2.1  2006/12/22 20:30:42  dcs
! conversion compiles.
!
! Revision 1.1.1.1  2005/06/14 14:59:02  cesrulib
! Beam Simulation Code
!
!
!........................................................................
!
 subroutine longitudinal_beta(ring, mode)

  use bmad_interface

  implicit none

  type (lat_struct) ring
  type (normal_modes_struct) mode

  real(rp) z_emit
  real(rp) a(6,6)
  real(rp) cosmu, mu, sinmu, beta_z

! compute beta, alpha in longitudinal phase space

   z_emit = mode%sig_z * mode%sige_e 
   ring%ele(0)%z%beta = mode%sig_z**2/z_emit

! to get alpha as well as beta 

   call transfer_matrix_calc (ring, a)
!  subroutine one_turn_matrix was removed from BMAD early Oct, 2004
!  call one_turn_matrix(ring,.true.,a)

  cosmu = 0.5*(a(5,5)+a(6,6))
  mu = acos(cosmu)
  sinmu = sin(mu)
  beta_z = abs(a(5,6)/sinmu)
  ring%ele(0)%z%alpha = -0.5*(a(5,5)-a(6,6))
  ring%ele(0)%z%gamma = (1+ring%ele(0)%z%alpha**2)/beta_z

  if(abs(beta_z-ring%ele(0)%z%beta)/beta_z > 0.1)then
    print *,' longitudinal beta: beta_z = ',beta_z,'   delta_beta_z = ', &
            beta_z-ring%ele(0)%z%beta
  endif

  return

  end 
