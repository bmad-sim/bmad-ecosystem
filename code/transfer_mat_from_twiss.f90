!+
! subroutine transfer_mat_from_twiss (twiss1, twiss2, mat)
!
! Subroutine to make a 2 x 2 transfer matrix from the twiss parameters
! at the end points.
!
! Modules Needed:
!   use bmad
!
! Input:
!   twiss1  -- Twiss_struct: Twiss parameters at the initial point.
!   twiss2  -- Twiss_struct: Twiss parameters at the end point.
!
! Output:
!   mat(2,2) -- Real(rdef): Transfer matrix between the two points.
!-

!$Id$
!$Log$
!Revision 1.3  2002/02/23 20:32:27  dcs
!Double/Single Real toggle added
!
!Revision 1.2  2001/09/27 18:31:59  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine transfer_mat_from_twiss (twiss1, twiss2, mat)

  use bmad

  implicit none
      
  type (twiss_struct) twiss1, twiss2

  real(rdef) mat(2,2), a1, a2, b1, b2, sin21, cos21

!

  sin21 = sin(twiss2.phi - twiss1.phi)
  cos21 = cos(twiss2.phi - twiss1.phi)
  b1 = twiss1.beta;  b2 = twiss2.beta
  a1 = twiss1.alpha; a2 = twiss2.alpha

	mat(1,1) = sqrt(b2/b1) * (cos21 + a1 * sin21)
	mat(1,2) = sqrt(b1 * b2) * sin21
	mat(2,1) = -((a2-a1) * cos21 + (1 + a1 * a2) * sin21) / (sqrt(b1 * b2))
	mat(2,2) = sqrt(b1/b2) * (cos21 - a2 * sin21)

end subroutine

