!+
! Subroutine twiss_to_1_turn_mat (twiss, phi, mat2)
!
! Subroutine to form the 2x2 1-turn transfer matrix from the twiss parameters.
!
! Modules needed:
!   use bmad
!
! Input:
!   twiss -- Twiss_struct: Structure holding the Twiss parameters.
!     %beta
!     %alpha
!     %gamma
!   phi   -- Real(rdef): Tune in radians.
!
! Output:
!   mat2(2,2) -- Real(rdef): 1-turn matrix.
!-

!$Id$
!$Log$
!Revision 1.3  2002/02/23 20:32:29  dcs
!Double/Single Real toggle added
!
!Revision 1.2  2001/09/27 18:32:00  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine twiss_to_1_turn_mat (twiss, phi, mat2)

  use bmad

  type (twiss_struct) twiss

  real(rdef) phi, mat2(2,2), c, s

!

  c = cos(phi)
  s = sin(phi)

  mat2(1,1) =  c + s * twiss%alpha
  mat2(1,2) =  s * twiss%beta
  mat2(2,1) = -s * twiss%gamma
  mat2(2,2) =  c - s * twiss%alpha

end subroutine
