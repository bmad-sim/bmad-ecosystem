!+         
! Subroutine one_turn_matrix (RING, MAT6)
!
! Subroutine to calculate the full 6X6 1 turn matrix
!
! Modules Needed:
!   use bmad
!
! Input:
!     RING   -- Ring_struct: Ring
!
! Output:
!    MAT6(6,6) --  Matrix Type           ring.param.symmetry
!                  --------              -----------
!                  1-turn                none$
!-

!$Id$
!$Log$
!Revision 1.5  2003/01/27 14:40:41  dcs
!bmad_version = 56
!
!Revision 1.4  2002/02/23 20:32:22  dcs
!Double/Single Real toggle added
!
!Revision 1.3  2002/01/08 21:44:42  dcs
!Aligned with VMS version  -- DCS
!
!Revision 1.2  2001/09/27 18:31:55  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"



subroutine one_turn_matrix (ring, mat6)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ring_struct)  ring

  real(rdef) mat6(6,6)
  integer i

  mat6 = ring%ele_(1)%mat6
  do i=2,ring%n_ele_use                 
    mat6 = matmul (ring%ele_(i)%mat6, mat6)
  end do
    

end subroutine
                                          
