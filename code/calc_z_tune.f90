!+
! Subroutine calc_z_tune (RING)
!
! Subroutine to calculate the synchrotron tune from the full 6X6 1 turn matrix
!
! Modules Needed:
!   use bmad
!
! Input:
!    RING  -- Ring_struct: Ring
!
! Output:
!    RING : ring%z%tune  synchrotron tune (radians)
!-

!$Id$
!$Log$
!Revision 1.3  2002/02/23 20:32:11  dcs
!Double/Single Real toggle added
!
!Revision 1.2  2001/09/27 18:31:48  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"



subroutine calc_z_tune ( ring)

  use bmad
  use nrtype
  use nr

  implicit none

  type (ring_struct) ring

  real(rdef) a(6,6), wr(6), wi(6), cos_z

  integer i
!

  call one_turn_matrix (ring, a)
  cos_z = (a(5,5) + a(6,6)) / (2 * (a(5,5)*a(6,6) - a(5,6)*a(6,5)))

  call balanc(a)
  call elmhes(a)
  call hqr(a,wr,wi)

! we need to find which eigen-value is closest to the z_tune

  i = minloc(abs(wr-cos_z), 1)
  ring%z%tune = -abs(atan2(wi(i),wr(i)))

end subroutine
