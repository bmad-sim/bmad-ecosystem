!+
! Subroutine track1_runge_kutta
!
! Dummy routine for custom elements. Will generate an error if called.
!-

!$Id$
!$Log$
!Revision 1.2  2001/09/27 18:31:58  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine track1_runge_kutta

  print *, 'ERROR IN TRACK1_RUNGE_KUTTA: THIS DUMMY ROUTINE SHOULD NOT HAVE'
  print *, '      BEEN CALLED IN THE FIRST PLACE.'
  call err_exit

end subroutine
