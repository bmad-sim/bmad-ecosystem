!+
! Subroutine custom_make_mat6
!
! Dummy routine for custom elements. Will generate an error if called.
!-

!$Id$
!$Log$
!Revision 1.2  2001/09/27 18:31:50  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine custom_make_mat6

  print *, 'ERROR IN CUSTOM_MAKE_MAT6: THIS DUMMY ROUTINE SHOULD NOT HAVE'
  print *, '      BEEN CALLED IN THE FIRST PLACE.'
  call err_exit

end subroutine
