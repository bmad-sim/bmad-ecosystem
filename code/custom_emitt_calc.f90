!+
! Subroutine custom_emitt_calc
!
! Dummy routine for custom elements. Will generate an error if called.
!-

!$Id$
!$Log$
!Revision 1.2  2001/09/27 18:31:50  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine custom_emitt_calc

  print *, 'ERROR IN CUSTOM_EMITT_CALC: THIS DUMMY ROUTINE SHOULD NOT HAVE'
  print *, '      BEEN CALLED IN THE FIRST PLACE.'
  call err_exit

end subroutine
