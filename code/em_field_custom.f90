!+
! Subroutine em_field_custom
!
! Dummy routine. Will generate an error if called.
!-

#include "CESR_platform.inc"

subroutine em_field_custom

  print *, 'ERROR IN EM_FIELD_CUSTOM: THIS DUMMY ROUTINE SHOULD NOT HAVE'
  print *, '      BEEN CALLED IN THE FIRST PLACE.'
  call err_exit

end subroutine
