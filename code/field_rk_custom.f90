!+
! Subroutine field_rk_custom
!
! Dummy routine. Will generate an error if called.
!-

#include "CESR_platform.inc"

subroutine field_rk_custom

  print *, 'ERROR IN FIELD_RK_CUSTOM: THIS DUMMY ROUTINE SHOULD NOT HAVE'
  print *, '      BEEN CALLED IN THE FIRST PLACE.'
  call err_exit

end subroutine
