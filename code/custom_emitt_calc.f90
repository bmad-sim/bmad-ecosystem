!+
! Subroutine custom_emitt_calc
!
! Dummy routine for custom elements. Will generate an error if called.
!-

subroutine custom_emitt_calc

  print *, 'ERROR IN CUSTOM_EMITT_CALC: THIS DUMMY ROUTINE SHOULD NOT HAVE'
  print *, '      BEEN CALLED IN THE FIRST PLACE.'
  call err_exit

end subroutine
