!+
! Subroutine custom_make_mat6
!
! Dummy routine for custom elements. Will generate an error if called.
!-

subroutine custom_make_mat6

  print *, 'ERROR IN CUSTOM_MAKE_MAT6: THIS DUMMY ROUTINE SHOULD NOT HAVE'
  print *, '      BEEN CALLED IN THE FIRST PLACE.'
  call err_exit

end subroutine
