!+
! Subroutine custom_track1
!
! Dummy routine for custom elements. Will generate an error if called.
!-

subroutine custom_track1

  print *, 'ERROR IN CUSTOM_TRACK1: THIS DUMMY ROUTINE SHOULD NOT HAVE'
  print *, '      BEEN CALLED IN THE FIRST PLACE.'
  call err_exit

end subroutine
