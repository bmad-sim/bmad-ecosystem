!+
! Subroutine track1_runge_kutta
!
! Dummy routine for custom elements. Will generate an error if called.
!-

subroutine track1_runge_kutta

  print *, 'ERROR IN TRACK1_RUNGE_KUTTA: THIS DUMMY ROUTINE SHOULD NOT HAVE'
  print *, '      BEEN CALLED IN THE FIRST PLACE.'
  call err_exit

end subroutine
