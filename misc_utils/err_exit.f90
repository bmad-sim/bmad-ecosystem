!+
! Subroutine err_exit
!
! Subroutine to first show the stack call list before exiting.
! This routine is typically used when a program detects an error condition.
!-

#include "CESR_platform.inc"

subroutine err_exit

  use precision_def

  implicit none

  integer i, ix(1)

#if defined (CESR_VMS)
  call lib$signal (%val(0))
#else
  print *
  print *, '!-------------------------------------------------'
  print *, 'ERR_EXIT: WILL BOMB PROGRAM TO GIVE TRACEBACK...'
  i = 0; ix = 0
  print *, 1/i, ix(i)
#endif

  stop

end subroutine
