!+
! Subroutine milli_sleep (milli_sec)
!
! Routine to pause the program for a given number of milli-seconds.
!
! Input:
!   milli_sec -- Integer: Number of milli-seconds to pause.
!-

subroutine milli_sleep (milli_sec)

use precision_def

implicit none

integer milli_sec

!

call program_sleep (real(1d-3 * milli_sec, dp))

end subroutine
