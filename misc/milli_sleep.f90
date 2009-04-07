!+
! Subroutine milli_sleep (milli_sec)
!
! Routine to pause the program for a given number of milli-seconds.
!
! Modules needed:
!   use sim_utils
!
! Input:
!   milli_sec -- Integer: Number of milli-seconds to pause.
!-

#include "CESR_platform.inc"

subroutine milli_sleep (milli_sec)

use precision_def

implicit none

integer milli_sec

!

#if defined (CESR_VMS)
call lib$wait(.001*float(milli_sec))
#else
call program_sleep (real(1d-3 * milli_sec, dp))
#endif

end subroutine
