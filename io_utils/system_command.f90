#include "CESR_platform.inc"
                                    
!+
! Subroutine system_command (line)
!
! Routine to execute an operating system command from within the program.
!
! Modules needed:
!   use cesr_utils
!
! Input:
!   line -- Character(*): Command to pass to the system shell.
!-

subroutine system_command (line)
implicit none
character(*) line
!

#if defined (CESR_VMS)
  call lib$spawn (line)
#else
  call system (line)
#endif

end subroutine
