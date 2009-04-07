!+
! Subroutine spawn_command (command)
!
! Subroutine to spawn a command to the Shell under unix or vms.
!
! Input:
!   command -- Character(*): The command to be spawned.
!-

#include "CESR_platform.inc"

subroutine spawn_command (command)

  implicit none

  character(*) command

!

#if defined (CESR_VMS)
  call lib$spawn (command)
#else
  call system (command)
#endif

end subroutine
