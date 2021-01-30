                                    
!+
! Subroutine system_command (line)
!
! Routine to execute an operating system command from within the program.
!
! Modules needed:
!   use sim_utils
!
! Input:
!   line -- Character(*): Command to pass to the system shell.
!-

subroutine system_command (line)
implicit none
character(*) line
!

call execute_command_line (line)

end subroutine
