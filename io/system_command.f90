!+
! Subroutine system_command (line, err_flag)
!
! Routine to execute an operating system command from within the program.
!
! Input:
!   line -- Character(*): Command to pass to the system shell.
!
! Output:
!   err_flag -- logical, optional: Set True if there is an error (bad command, etc.).
!-

subroutine system_command (line, err_flag)

implicit none
character(*) line
integer ix1, ix2
logical, optional :: err_flag

! Note: Return values of exitstat and cmdstat if there is an error are not standardized.

call execute_command_line (line, exitstat = ix1, cmdstat=ix2)
if (present(err_flag)) err_flag = (ix1 /= 0 .or. ix2 /= 0)

end subroutine
