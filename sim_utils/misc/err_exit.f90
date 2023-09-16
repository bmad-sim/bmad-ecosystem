!+
! Subroutine err_exit(err_str)
!
! Subroutine to first show the stack call list before exiting.
! This routine is typically used when a program detects an error condition.
!
! Input:
!   err_str   -- logical, optional: If present, print string before stopping.
!-

subroutine err_exit(err_str)

use output_mod, dummy => err_exit

implicit none

character(*), optional :: err_str
integer i, ix(1)
logical op
character(*), parameter :: r_name = 'err_exit'

! Close open files. The range [30, 60] is the range used by lunget.

do i = 30, 60
  inquire (i, opened = op)
  if (op) close(i)
enddo

if (present(err_str)) then
  call out_io (s_abort$, r_name, err_str)
endif

!

print *
print *, '!-------------------------------------------------'
print *, 'ERR_EXIT: WILL BOMB PROGRAM TO GIVE TRACEBACK...'
i = 0; ix = 0
print *, 1/i, ix(i)

stop

end subroutine
