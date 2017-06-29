!+
! Subroutine err_exit
!
! Subroutine to first show the stack call list before exiting.
! This routine is typically used when a program detects an error condition.
!-

subroutine err_exit

implicit none

integer i, ix(1)
logical op

!

do i = 30, 60
  inquire (i, opened = op)
  if (op) close(op)
enddo

!

print *
print *, '!-------------------------------------------------'
print *, 'ERR_EXIT: WILL BOMB PROGRAM TO GIVE TRACEBACK...'
i = 0; ix = 0
print *, 1/i, ix(i)

stop

end subroutine
