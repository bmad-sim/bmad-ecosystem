!+
! Function lunget()
!
! Function to return a free file unit number to be used with an open statement.
!
! Modules needed:
!   use sim_utils
!
! Output:
!   lunget -- Integer: Free file number between 30 and 60
!
! Example:
!   iu = lunget()
!   open (iu, file = ... )
!-

integer function lunget()

  implicit none

  logical op
  integer lun
  lunget=0                ! happy compiler
  do lun=30,60
    inquire(lun,opened=op)
    if(.not.op) then
      lunget=lun
      return
    endif
  enddo
  print *, ' LUNGET: no logical unit openable between 30 and 60 ! '
  stop 

end function



