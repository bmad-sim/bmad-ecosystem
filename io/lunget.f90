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
  integer lun, i
  lunget=0                ! happy compiler
  do i=30,60
    lun = i
    inquire(lun,opened=op)
    if(.not.op) then
      lunget=i
      return
    endif
  enddo
  print *, ' LUNGET: no logical unit openable between 30 and 60 ! '
  stop 

end function



