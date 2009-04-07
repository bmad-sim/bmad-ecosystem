!+
! Subroutine type_this_file (FILENAME)
!
! Subroutine to type out a file to the screen.
!
! Example:
!     call type_this_file ('HELP_FILE.HELP')
!-

#include "CESR_platform.inc"

subroutine type_this_file(filename)

  use precision_def

  implicit none

  integer lun, lunget

  character(*) filename
  character(132) line

!

  lun = lunget()
  open (unit = lun, file = filename, status = 'old', &
                                             action = 'READ', err = 9000)

  do 
    read (lun, '(a)', end = 8000) line
    print '(1x, a)', trim(line)
  enddo

8000  continue
  close (unit = lun)
  return

9000  print *, 'ERROR IN TYPE_THIS_FILE: CANNOT FIND: ', filename

end subroutine
 
