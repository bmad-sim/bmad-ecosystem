!+
! Subroutine type_this_file (FILENAME)
!
! Subroutine to type out a file to the screen.
!
! Example:
!     call type_this_file ('HELP_FILE.HELP')
!-

subroutine type_this_file(filename)

use filename_mod, except => type_this_file

implicit none

integer lun

character(*) filename
character(200) f_name, line

!

lun = lunget()
call fullfilename (filename, f_name)
open (unit = lun, file = f_name, status = 'old', action = 'READ', err = 9000)

do 
  read (lun, '(a)', end = 8000) line
  print '(1x, a)', trim(line)
enddo

8000  continue
close (unit = lun)
return

9000  print *, 'ERROR IN TYPE_THIS_FILE: CANNOT FIND: ', filename

end subroutine
 
