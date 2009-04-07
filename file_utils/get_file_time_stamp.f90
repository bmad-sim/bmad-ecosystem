!+
! Subroutine get_file_time_stamp (file, time_stamp)
!
! Routine to get the "last modified" time stamp for a file.
!
! Modules needed:
!   use cesr_utils
!
! Input:
!   file        -- Character(*): File name.
!
! Output:
!   time_stamp  -- Character(23): Modification time: "dd-mmm-yyyy hh:mm:ss.xx"
!-

#include "CESR_platform.inc"

subroutine get_file_time_stamp (file, time_stamp)

use cesr_interface, except => get_file_time_stamp

implicit none

character(*) file, time_stamp
character(60) line
integer iu, ix, ios

!

time_stamp = ''  ! In case of error

#if defined (CESR_VMS)

call lib$spawn ('directory/full/out=file_date.temp ' // file)

iu = lunget()
open (iu, file = 'file_date.temp', status = 'old')

do
  read (iu, '(a)', iostat = ios) line
  if (ios /= 0) exit
  if (line(1:8) == 'Revised:') then
    call string_trim (line(10:), line, ix)
    time_stamp = line(1:23)
    exit
  endif
enddo
close (iu)

call lib$spawn ('delete file_date.temp;*')

#else
  print *, 'ERROR: GET_FILE_TIME_STAMP NOT IMPLEMENTED FOR UNIX.'
  call err_exit
#endif

end subroutine
