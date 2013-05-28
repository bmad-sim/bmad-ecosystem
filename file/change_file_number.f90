!+
! Subroutine CHANGE_FILE_NUMBER (FILE_NAME, CHANGE)
!
!
! Subroutine to:
!   1) Open FILE_NAME and read NUMBER
!   2) Change NUMBER to: NUMBER+CHANGE
!   3) Put changed NUMBER back into FILE_NAME
!
! Input:
!   FILE_NAME -- Character(*): File name.
!   CHANGE    -- Integer: Change in stored number.
!-

subroutine change_file_number (file_name, change)

  use precision_def

  implicit none

  character(*) file_name 
  integer lun, lunget, number, change, ios

! open file

  lun = lunget()
  open (unit = lun, file = file_name, status = 'old', iostat = ios)

  if (ios .ne. 0) then   ! error
    print *, 'ERROR IN CHANGE_FILE_NUMBER: CANNOT OPEN FILE: ', file_name
    stop
  endif

! read number

  read (lun, *, iostat = ios) number
  if (ios .ne. 0) then  ! error
    print *, 'ERROR IN CHANGE_FILE_NUMBER: CANNOT READ NUMBER IN FILE: '
    print *, '     ', file_name
    stop
  endif

! change number

  number = number + change
  if (number .lt. 0) number = 0
  rewind (unit = lun)
  write (lun, *) number
  close (unit = lun)

end subroutine
