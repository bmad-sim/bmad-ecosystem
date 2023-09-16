!+
! Subroutine increment_file_number (file_name, digits, number, cnumber)
!
!
! Subroutine to:
!   1) Open FILE_NAME and read NUMBER
!   2) Increment NUMBER to: mod(NUMBER+1, 10**digits)
!   3) Put incremented NUMBER back into FILE_NAME
!
! Input:
!   FILE_NAME -- Character: File name.
!   DIGITS    -- Integer: Integer number of digits for NUM_STRING.
!
! Output:
!   NUMBER  -- Integer: Incremented number.
!   CNUMBER -- Character: Encoded number.
!
! Example:
!   with number in file_name = 34
!   call increment_file_number (file_name, 5, number, cnumber)
! Gives:
!   NUMBER = 35
!   CNUMBER(1:5) = '00035'
!-

subroutine increment_file_number (file_name, digits, number, cnumber)

  implicit none

  character(*) file_name, cnumber
  character(16) fmt
  integer lun, lunget, digits, number, ios

! init check

  if (digits < 1) then
    print *, 'ERROR IN INCREMENT_FILE_NUMBER: "DIGITS" MUST BE POSITIVE:', &
                                                                    digits
    stop
  endif

! open file

  lun = lunget()
  open (unit = lun, file = file_name, status = 'old', iostat = ios)

  if (ios /= 0) then   ! error
    print *, 'ERROR IN INCREMENT_FILE_NUMBER: CANNOT OPEN FILE: ', file_name
    stop
  endif

! read number

  read (lun, *, iostat = ios) number
  if (ios /= 0) then  ! error
    print *, 'ERROR IN INCREMENT_FILE_NUMBER: CANNOT READ NUMBER IN FILE: '
    print *, '     ', file_name
    stop
  endif

! increment number and encode

  number = number + 1
  if (number >= 10**digits) number = 1
  rewind (unit = lun)
  write (lun, *) number
  close (unit = lun)
  write (fmt, '(a, i2.2, a, i2.2, a)') '(i', digits, '.', digits, ')'
  write (cnumber, fmt) number

  return

  end
