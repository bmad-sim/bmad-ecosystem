!+
! Subroutine GET_FILE_NUMBER (FILE_NAME, CNUM_IN, NUM_OUT, ERR_FLAG)
!
!
! Subroutine to:
!     A) If CNUM_IN is positive then set NUM_OUT = CNUM_IN
!     B) If CNUM_IN is not positive then
!       1) Open FILE_NAME and read FILE_NUMBER
!       2) Set NUM_OUT = FILE_NUMBER - CNUM_IN
!
! Input:
!     FILE_NAME -- Character: File name.
!     CNUM_IN   -- Character: Encoded input number.
!
! Output:
!     NUM_OUT  -- Integer: Num_out from FILE_NAME.
!     ERR_FLAG -- Logical: Set true if decode error with CNUM_IN
!
! Example:
!     with num_out in file = 34
!     call get_file_number (file_name, '-3', num_out)
! Gives:
!     num_out = 31
!-

subroutine get_file_number (file_name, cnum_in, num_out, err_flag)

  use precision_def

  implicit none

  character(*) file_name, cnum_in
  integer lun, lunget, ios, num_in, num_out
  logical err_flag

! read cnum_in

  err_flag = .false.

  read (cnum_in, *, iostat = ios) num_in
  if (ios /= 0) then
    print *, 'ERROR IN GET_FILE_NUMBER: ERROR DECODEING: ', cnum_in
    err_flag = .true.
    return
  endif

! if num_in is positive then num_out = num_in

  if (num_in > 0) then
    num_out = num_in
    return
  endif

! else num_out is relative to number in file

  lun = lunget()
  open (unit = lun, file = file_name, status = 'old', &
                                      action = 'read', iostat = ios)
  if (ios /= 0) then
    print *, 'ERROR IN GET_FILE_NUMBER: FILE DOES NOT EXIST: ', file_name
    stop
  endif

  read (lun, *, iostat = ios) num_out
  close (unit = lun)
  if (ios /= 0) then
    print *, 'ERROR IN GET_FILE_NUM_OUT: CANNOT READ NUM_OUT IN FILE: '
    print *, '     ', file_name
    stop
  endif

  num_out = num_out + num_in

end subroutine
