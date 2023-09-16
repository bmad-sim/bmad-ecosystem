!+
! Subroutine calc_file_number (file_name, num_in, num_out, err_flag)
!
! Subroutine to:
!     A) If num_in is positive then set NUM_OUT = NUM_IN
!     B) If num_in is not positive then
!       1) Open file_name and read file_number
!       2) Set num_out = file_number - num_in
!
! Input:
!     file_name -- Character: File name.
!     num_in    -- Integer: Input number.
!
! Output:
!     num_out  -- Integer: Num_out from file_name.
!     err_flag -- Logical: Set true if decode error with num_in
!
! Example:
!     with num_out in file = 34
!     call calc_file_number (file_name, -3, num_out)
! Gives:
!     num_out = 31
!-

subroutine calc_file_number (file_name, num_in, num_out, err_flag)

  use sim_utils, except => calc_file_number

  implicit none


  character(*) file_name
  character(200) f_name
  integer lun, ios, num_in, num_out, num_in_file
  logical err_flag

! if num_in is positive then num_out = num_in

  err_flag = .false.

  if (num_in > 0) then
    num_out = num_in
    return
  endif

! else num_out is relative to number in file

  lun = lunget()
  call fullfilename (file_name, f_name)
  open (unit = lun, file = f_name, status = 'old',  &
                                       action = 'read', iostat = ios)
  if (ios /= 0) then
    print *, 'ERROR IN CALC_FILE_NUMBER: FILE DOES NOT EXIST: '
    print *, '      ', trim(f_name)
    if (global_com%exit_on_error) call err_exit
  endif

  read (lun, *, iostat = ios) num_in_file
  close (unit = lun)
  if (ios /= 0) then
    print *, 'ERROR IN CALC_FILE_NUM_OUT: CANNOT READ NUM_OUT IN FILE: '
    print *, '     ', trim(f_name)
    if (global_com%exit_on_error) call err_exit
  endif

  num_out = num_in_file + num_in

end subroutine
