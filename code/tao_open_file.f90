!+
! Subroutine tao_open_file (logical_dir, file_name, iunit, full_file_name)
!
! Subroutine to open a file for reading.
! This subroutine will first look for a file in the current directory before
! it looks in the logical_dir directory.
!
! Input:
!   logical_dir -- Character(*): Logical directory.
!   file_name   -- Character(*): File name.
!
! Output:
!   iunit          -- Integer: Logical unit number. Set to 0 if file not openable.
!   full_file_name -- Character(*): File name of found file.
!-

subroutine tao_open_file (logical_dir, file_name, iunit, full_file_name)

  use tao_mod

  implicit none

  character(*) logical_dir, file_name, full_file_name
  character(20) :: r_name = 'tao_open_file'

  integer iunit, ios
  logical valid

! open file

  iunit = lunget()
  full_file_name = file_name
  open (iunit, file = full_file_name, status = 'old', action = 'READ', &
                                                       iostat = ios)

! If we cannot open a file then try the the logical_dir 

  if (ios /= 0) then

    call fullfilename (trim(logical_dir) // ':' // file_name, full_file_name, valid)
    if (valid) then
      open (iunit, file = full_file_name, status = 'old', &
                                      action = 'READ', iostat = ios)
    endif

    ! If still nothing then this is an error.

    if (ios /= 0) then
      if (valid) then
         call out_io (s_blank$, r_name, 'File not found: ' // file_name, &
                                       '           Nor: ' // full_file_name)
      else
         call out_io (s_blank$, r_name, 'File not found: ' // file_name)
      endif
      iunit = 0
    endif

  endif

end subroutine
