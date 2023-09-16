!+
! Subroutine get_file_time_stamp (file, time_stamp)
!
! Routine to get the "last modified" time stamp for a file.
!
! Input:
!   file        -- Character(*): File name.
!
! Output:
!   time_stamp  -- Character(23): Modification time: "dd-mmm-yyyy hh:mm:ss.xx"
!-

subroutine get_file_time_stamp (file, time_stamp)

use sim_utils_interface, except => get_file_time_stamp

implicit none

character(*) file, time_stamp
character(60) line
integer iu, ix, ios

!

time_stamp = ''  ! In case of error

print *, 'ERROR: GET_FILE_TIME_STAMP NOT IMPLEMENTED FOR UNIX.'
if (global_com%exit_on_error) call err_exit

end subroutine
