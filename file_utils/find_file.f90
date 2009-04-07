!+
! Subroutine find_file (file_in, found, file_out, dirs)
!
! Subroutine to find a file. 
! If file_in does not have a specific directory specification then 
! the directories given by dir will also be searched.
!
! Modules needed:
!   use sim_utils
!
! Input:
!   file_in    -- Character(*): Input file name
!   dirs(:)    -- Character(*): Array of directories to search.
!                   Include ending '/', ']', or ':'.
!   
! Output:
!   found      -- Logical: Set True if a file is found. False otherwise.
!   file_out   -- Character(*): Full name of a file.
!                   If file is not found then file_out = file_in 
!-

#include "CESR_platform.inc"

subroutine find_file (file_in, found, file_out, dirs)

  use sim_utils, except => find_file

  implicit none

  character(*) file_in, dirs(:), file_out
  character(200) file_name, file_name_in

  integer i

  logical found, err

! Try without appending any dir

  file_name_in = file_in ! in case file_out == file_in

  call fullfilename (file_name_in, file_name, found)
  if (found) then
    inquire (file = file_name, exist = found, name = file_out)
    if (found) return
  endif

!

  do i = 1, size(dirs)
    call append_subdirectory (dirs(i), file_name_in, file_name, err)
    if (err) cycle
    call fullfilename (file_name, file_name, found)
    if (.not. found) cycle
    inquire (file = file_name, exist = found, name = file_out)
    if (found) return
  enddo  

  file_out = file_in ! file has not been found
  found = .false.

end subroutine
