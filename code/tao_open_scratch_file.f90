!+
! Function tao_open_scratch_file (err) result (iu)
!
! Routine to open a scratch file.
! To close the file use:
!   close (iu, status = 'delete')
!
! Output:
!   err       -- logical: Set True if there is an error. False otherwise.
!   iu        -- integer: File handle unit number.
!- 

Function tao_open_scratch_file (err) result (iu)

use sim_utils

implicit none

integer iu, ios
logical err, valid
character(100) file_name
character(*), parameter :: r_name = 'tao_open_scratch_file'


! Try to open a scratch file.

err = .false.
iu = lunget()
open (iu, status = 'scratch', iostat = ios)
if (ios == 0) return

! The alternative is to open a regular file in the home directory

call fullfilename('$HOME/scratch.file', file_name, valid)
if (.not. valid) then
  err = .true.
  call out_io (s_error$, r_name, 'CANNOT OPEN A SCRATCH FILE!')
  return
endif

open (iu, file = file_name, iostat = ios)
if (ios /= 0) then
  err = .true.
  call out_io (s_error$, r_name, 'CANNOT OPEN A SCRATCH FILE!')
  return
endif

end function
