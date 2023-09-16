module directory_mod

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function dir_list (dir_name, file_list) result (ok)
!
! Routine to return a list of files/subdirectories within a given directory.
!
! Input:
!   dir_name     -- character(*): Directory name.
!
! Output:
!   file_list(:) -- character(*), allocatable: A list of files/subdirectories.
!                     Files "." and ".." will not appear in the list.
!                     File names will not contain the path name.
!   ok           -- logical: True if everything is OK.
!-

function dir_list (dir_name, file_list) result (ok)

use re_allocate_mod

implicit none

character(*) dir_name
character(*), allocatable :: file_list(:)
character(200) file_name
integer i, n
logical ok

!

ok = dir_open(dir_name); 
if (.not. ok) then
  if (allocated(file_list)) deallocate(file_list)
  return
endif

n = 0
do
  if (.not. dir_read(file_name)) exit
  n = n + 1
enddo

call re_allocate (file_list, n)

ok = dir_open(dir_name); 
do i = 1, n
  ok = dir_read(file_list(i))
enddo
call dir_close()

end function dir_list

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function dir_open (dir_name) result (opened)
!
! Routine to open a directory to obtain a list of its files.
! Use this routine with dir_read and dir_close.
!
! Note: Subdirectories will not be searched.
!
! Note: calling this routine when there is already an open directory
!   will close the existing open directory.
!
! Input:
!   dir_name -- character(*): Directory name.
!
! Output:
!   opened -- logical: True if successful. False otherwise.
!-

function dir_open (dir_name) result (opened)

use fortran_cpp_utils

implicit none

character(*) dir_name
integer ok
logical opened

!

call open_dir (c_string(dir_name), ok)
opened = f_logic(ok)

end function dir_open

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function dir_read (file_name) result (valid)
!
! Routine to get the names of the files in a directory.
! Use this routine with dir_open and dir_close.
!
! After calling dir_open, call this routine repeatedly until
! valid = False. Each call will return the next file name in
! the list of files in the directory.
!
! Note: When there are no more files on the list and valid is
!   set False, dir_close will be automatically called.
!
! Note: File names will *not* contain the directory path.
!
! Note: Files "." and ".." will *not* appear in the list. 
!
! Output:
!   file_name -- Character(*): Name of a file.
!   valid     -- Logical: True if file_name is a valid file name.
!-

function dir_read (file_name) result (valid)

implicit none

character(*) file_name
logical valid

!

do
  call read_dir(file_name, len(file_name), valid)
  if (.not. valid) return
  if (file_name == ".") cycle
  if (file_name == "..") cycle
  return
enddo

end function dir_read

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine dir_close () 
!
! Routine to close a directory that was opened with dir_open.
! Also see dir_read.
! This routine only needs to be called. when the directory file
! list is not fully read.
!

subroutine dir_close () 

implicit none


!

call close_dir()

end subroutine dir_close

end module
