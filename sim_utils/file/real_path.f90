!+
! Function real_path (path_in, path_out) result (is_ok)
!
! Routine to expand all symbolic links and resolves references to "/./", "/../" and 
! extra '/' characters in path_in to produce a canonicalized absolute pathname.
!
! Note: If path_in does not represent an actual path, there can be problems.
!
! This is a wrapper for the C realpath routine.
!
! Input:
!   path_in     -- character(*): Input path
!
! Output:
!   path_out    -- character(*): Canonical path.
!   is_ok       -- logical: True if path converted without problems.

function real_path (path_in, path_out) result (is_ok)

use, intrinsic :: iso_c_binding
implicit none

integer, parameter :: PATH_MAX=1024 ! Would be safer to get from the system

character(*) path_in, path_out
logical is_ok

integer :: n
character(len=1) :: a(PATH_MAX)
type(c_ptr) :: ptr

! Fortran interface to C function, realpath()

interface
   function realpath(path,resolved_path) bind(c,name="realpath")
   use, intrinsic :: iso_c_binding
   type(c_ptr) :: realpath
   character(len=1,kind=c_char), intent(in) :: path(*)
   character(len=1,kind=c_char), intent(out) :: resolved_path(*)
   end function realpath
end interface

!

a = ''
ptr = realpath(trim(path_in) // C_NULL_CHAR, a)

is_ok = c_associated(ptr)

path_out = ''
do n = 1, PATH_MAX
  if (iachar(a(n)) == 0) exit
  path_out(n:n) = a(n)
end do

end function 
