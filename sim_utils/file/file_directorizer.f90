!+
! Subroutine file_directorizer (in_file, out_file, directory, add_switch)
!
! routine to add/replace a directory specification to a file name.
!
! Input:
!   in_file     - Character(*) string: Input file name
!   directory   - Character(*) string: directory specification to be added.
!   add_switch  - Logical:
!                       if .true.  directory spec is always added
!                       if .false. directory spec is only added if FILE
!                                doesn't allready contain a directory spec
!
! Output:
!   out_file    - Character(*) string: Leading blanks are deleted and
!                           the directory spec is added if file
!                           name does not already contain a directory spec.
!
! Example:
!   In program:
!     infile = '/dcs/garbage/thisfile.foo'
!     call file_directorizer (infile, outfile, 'abc/def', .true.)
!   Result:
!     outfile = 'abc/def/thisfile.foo'
!
! Example:
!   In program:
!     infile = '/dcs/garbage/thisfile.foo'
!     call file_directorizer (infile, outfile, 'abc/def', .false.)
!   Result:
!     outfile = '/dcs/garbage/thisfile.foo'
!
! Example:
!   In program:
!     infile = 'thisfile.foo'
!     call file_directorizer (infile, outfile, 'abc/def', .false.)
!   Result:
!     outfile = 'abc/def/thisfile.foo'
!-

subroutine file_directorizer (in_file, out_file, directory, add_switch)

use precision_def

implicit none

integer ilen, i, ix, dix

logical add_switch

character(*) in_file, out_file, directory

character(1) :: dir_str = '/'

! 

call string_trim (in_file, out_file, ilen) ! trim leading blanks
out_file = out_file(1:ilen)                 ! trim trailing words

dix = 0
do i = ilen, 1, -1
  if (index(dir_str, out_file(i:i)) /= 0) then
    dix = i
    exit
  endif
enddo

! if in_file doesn't already contain a directory then add directory independent
! of the value of add_switch. If in_file contains a directory then only add a
! directory if add_switch is .true.

if (.not. add_switch .and. dix /= 0) return

if (len_trim(directory) == 0) then
  out_file = out_file(dix+1:)
else
  ix = len_trim(directory)
  if (dir_str == '/' .and. directory(ix:ix) /= '/') then
    out_file = trim(directory) // '/' // out_file(dix+1:)  ! add directory
  else
    out_file = trim(directory) // out_file(dix+1:)  ! add directory
  endif
endif

end subroutine

