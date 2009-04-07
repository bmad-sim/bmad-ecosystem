!+
! Subroutine FILE_DIRECTORIZER (IN_FILE, OUT_FILE, DIRECTORY, ADD_SWITCH)
!
! routine to add/replace a directory specification to a file name.
!
! Input:
!     IN_FILE     - Character string: Input file name
!     DIRECTORY   - Character string: directory specification to
!                     be added (DO INCLUDE THE SQUARE BRACKETS '[ ]').
!     ADD_SWITCH  - Logical:
!                       if .true.  directory spec is always added
!                       if .false. directory spec is only added if FILE
!                                doesn't allready contain a directory spec
!
! Output:
!     OUT_FILE    - Character string: Leading blanks are deleted and
!                           the directory spec is added if file
!                           name does not already contain a directory spec.
!
! Example:
!   In program:
!     infile = '[dcs.garbage]thisfile.foo'
!     call file_directorizer (infile, outfile, '[abc.def]', .true.)
!   Result:
!     outfile = '[abc.def]thisfile.foo'
!
! Example:
!   In program:
!     infile = '[dcs.garbage]thisfile.foo'
!     call file_directorizer (infile, outfile, '[abc.def]', .false.)
!   Result:
!     outfile = '[dcs.garbage]thisfile.foo'
!
! Example:
!   In program:
!     infile = 'thisfile.foo'
!     call file_directorizer (infile, outfile, '[abc.def]', .false.)
!   Result:
!     outfile = '[abc.def]thisfile.foo'
!-

#include "CESR_platform.inc"

subroutine file_directorizer (in_file, out_file, directory, add_switch)

  use precision_def

  implicit none

  integer ilen, dir_len, i, dix

  logical add_switch

  character(*) in_file, out_file, directory
  character(80) directory2

#if defined (CESR_VMS)
  character(2) :: dir_str = ']:'
#endif

#if defined(CESR_UNIX) || defined(CESR_WINCVF) || defined(CESR_LINUX)
  character(1) :: dir_str = '/'
#endif

! DIRECTORY2 is used so that DIRECTORY is not changed by program.
! Trim directory of leading blanks.

  directory2 = directory
  call string_trim (directory2, directory2, dir_len)

  out_file = in_file
  call string_trim (out_file, out_file, ilen) ! trim leading blanks
  out_file = out_file(1:ilen)      ! trim trailing words

  dix = 0
  do i = ilen, 1, -1
    if (index(dir_str, out_file(i:i)) /= 0) then
      dix = i
      exit
    endif
  enddo

! if IN_FILE doesn't already contain a directory then add DIRECTORY independent
! of the value of ADD_SWITCH. If IN_FILE contains a directory then only add a
! directory if ADD_SWITCH is .true.

  if (.not. add_switch .and. dix /= 0) return

  if (dir_len == 0) then
    out_file = out_file(dix+1:)
  else
    out_file = directory2(:dir_len) // out_file(dix+1:)  ! add directory
  endif

end subroutine

