!+
! Subroutine FILE_SUFFIXER (IN_FILE_NAME, OUT_FILE_NAME, SUFFIX, ADD_SWITCH)
!
! Routine to add/replace a suffix to a file name.
!
! Input:
!       IN_FILE_NAME  - Character string: Input file name
!       SUFFIX        - Character string: suffix to be added
!       ADD_SWITCH    - Logical: if .true. suffix is always added
!                                if .false. suffix is only added if FILE_NAME
!                                  doesn't allready contain a suffix
!
! Output:
!       OUT_FILE_NAME - Character string: Leading blanks are deleted and
!                       the suffix is added if file
!                       name does not already contain a suffix.
!
! Example:
!   In program:
!     infile = '[dcs.garbage]thisfile.foo'
!     call file_suffixer (infile, outfile, 'tum', .true.)
!   Result:
!     outfile = '[dcs.garbage]thisfile.tum'
!
! Example:
!   In program:
!     infile = '[dcs.garbage]thisfile.foo'
!     call file_suffixer (infile, outfile, 'tum', .false.)
!   Result:
!     outfile = '[dcs.garbage]thisfile.foo'
!
! Example:
!   In program:
!     infile = '[dcs.garbage]thisfile'
!     call file_suffixer (infile, outfile, 'tum', .false.)
!   Result:
!     outfile = '[dcs.garbage]thisfile.tum'
!
!-

#include "CESR_platform.inc"

subroutine file_suffixer (in_file_name, out_file_name, suffix, add_switch)

  use precision_def

  implicit none

  integer i, isl, ifl, pix

  logical add_switch

  character(*) in_file_name, out_file_name, suffix
  character(80) suffix2
#if defined (CESR_VMS)
  character(1) :: dir_str = ']'
#else
  character(1) :: dir_str = '/'
#endif

! SUFFIX2 is used so that SUFFIX is not changed by program

  suffix2 = suffix
  call string_trim (suffix2, suffix2, isl)  ! trim suffix of leading blanks
  suffix2 = suffix2(1:isl)                  ! trim trailing words
  if (suffix2(1:1) == '.') then           ! trim period if needed
    suffix2 = suffix2(2:)
  endif
  out_file_name = in_file_name
  call string_trim (out_file_name, out_file_name, ifl) ! trim leading blanks
  out_file_name = out_file_name(1:ifl)                 ! trim trailing words

! Locate period in file name if it exists. 

  pix = ifl+1
  do i = ifl, 1, -1
    if (out_file_name(i:i) == '.') then
      pix = i
      exit
    endif
    if (out_file_name(i:i) == dir_str) exit   ! no period found
  enddo

! if IN_FILE_NAME doesn't allready contain a suffix then add SUFFIX independent
! of the value of ADD_SWITCH. If IN_FILE_NAME contains a suffix then only add a
! suffix if ADD_SWITCH is .true.

  if (.not. add_switch .and. pix /= ifl+1) return

  if ((pix + isl) > len(out_file_name)) then
    print *, 'ERROR IN FILE_SUFFIXER: OUTPUT FILE NAME CHARACTER ARRAY'
    print *, '              FROM CALLING ROUTINE IS TOO SHORT!'
    call err_exit
  endif
  out_file_name(pix:) = '.' // suffix2   ! add suffix

end subroutine

