!+
! Subroutine file_suffixer (in_file_name, out_file_name, suffix, add_switch)
!
! Routine to add/replace a suffix to a file name.
!
! Input:
!   in_file_name  - Character string: Input file name
!   suffix        - Character string: suffix to be added
!   add_switch    - Logical: if .true. suffix is always added
!                            if .false. suffix is only added if FILE_NAME
!                              doesn't allready contain a suffix
!
! Output:
!   out_file_name - Character string: Leading blanks are deleted and
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

subroutine file_suffixer (in_file_name, out_file_name, suffix, add_switch)

use precision_def

implicit none

integer i, isl, ifl, pix

logical add_switch

character(*) in_file_name, out_file_name, suffix
character(80) suffix2
character(1) :: dir_str = '/'

! Suffix2 is used so that SUFFIX is not changed by program

suffix2 = suffix
call string_trim (suffix2, suffix2, isl)  ! trim suffix of leading blanks
suffix2 = suffix2(1:isl)                  ! trim trailing words
if (suffix2(1:1) == '.') then             ! trim period if needed
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

! if in_file_name doesn't allready contain a suffix then add suffix independent
! of the value of add_switch. If in_file_name contains a suffix then only add a
! suffix if add_switch is .true.

if (.not. add_switch .and. pix /= ifl+1) return

if ((pix + isl) > len(out_file_name)) then
  print *, 'ERROR IN FILE_SUFFIXER: OUTPUT FILE NAME CHARACTER ARRAY'
  print *, '              FROM CALLING ROUTINE IS TOO SHORT!'
  if (global_com%exit_on_error) call err_exit
endif
out_file_name(pix:) = '.' // suffix2   ! add suffix

end subroutine

