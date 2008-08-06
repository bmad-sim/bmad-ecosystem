!+
! Program to convert a Bmad file to an XSIF file and a MAD file
!
! Usage:
!   bmad_to_leo bmad_file_name
!
! The bmad_file_name will have a '.bmad' appended to the name if there
! is no '.' in the original name.
!
! The XSIF file name will be the bmad_file_name with the '.bmad' suffix
! (or whatever suffix is there) replaced by a '.xsif' suffix.
!
! The MAD file name will be the bmad_file_name with the '.bmad' suffix
! (or whatever suffix is there) replaced by a '.mad' suffix.
!-

program bmad_to_mad_and_xsif

use bmad
use write_lat_file_mod

implicit none

type (lat_struct) lat

integer n_arg
character(80) file_name

!

n_arg = cesr_iargc()

if (n_arg == 0) then
  write (*, '(a)', advance = 'NO') 'Bmad file name: '
  read (*, '(a)') file_name

elseif (n_arg == 1) then
  call cesr_getarg (1, file_name)

else
  print *, 'Usage: bmad_to_mad_and_xsif bmad_file_name'
  stop
endif

!

call file_suffixer (file_name, file_name, 'bmad', .false.)
call bmad_parser (file_name, lat)

call file_suffixer (file_name, file_name, 'xsif', .true.)
call bmad_to_xsif (file_name, lat)

call file_suffixer (file_name, file_name, 'mad', .true.)
call bmad_to_mad (file_name, lat)

end program
