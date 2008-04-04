!+
! Subroutine lattice_to_bmad_file_name (lattice, bmad_file_name)
!
! Subroutine to convert a lattice name to the appropriate bmad file name.
!
! Input:
!   lattice -- Character(40): CESR Lattice name.
!
! Output:
!  bmad_file_name -- Character(*): Appropriate name of the bmad file
!-

subroutine lattice_to_bmad_file_name (lattice, bmad_file_name)

  use cesr_utils

  implicit none

  character(*) lattice, bmad_file_name
  character(len(lattice)) lat
  logical is_there

!

  lat = lattice
  call downcase_string(lat)

  bmad_file_name = 'BMAD_LAT:' // trim(lat) // '.lat'
  call FullFileName(bmad_file_name, bmad_file_name)

  inquire (file = bmad_file_name, exist = is_there)
  if (.not. is_there) then
    bmad_file_name = 'BMAD_LAT:bmad_' // trim(lat) // '.lat'
    call FullFileName(bmad_file_name, bmad_file_name)
  endif
end subroutine
