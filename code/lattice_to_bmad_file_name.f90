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

  implicit none

  character(*) lattice, bmad_file_name

!

  if (index(lattice, '.') == 0) then
    bmad_file_name = 'BMAD_LAT:bmad_' // trim(lattice) // '.lat'
  else
    bmad_file_name = 'BMAD_LAT:bmad_' // trim(lattice) 
  endif

  call FullFileName(bmad_file_name, bmad_file_name)

end subroutine
