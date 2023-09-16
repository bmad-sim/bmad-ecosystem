!+
! Program lattice_cleaner
!
! Program to "cleanup" a lattice. That is, transform a lattice to remove
! unwanted elements, etc. 
!
! Calling syntax:
!   lattice_cleaner <input-bmad-lattice-file-name>
!
! This program will produce a cleaned lattice file.
! 
! To save development time, rather than relying on an external file to 
! specify the lattice transformation rules, 
! this program can be modified directly to customize the transformation.
!- 


program lattice_cleaner

use bmad

implicit none

type (lat_struct) lat
type (ele_struct), pointer :: ele

integer i, ie, ix, n_loc

character(200) lat_file

logical err

! Get input lattice file name

if (command_argument_count() == 0) then
  print *, 'Command line synrax:'
  print *, '  lattic_cleaner <input-bmad-lattice-file-name>'
  stop
endif

if (command_argument_count() > 1) then
  print *, 'Too much stuff on the command line! Stopping here.'
  stop
endif

call get_command_argument (1, lat_file)

! Read in the lattice

call bmad_parser (lat_file, lat)

! Transformations:
!   1) Remove all marker elements.
!   2) Remove all multipole elements that do not have any non-zero multipole values.

do i = 1, lat%n_ele_max
  ele => lat%ele(i)

  if (ele%key == marker$) ele%ix_ele = -1  ! -1 -> mark for delection

  if (ele%key == multipole$ .or. ele%key == ab_multipole$) then
    if (all(ele%a_pole == 0) .and. all(ele%b_pole == 0)) ele%ix_ele = -1
  endif
enddo

call remove_eles_from_lat(lat, .true.)

! And write the transformed lattice file.

lat_file = 'clean_' // lat_file
call write_bmad_lattice_file (lat_file, lat)

print *, 'Cleaned lattice file: ', trim(lat_file)

end program
