program lattice_cleaner

use bmad
use write_lat_file_mod

implicit none

type (lat_struct) lat
type (ele_struct), pointer :: ele

integer i, ie, ix, n_loc

character(200) lat_file

logical err

!

call cesr_getarg (1, lat_file)


if (cesr_iargc() > 1) then
  print *, 'Too much stuff on the command line! Stopping here.'
  stop
endif

call bmad_parser (lat_file, lat)

! Mark elements for deletion.

do i = 1, lat%n_ele_max
  ele => lat%ele(i)

  if (ele%key == marker$) ele%key = -1  ! -1 -> mark for delection

  if (ele%key == multipole$ .or. ele%key == ab_multipole$) then
    if (all(ele%a_pole == 0) .and. all(ele%b_pole == 0)) ele%key = -1
  endif
enddo

! Delete 

call remove_eles_from_lat(lat, .true.)

! And write

lat_file = 'clean_' // lat_file
call write_bmad_lattice_file (lat_file, lat)

print *, 'Cleaned lattice file: ', trim(lat_file)

end program
