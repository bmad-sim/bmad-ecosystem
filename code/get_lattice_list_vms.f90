!+
! Subroutine get_lattice_list_vms (lat_list, num_lats, directory)
!
! Subroutine to get the names of the lattices of the form:
!     directory // BMAD_<lattice_name>.LAT
! or if there is no .LAT extension then a lattice file cam be of the form:
!     directory // BMAD_<lattice_name>
!
! Input:
!   directory  -- Character*(*): Directory to use. E.g: "U:[CESR.BMAD.LAT]"
!
! Output:
!   lat_list(*) -- Character*40: List of lattice names.
!   num_lats    -- Integer: Number of lattices found.
!-

#include "CESR_platform.inc"

subroutine get_lattice_list_vms (lat_list, num_lats, directory)

! keep compiler happy

  use cesr_utils

  implicit none

  character*(*) directory
  character*40 lat_list(*)
  character*200 lat_file

  integer num_lats, ios, context, ix, ixx, lib$find_file
  integer i, stat

  character*40 match_file

#ifdef CESR_VMS

  include '($ssdef)'
  include '($rmsdef)'

! get twiss file names for matching files 

  context = 0
  match_file = trim(directory) // 'BMAD_*.*'

  do i = 1, 10000

    stat = lib$find_file (match_file, lat_file, context, , , ios, 0)
    call str_upcase (lat_file, lat_file)

    if (stat) then
      ix = index(lat_file, ']') + 1       ! strip [...] prefix
      ixx = index(lat_file, ';') - 1      ! strip version number suffix
      lat_list(i) = lat_file(ix:ixx)
    else if (stat == rms$_nmf .or. stat == rms$_fnf) then
      num_lats = i - 1
      return
    else
      print *, 'FIND FILE ERROR:', stat
      call lib$signal(%val(stat))
      call err_exit
    endif

  enddo

  print *, 'GET_LATTICE_LIST: INTERNAL ERROR!'
  call err_exit

#endif

end subroutine
