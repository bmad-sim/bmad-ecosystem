!+
! Subroutine get_lattice_list (lat_list, num_lats, directory)
!
! Subroutine to get the names of the lattices of the form:
!     directory // bmad_<lattice_name>.lat
! or if there is no .lat extension then a lattice file cam be of the form:
!     directory // bmad_<lattice_name>
!
! Input:
!   directory  -- Character(*): Directory to use. E.g: "U:[CESR.BMAD.LAT]"
!
! Output:
!   lat_list(*) -- Character(40): List of lattice names.
!   num_lats    -- Integer: Number of lattices found.
!-

#include "CESR_platform.inc"

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
! VMS version

#ifdef CESR_VMS

subroutine get_lattice_list (lat_list, num_lats, directory)

  use cesr_utils

  implicit none

  integer num_lats
  integer i, ix

  character(*) directory
  character(200) directory2
  character(40) lat_list(:)

  integer ios, context, ixx, lib$find_file
  integer stat

  character(40) match_file
  character(200) lat_file

!

  include '($ssdef)'
  include '($rmsdef)'

  call fullfilename (directory, directory2)
  match_file = trim(directory2) // 'BMAD_*.*'

! get twiss file names for matching files 

  context = 0

  i = 0
  do 

    i = i + 1
    stat = lib$find_file (match_file, lat_file, context, , , ios, 0)
    call str_upcase (lat_file, lat_file)

    if (stat) then
      if (i > size(lat_list)) then
        print *, 'ERROR IN GET_LATTICE_LIST_VMS: NUMBER OF LATTICES EXCEEDS ARRAY SIZE!'
        call err_exit
      endif
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

! Strip off beginning "bmad_" and endding ".lat"

  do i = 1, num_lats
    lat_list(i) = lat_list(i)(6:)
    ix = index(lat_list(i), '.lat')
    if (ix /= 0) lat_list(i) = lat_list(i)(1:ix-1)
    ix = index(lat_list(i), '.LAT')
    if (ix /= 0) lat_list(i) = lat_list(i)(1:ix-1)
  enddo

end subroutine get_lattice_list


!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
! unix version

#else

subroutine get_lattice_list (lat_list, num_lats, directory)

  use cesr_utils

  implicit none

  character(*) directory
  character(200) directory2
  character(40) lat_list(:)

  integer num_lats, ios, context, ix, ixx, lib$find_file
  integer i, stat

  character(40) match_file

!

  print *, 'ERROR: GET_LATTICE_LIST NOT YET IMPLEMENTED FOR UNIX!'
  call err_exit

end subroutine

#endif
