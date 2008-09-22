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
!                 If Unix style, directory must end with a "/".
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

#if defined (CESR_VMS)

subroutine get_lattice_list (lat_list, num_lats, directory)

  use cesr_utils

  implicit none

  integer num_lats
  integer i, ix

  character(*) directory
  character(*) lat_list(:)
  character(200) directory2

  integer ios, context, ixx, lib$find_file
  integer stat

  character(40) match_file
  character(200) lat_file, l_file

!

  include '($ssdef)'
  include '($rmsdef)'

  call fullfilename (directory, directory2)
  match_file = trim(directory2) // '*.LAT'

! get twiss file names for matching files 

  context = 0

  i = 0
  do 

    stat = lib$find_file (match_file, lat_file, context, , , ios, 0)
    call str_upcase (lat_file, lat_file)

    if (stat) then
      ix = index(lat_file, ']') + 1       ! strip [...] prefix
      ixx = index(lat_file, ';') - 1      ! strip version number suffix
      l_file = lat_file(ix:ixx)
      if (l_file(1:8) == 'DIGESTED') cycle
      if (index(l_file, '.DIR') /= 0) cycle  ! Ignore layout subdirectory.
      if (i > size(lat_list)) then
        print *, 'ERROR IN GET_LATTICE_LIST: NUMBER OF LATTICES EXCEEDS ARRAY SIZE!'
        call err_exit
      endif
      ix = index(l_file, '.LAT')
      if (ix /= 0) l_file = l_file(1:ix-1)
      ix = index(l_file, '.lat')
      if (ix /= 0) l_file = l_file(1:ix-1)
      if (l_file(1:5) == 'BMAD_') l_file = l_file(6:)
      i = i + 1
      lat_list(i) = l_file
    else if (stat == rms$_nmf .or. stat == rms$_fnf) then
      num_lats = i - 1
      exit
    else
      print *, 'FIND FILE ERROR:', stat
      call lib$signal(%val(stat))
      call err_exit
    endif

  enddo

end subroutine get_lattice_list

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
! unix version

#else

subroutine get_lattice_list (lat_list, num_lats, directory)

  use dcslib
  use directory_mod

  implicit none

  integer num_lats
  integer i, ix

  character(*) directory
  character(*) lat_list(:)
  character(200) directory2
  character(40) match_file
  character(200) lat_file

  integer, automatic :: indx(size(lat_list))

  logical ok

  ! get twiss file names for matching files 

  call fullfilename (directory, directory2)
  match_file = '*.lat'

  ok = dir_open(directory2)

  num_lats = 0
  do 

    ok = dir_read (lat_file)
    if (.not. ok) exit

    if (.not. match_wild(lat_file, match_file)) cycle
    if (index(lat_file, ';') /= 0) cycle   ! ignore files with version numbers
    if (lat_file(1:8) == 'digested') cycle

    ! Old style: file = bmad_nnn.lat --> lattice name = nnn
    ! New style: file = nnn.lat      --> lattice name = nnn

    ix = index(lat_file, '.lat')   ! strip off .lat
    lat_file = lat_file(1:ix-1)

    ! Strip of beginning "bmad_" if it is there
    if (lat_file(1:5) /= 'bmad_') cycle
    lat_file = lat_file(6:)  

    num_lats = num_lats + 1
    if (num_lats > size(lat_list)) then
      print *, 'ERROR IN GET_LATTICE_LIST: NUMBER OF LATTICES EXCEEDS ARRAY SIZE!'
      call err_exit
    endif

     lat_list(num_lats) = lat_file  ! Strip of beginning "bmad_"

  enddo

  ! Sort in alphabetical order

  call indexx (lat_list(1:num_lats), indx(1:num_lats))
  lat_list(1:num_lats) = lat_list(indx(1:num_lats))

end subroutine get_lattice_list

#endif
