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

!$Id$
!$Log$
!Revision 1.9  2003/01/27 14:40:35  dcs
!bmad_version = 56
!
!Revision 1.8  2002/02/23 20:32:16  dcs
!Double/Single Real toggle added
!
!Revision 1.7  2002/01/11 21:16:00  dcs
!Bug fix
!
!Revision 1.6  2002/01/11 16:58:20  dcs
!Minor bug fix
!
!Revision 1.5  2002/01/11 16:32:13  cesrulib
!Fixed typos
!
!Revision 1.4  2001/10/08 17:18:14  rwh24
!DCS changes to f90 files.
!Bug fixes to c file.
!
!Revision 1.3  2001/10/05 18:23:57  rwh24
!Bug Fixes
!
!Revision 1.2  2001/10/02 18:49:12  rwh24
!More compatibility updates; also added many explicit variable declarations.
!
!Revision 1.1  2001/09/27 18:33:14  rwh24
!UNIX compatibility updates
!
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
      type *, 'FIND FILE ERROR:', stat
      call lib$signal(%val(stat))
      call err_exit
    endif

  enddo

  type *, 'GET_LATTICE_LIST: INTERNAL ERROR!'
  call err_exit

#endif

end subroutine








