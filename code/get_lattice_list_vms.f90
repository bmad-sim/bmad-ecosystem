!+
! Subroutine get_lattice_list (lat_list, num_lats, directory)
!
! Subroutine to get the names of the lattices of the form:
!     directory // BMAD_*.*
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
!Revision 1.2  2001/10/02 18:49:12  rwh24
!More compatibility updates; also added many explicit variable declarations.
!
!Revision 1.1  2001/09/27 18:33:14  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"

subroutine get_lattice_list (lat_list, num_lats, directory)

  character*(*) directory
  character*40 lat_list(*)
  integer num_lats

!keep compiler happy
#ifdef CESR_VMS
  use cesr_utils

  implicit none

  integer num_lats, ios, context, ix, ixx, lib$find_file
  integer i, stat

  character*40 match_file

#include '($ssdef)'
#include '($rmsdef)'

! get twiss file names for matching files 

  context = 0
  match_file = trim(directory) // 'BMAD_*.*'

  do i = 1, 10000

    stat = lib$find_file (match_file, lat_file, context, , , ios, 0)
    call str_upcase (lat_file, lat_file)

    if (stat) then
      ix = index(lat_file, ']BMAD') + 6   ! strip [~]BMAD_ prefix
      ixx = index(lat_file, ';') - 1      ! strip version number suffix
      lat_list(i) = lat_file(ix:ixx)
    else if (stat == rms$_nmf) then
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
  end



