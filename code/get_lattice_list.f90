!+
! Subroutine get_lattice_list (lat_list, num_lats, directory)
!
! Subroutine to choose between VMS and UNIX versions
!    
!
! Input:
!   directory  -- Character*(*): Directory to use. E.g: "U:[CESR.BMAD.LAT]"
!
! Output:
!   lat_list(*) -- Character*40: List of lattice names.
!   num_lats    -- Integer: Number of lattices found.
!-

#include "CESR_platform.inc"

subroutine get_lattice_list (lat_list, num_lats, directory)

  implicit none

  integer num_lats
  integer i, ix

  character*(*) directory
  character*200 directory2
  character*40 lat_list(:)

!

  call fullfilename (directory, directory2)

#ifdef CESR_VMS
  call get_lattice_list_vms(lat_list, num_lats, directory2)
#endif

#ifdef CESR_UNIX
  call get_lattice_list_unix(lat_list, num_lats, directory2)
#endif
  
! Strip off beginning "bmad_" and endding ".lat"

  do i = 1, num_lats
    lat_list(i) = lat_list(i)(6:)
    ix = index(lat_list(i), '.lat')
    if (ix /= 0) lat_list(i) = lat_list(i)(1:ix-1)
    ix = index(lat_list(i), '.LAT')
    if (ix /= 0) lat_list(i) = lat_list(i)(1:ix-1)
  enddo

end subroutine get_lattice_list
