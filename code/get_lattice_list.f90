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


subroutine get_lattice_list (lat_list, num_lats, directory)

  implicit none

  integer num_lats

  character*(*) directory
  character*40 lat_list(*)

#ifdef CESR_VMS
  call get_lattice_list_vms(lat_list, num_lats, directory)
#endif

#ifdef CESR_UNIX
  call get_lattice_list_unix(lat_list, num_lats, directory)
#endif
  
end subroutine get_lattice_list
