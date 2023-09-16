MODULE touschek_background_mod

USE touschek_background_mpi_mod

TYPE orbit_data_struct
  TYPE(coord_struct), ALLOCATABLE :: orb(:)
END TYPE orbit_data_struct

CONTAINS



!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function raw_touschek_analysis(lat, slix, cur_ne, trackresults, orbit_data)
!
! Function to determine the index of the next tracking element. 
! To be used when tracking forwards or backwards.
!
! Module needed:
!   use track1_mod
!
! Input:
!   orb         -- Coord_struct: Starting coords.
!   lat         -- Lat_struct
!
! Output:
!   index       -- Integer: index of the next element
!
!-
subroutine raw_touschek_analysis(ele, cur_ne, trackresults, orbit_data)

type (lat_struct) :: lat
type (ele_struct) :: ele 
type (trackresult_struct) :: trackresults(:)
type (orbit_data_struct)  :: orbit_data(:)
type (coord_struct) :: end_orb

real(rp) :: cur_ne
real(rp) :: temp(6)

integer :: slix, particle_id

!

!do particle_id = 1, size(orbit_data)

  !end_orb = orbit_data(particle_id)%orb(trackresults%slixlost)
!  end_orb = orbit_data(1)%orb(1)
!  end_orb%p0c = ele%value(p0c$)

!end do


  

end subroutine raw_touschek_analysis

END MODULE touschek_background_mod
