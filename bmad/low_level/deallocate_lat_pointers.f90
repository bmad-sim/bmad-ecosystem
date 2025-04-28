!+
! Subroutine deallocate_lat_pointers (lat)
!
! Subroutine to deallocate the pointers in a lat.
!
! Input:
!   lat -- lat_struct: Lat with pointers.
!
! Output:
!   lat -- lat_struct: Lat with deallocated pointers.
!-

subroutine deallocate_lat_pointers (lat)

use bmad_routine_interface, dummy => deallocate_lat_pointers

implicit none

type (lat_struct) lat
integer i

!

call deallocate_ele_pointers (lat%ele_init)

if (allocated(lat%control))    deallocate(lat%control)
if (allocated(lat%ic))         deallocate(lat%ic)
if (allocated(lat%custom))     deallocate(lat%custom)
if (allocated(lat%print_str))  deallocate(lat%print_str)

!

if (allocated (lat%branch)) then
  do i = 0, ubound(lat%branch, 1)
    call deallocate_ele_array_pointers (lat%branch(i)%ele)
    call unlink_wall3d (lat%branch(i)%wall3d)
  enddo
  deallocate (lat%branch)
endif

!

nullify(lat%n_ele_track)
nullify(lat%n_ele_max)

lat%nametable%n_min = 0
lat%nametable%n_max = -1

end subroutine deallocate_lat_pointers
