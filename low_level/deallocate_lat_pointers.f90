!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
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

if (associated (lat%ele)) then
  call deallocate_ele_array_pointers (lat%ele)
  call deallocate_ele_pointers (lat%ele_init)
endif

if (allocated(lat%control))  deallocate (lat%control)
if (allocated(lat%ic))       deallocate (lat%ic)

! Do not need to deallocate stuff in lat%branch(0) since
! these pointers have been deallocated above.

if (allocated (lat%branch)) then
  call unlink_wall3d (lat%branch(0)%wall3d)

  do i = 1, ubound(lat%branch, 1)
    call deallocate_ele_array_pointers (lat%branch(i)%ele)
    deallocate (lat%branch(i)%param, lat%branch(i)%a, lat%branch(i)%b, lat%branch(i)%z)
    call unlink_wall3d (lat%branch(i)%wall3d)
  enddo
  deallocate (lat%branch)
endif

!

lat%n_ele_track  = -1
lat%n_ele_max  = -1

end subroutine deallocate_lat_pointers

