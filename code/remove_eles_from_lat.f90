!+
! Subroutine remove_eles_from_lat (lat)
!
! Subroutine to compress the ele(:), control(:), and ic(:) arrays to remove
! elements no longer used. Note: to mark an element for removal use:
!     lat%ele(i)%key = -1
!
! Modules Needed:
!   use bmad
!
! Input:
!     lat -- lat_struct: Lattice to compress.
!
! Output:
!     lat -- lat_struct: Compressed lattice.
!-

#include "CESR_platform.inc"

subroutine remove_eles_from_lat (lat)

use bmad_struct
use bmad_interface

implicit none
                         
type (lat_struct), target :: lat
type (ele_struct), pointer :: ele

integer i, j, ix, i1, i2
integer, allocatable :: ixa(:), ic(:), control(:)

! Allocate

allocate (ixa(lat%n_ele_max))           ! ixa(old_ele_index) = new_ele_index
allocate (control(lat%n_control_max))
allocate (ic(lat%n_ic_max))

control = 0
ic = 0

! Mark entries in control and ic arrays for deletion.

do i = 1, lat%n_control_max
  if (lat%ele(lat%control(i)%ix_slave)%key == -1) control(i) = -1
  if (lat%ele(lat%control(i)%ix_lord)%key == -1) control(i) = -1
enddo

do i = 1, lat%n_ic_max
  if (control(lat%ic(i)) == -1) ic(i) = -1
enddo

! Compress lat%ele(:) array and fill in ixa(:) array.

i2 = 0
do i = 1, lat%n_ele_max
  ele => lat%ele(i)
  if (ele%key == -1) then
    ixa(i) = -1
  else
    i2 = i2 + 1
    ixa(i) = i2
    if (i2 /= i) lat%ele(i2) = ele
  endif
  if (i == lat%n_ele_track) then
     lat%n_ele_track = i2
  endif
enddo

if (i2 == lat%n_ele_max) then
  deallocate (ixa, control, ic)
  return
endif

do i = i2+1, lat%n_ele_max
  call init_ele(lat%ele(i))
enddo

lat%n_ele_max = i2

! Compress lat%control() array and correct %ix_lord and %ix_slave pointers.

i2 = 0
do i = 1, lat%n_control_max
  if (control(i) == -1) cycle
  i2 = i2 + 1
  control(i) = i2
  lat%control(i2) = lat%control(i)
  lat%control(i2)%ix_lord = ixa(lat%control(i)%ix_lord)
  lat%control(i2)%ix_slave = ixa(lat%control(i)%ix_slave)
enddo

lat%n_control_max = i2

! Correct ele%ix1_slave, etc.

do i = 1, lat%n_ele_max
  ele => lat%ele(i)
  if (ele%ix1_slave < 1) cycle
  i1 = ele%ix1_slave; i2 = ele%ix2_slave
  ele%ix1_slave = 0  ! Assume no slaves
  ele%ix2_slave = -1
  ele%n_slave = 0
  do j = i1, i2
    if (control(j) /= -1 .and. ele%ix1_slave == 0) ele%ix1_slave = control(j)
    if (control(j) /= -1) ele%ix2_slave = control(j)
  enddo
  ele%n_slave = ele%ix2_slave - ele%ix1_slave + 1
enddo

! Compress lat%ic() array

i2 = 0
do i = 1, lat%n_ic_max
  if (ic(i) == -1) cycle
  i2 = i2 + 1
  ic(i) = i2
  lat%ic(i2) = control(lat%ic(i))
enddo

lat%n_ic_max = i2

! Correct ele%ic1_lord, etc.

do i = 1, lat%n_ele_max
  ele => lat%ele(i)
  if (ele%ic1_lord < 1) cycle

  i1 = ele%ic1_lord; i2 = ele%ic2_lord
  ele%ic1_lord = 0  ! Assume no lords
  ele%ic2_lord = -1
  ele%n_lord = 0
  do j = i1, i2
    if (ic(j) /= -1 .and. ele%ic1_lord == 0) ele%ic1_lord = ic(j)
    if (ic(j) /= -1) ele%ic2_lord = ic(j)
  enddo
  ele%n_lord = ele%ic2_lord - ele%ic1_lord + 1
  if (ele%control_type == super_slave$ .and. ele%n_lord == 0) ele%control_type = free$
enddo

! deallocate and do a check

deallocate (ixa, control, ic)

call check_lat_controls (lat, .true.)

end subroutine
          
