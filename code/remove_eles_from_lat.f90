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

use bmad_utils_mod

implicit none
                         
type (lat_struct), target :: lat
type (ele_struct), pointer :: ele

integer i, j, ix, i2
integer, allocatable :: ixa(:), ic(:), control(:)

! Compess lat%ele(:) array and fill in ixa(:) array.

allocate (ixa(lat%n_ele_max))           ! ixa(old_ele_index) = new_ele_index
allocate (control(lat%n_control_max))
allocate (ic(lat%n_ic_max))

control = 0
ic = 0

i2 = 0
do i = 1, lat%n_ele_max
  ele => lat%ele(i)
  if (ele%key == -1) then
    ixa(i) = -1
    control(ele%ix1_slave:ele%ix2_slave) = -1
    ic(ele%ic1_lord:ele%ic2_lord) = -1
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

! Renumber lat%control()%ix_ele  

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

do i = 1, lat%n_ele_max
  ele => lat%ele(i)
  if (ele%ix1_slave > 0) ele%ix1_slave = control(ele%ix1_slave)
  if (ele%ix2_slave > 0) ele%ix2_slave = control(ele%ix2_slave)
enddo

! Renumber lat%ic() array

i2 = 0
do i = 1, lat%n_ic_max
  if (ic(i) == -1) cycle
  i2 = i2 + 1
  ic(i) = i2
  lat%ic(i2) = control(lat%ic(i))
enddo

lat%n_ic_max = i2

do i = 1, lat%n_ele_max
  ele => lat%ele(i)
  if (ele%ic1_lord > 0) ele%ic1_lord = ic(ele%ic1_lord)
  if (ele%ic2_lord > 0) ele%ic2_lord = ic(ele%ic2_lord)
enddo

! deallocate and do a check

deallocate (ixa, control, ic)

call check_lat_controls (lat, .true.)

end subroutine
          
