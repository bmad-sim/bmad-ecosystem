!+
! Subroutine Insert_Element (ring, insert_ele, insert_index)
!
! Insert a new element into the regular part of the ring structure.
! The existing elements from INSERT_INDEX to N_ELE_MAX get shoved up 1
! in the element array.
!
! Note: If the insert element has nonzero length then the S positions
! of all the elements will be recalculated.
!
! Modules Needed:
!   use bmad
!
! Input:
!     ring         -- ring_struct: ring that will be modified
!     insert_ele   -- ele_struct: element to insert into the ring
!     insert_index -- integer: ring index where the new element is inserted.
!
! Output:
!     ring -- ring_struct: ring with new element inserted
!
!-

!$Id$
!$Log$
!Revision 1.5  2003/05/02 15:44:00  dcs
!F90 standard conforming changes.
!
!Revision 1.4  2003/01/27 14:40:35  dcs
!bmad_version = 56
!
!Revision 1.3  2002/02/23 20:32:17  dcs
!Double/Single Real toggle added
!
!Revision 1.2  2001/09/27 18:31:52  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine insert_element (ring, insert_ele, insert_index)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ring_struct)  ring
  type (ele_struct)  insert_ele
  integer insert_index, index

!

  ring%n_ele_max = ring%n_ele_max + 1
  if (ring%n_ele_max > n_ele_maxx) then
    print *, 'ERROR IN INSERT_ELEMENT: NOT ENOUGH RING ELEMENTS!!!'
    print *, '      YOU NEED TO INCREASE N_ELE_MAXX IN BMAD_STRUCT!!!'
    call err_exit
  endif

  do index = ring%n_ele_max-1, insert_index, -1
    ring%ele_(index+1) = ring%ele_(index)
  enddo
  ring%ele_(insert_index) = insert_ele

  do index = 1, ring%n_control_array
    if (ring%control_(index)%ix_slave >= insert_index)  &
           ring%control_(index)%ix_slave = ring%control_(index)%ix_slave + 1
    if (ring%control_(index)%ix_lord >= insert_index)  &
           ring%control_(index)%ix_lord = ring%control_(index)%ix_lord + 1
  enddo

  if (insert_index <= ring%n_ele_ring+1) then
    ring%n_ele_ring = ring%n_ele_ring + 1
    ring%n_ele_use = ring%n_ele_use + 1
  else
    print *, 'WARNING FROM INSERT_ELEMENT: YOU ARE INSERTING AN ELEMENT'
    print *, '        *NOT* INTO THE REGULAR PART OF THE RING STRUCTURE!'
    print *, '        ELEMENT: ', insert_ele%name
  endif

  if (insert_index <= ring%n_ele_symm)  &
                                     ring%n_ele_symm = ring%n_ele_symm + 1

  if (insert_ele%value(l$) /= 0) call s_calc(ring)

end
