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

#include "CESR_platform.inc"

subroutine insert_element (ring, insert_ele, insert_index)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ring_struct)  ring
  type (ele_struct)  insert_ele
  integer insert_index, index

! transfer_ele is fast since re reuse storage.

  ring%n_ele_max = ring%n_ele_max + 1
  if (ring%n_ele_max > ring%n_ele_maxx) call allocate_ring_ele_(ring)

  do index = ring%n_ele_max-1, insert_index, -1
    call transfer_ele (ring%ele_(index), ring%ele_(index+1))
  enddo
  
! ring%ele_(insert_index) pointers need to be nullified since they now point to
! the same memory as ring%ele_(insert_index+1)

  call deallocate_ele_pointers (ring%ele_(insert_index), nullify_only = .true.)
  ring%ele_(insert_index) = insert_ele

! correct the control info

  do index = 1, ring%n_control_max
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
