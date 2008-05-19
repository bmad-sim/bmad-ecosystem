!+
! Subroutine Insert_Element (lat, insert_ele, insert_index)
!
! Insert a new element into the tracking part of the lattice list.
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
!     lat         -- lat_struct: lat that will be modified
!     insert_ele   -- ele_struct: element to insert into the lat
!     insert_index -- integer: lat index where the new element is inserted.
!
! Output:
!     lat -- lat_struct: lat with new element inserted
!-

#include "CESR_platform.inc"

subroutine insert_element (lat, insert_ele, insert_index)

  use bmad_struct
  use bmad_interface, except_dummy => insert_element

  implicit none

  type (lat_struct), target :: lat
  type (ele_struct)  insert_ele
  type (branch_struct), pointer :: line
  integer insert_index, ix

! transfer_ele is fast since re reuse storage.

  lat%n_ele_max = lat%n_ele_max + 1
  if (lat%n_ele_max > ubound(lat%ele, 1)) call allocate_ele_array(lat%ele)

  do ix = lat%n_ele_max-1, insert_index, -1
    call transfer_ele (lat%ele(ix), lat%ele(ix+1))
    lat%ele(ix+1)%ix_ele = ix+1
  enddo
  
! lat%ele(insert_index) pointers need to be nullified since they now point to
! the same memory as lat%ele(insert_index+1)

  call deallocate_ele_pointers (lat%ele(insert_index), nullify_only = .true.)
  lat%ele(insert_index) = insert_ele
  lat%ele(insert_index)%ix_ele = insert_index

! correct the control info

  do ix = 1, lat%n_control_max
    if (lat%control(ix)%ix_slave >= insert_index)  &
           lat%control(ix)%ix_slave = lat%control(ix)%ix_slave + 1
    if (lat%control(ix)%ix_lord >= insert_index)  &
           lat%control(ix)%ix_lord = lat%control(ix)%ix_lord + 1
  enddo

  if (insert_index <= lat%n_ele_track + 1) then
    lat%n_ele_track  = lat%n_ele_track + 1
    lat%n_ele_track = lat%n_ele_track
  else
    print *, 'WARNING FROM INSERT_ELEMENT: YOU ARE INSERTING AN ELEMENT'
    print *, '        *NOT* INTO THE REGULAR PART OF THE LAT STRUCTURE!'
    print *, '        ELEMENT: ', insert_ele%name
  endif

  if (allocated(lat%branch)) then
    do ix = 1, size(lat%branch)
      line => lat%branch(ix)
      if (line%ix_from_ele >= insert_index .and. line%ix_from_line == 0) &
                                   line%ix_from_ele = line%ix_from_ele + 1
    enddo
  endif

  if (insert_ele%value(l$) /= 0) call s_calc(lat)

end subroutine
