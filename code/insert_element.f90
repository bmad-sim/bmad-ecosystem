!+
! Subroutine Insert_Element (lat, insert_ele, insert_index, ix_branch)
!
! Insert a new element into the tracking part of the lattice list.
! The existing elements from insert_index to n_ele_max get shoved up 1
! in the element array.
!
! Note: If the insert element has nonzero length then the S positions
! of all the elements will be recalculated.
!
! Modules Needed:
!   use bmad
!
! Input:
!     lat          -- lat_struct: lattice that will be modified
!     insert_ele   -- ele_struct: element to insert into the lat
!     insert_index -- integer: lat index where the new element is inserted.
!     ix_branch    -- Integer, optional :: branch index for the insertion. Default = 0.
!
! Output:
!     lat -- lat_struct: lattice with new element inserted
!-

subroutine insert_element (lat, insert_ele, insert_index, ix_branch)

use bmad_struct
use bmad_interface, except_dummy => insert_element

implicit none

type (lat_struct), target :: lat
type (ele_struct)  insert_ele
type (branch_struct), pointer :: branch, branch2
type (control_struct), pointer :: con

integer insert_index, ix, ix_br
integer, optional :: ix_branch

! transfer_ele is fast since it re reuse storage.

ix_br = integer_option(0, ix_branch)
branch => lat%branch(ix_br)

branch%n_ele_max = branch%n_ele_max + 1
if (branch%n_ele_max > ubound(branch%ele, 1)) call allocate_lat_ele_array(lat, ix_branch = ix_br)

do ix = branch%n_ele_max-1, insert_index, -1
  call transfer_ele (branch%ele(ix), branch%ele(ix+1))
  branch%ele(ix+1)%ix_ele = ix+1
enddo

! branch%ele(insert_index) pointers need to be nullified since they now point to
! the same memory as branch%ele(insert_index+1)

call deallocate_ele_pointers (branch%ele(insert_index), nullify_only = .true.)
branch%ele(insert_index) = insert_ele
branch%ele(insert_index)%ix_ele    = insert_index
branch%ele(insert_index)%ix_branch = ix_br

! Correct the control info

do ix = 1, lat%n_control_max
  con => lat%control(ix)
  if (con%ix_slave >= insert_index .and. con%ix_branch == ix_br)  con%ix_slave = con%ix_slave + 1
  if (con%ix_lord >= insert_index .and. ix_br == 0) con%ix_lord = con%ix_lord + 1
enddo

if (insert_index <= branch%n_ele_track + 1) then
  branch%n_ele_track = branch%n_ele_track + 1
  branch%n_ele_track = branch%n_ele_track
else
  print *, 'WARNING FROM INSERT_ELEMENT: YOU ARE INSERTING AN ELEMENT'
  print *, '        *NOT* INTO THE TRACKING PART OF THE LATTICE!'
  print *, '        ELEMENT: ', insert_ele%name
endif

do ix = 1, ubound(lat%branch, 1)
  branch2 => lat%branch(ix)
  if (branch2%ix_from_ele >= insert_index .and. branch2%ix_from_branch == ix_br) &
                                        branch2%ix_from_ele = branch2%ix_from_ele + 1
enddo

if (insert_ele%value(l$) /= 0) call s_calc(lat)

end subroutine
