!+
! Subroutine Insert_Element (lat, insert_ele, insert_index, ix_branch)
!
! Insert a new element into the tracking part of the lattice list.
! The existing elements from insert_index to n_ele_max get shoved up 1
! in the element array.
!
! Note: This routine is not for creating new control elements. For creating
!   new control elements use the routine new_control.
! Note: Bookkeeping like recalculating s-positions, reference energy, etc. is not done by this routine.
! Note: set_flags_for_changed_attribute is called for the inserted element.
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
use bookkeeper_mod

implicit none

type (lat_struct), target :: lat
type (ele_struct)  insert_ele
type (ele_struct), pointer :: inserted_ele, ele0
type (branch_struct), pointer :: branch, branch2
type (control_struct), pointer :: con

integer insert_index, ix, ix_br
integer, optional :: ix_branch

character(16), parameter :: r_name = 'insert_element'

logical err_flag

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

inserted_ele => branch%ele(insert_index)
call deallocate_ele_pointers (inserted_ele, nullify_only = .true.)
inserted_ele = insert_ele
inserted_ele%ix_ele    = insert_index
inserted_ele%ix_branch = ix_br
inserted_ele%lat => lat

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
  call out_io (s_warn$, r_name, &
                  'YOU ARE INSERTING AN ELEMENT *NOT* INTO THE TRACKING PART OF THE LATTICE!', &
                  'ELEMENT: ' // insert_ele%name)
endif

do ix = 1, ubound(lat%branch, 1)
  branch2 => lat%branch(ix)
  if (branch2%ix_from_ele >= insert_index .and. branch2%ix_from_branch == ix_br) &
                                        branch2%ix_from_ele = branch2%ix_from_ele + 1
enddo

call set_flags_for_changed_attribute (lat, inserted_ele)

! The reference energy bookkeeping will not be flagged for elements where the 
! reference energy is not changing. Do the bookkeeping locally.

ele0 => branch%ele(insert_index-1)
if (ele0%value(e_tot$) /= 0) then
  call attribute_bookkeeper (inserted_ele, branch%param)
  call ele_compute_ref_energy_and_time (inserted_ele, branch%param, &
                      ele0%value(e_tot$), ele0%value(p0c$), ele0%ref_time, err_flag)
endif

end subroutine
