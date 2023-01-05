!+
! Subroutine insert_element (lat, insert_ele, ix_ele, ix_branch, orbit)
!
! Insert a new element into the tracking part of the lattice list.
! The existing elements from ix_ele to n_ele_max get shoved up 1 in the element array.
!
! Additionally, the length of an optional orbit array can be correspondingly increased
! by addition of an element at the corresponding point in the orbit array.
!
! Notes:
! 1) This routine is not for creating new control elements. For creating
!   new control elements use the routine new_control.
! 2) Bookkeeping like recalculating s-positions, reference energy, etc. is not done by this routine.
! 3) set_flags_for_changed_attribute is called for the inserted element.
!
! Input:
!   lat          -- lat_struct: lattice that will be modified
!   insert_ele   -- ele_struct: element to insert into the lat
!   ix_ele       -- integer: branch%ele(:) index where the new element is inserted.
!   ix_branch    -- Integer, optional :: branch index for the insertion. Default = 0.
!   orbit(:)     -- coord_struct, optional, allocatable: orbit array to enlarge.
!
! Output:
!   lat       -- lat_struct: lattice with new element inserted
!   orbit(:)  -- coord_struct, optional, allocatable: Enlarged orbit array.
!-

subroutine insert_element (lat, insert_ele, ix_ele, ix_branch, orbit)

use bmad_interface, except_dummy => insert_element

implicit none

type (lat_struct), target :: lat
type (ele_struct)  insert_ele, insert_ele_copy
type (ele_struct), pointer :: inserted_ele, ele0, ele2
type (branch_struct), pointer :: branch, branch2
type (control_struct), pointer :: con
type (lat_ele_loc_struct), pointer :: loc
type (coord_struct), optional, allocatable :: orbit(:)

integer ix_ele, ix, ix_br, i, j
integer, optional :: ix_branch

character(16), parameter :: r_name = 'insert_element'

logical err_flag

! transfer_ele is fast since it reuses storage.

ix_br = integer_option(0, ix_branch)
branch => lat%branch(ix_br)
insert_ele_copy = insert_ele   ! In case insert_ele is an element in the lattice
call upcase_string(insert_ele_copy%name)

branch%n_ele_max = branch%n_ele_max + 1
if (branch%n_ele_max > ubound(branch%ele, 1)) call allocate_lat_ele_array(lat, ix_branch = ix_br)

do ix = branch%n_ele_max-1, ix_ele, -1
  call transfer_ele (branch%ele(ix), branch%ele(ix+1))
  branch%ele(ix+1)%ix_ele = ix+1
enddo

! Enlarge orbit array

if (present(orbit)) then
  call reallocate_coord(orbit, branch%n_ele_max)
  orbit(ix_ele:branch%n_ele_max) = orbit(ix_ele-1:branch%n_ele_max-1)
endif

! branch%ele(ix_ele) pointers need to be nullified since they now point to
! the same memory as branch%ele(ix_ele+1)

inserted_ele => branch%ele(ix_ele)
call deallocate_ele_pointers (inserted_ele, nullify_only = .true.)
inserted_ele = insert_ele_copy
inserted_ele%ix_ele    = ix_ele
inserted_ele%ix_branch = ix_br
inserted_ele%branch => branch

! Correct the control info

do ix = 1, lat%n_control_max
  con => lat%control(ix)
  if (con%slave%ix_ele >= ix_ele .and. con%slave%ix_branch == ix_br)  con%slave%ix_ele = con%slave%ix_ele + 1
  if (con%lord%ix_ele >= ix_ele .and. ix_br == 0) con%lord%ix_ele = con%lord%ix_ele + 1
enddo

if (ix_ele <= branch%n_ele_track + 1) then
  branch%n_ele_track = branch%n_ele_track + 1
  branch%n_ele_track = branch%n_ele_track
else
  call out_io (s_warn$, r_name, &
                  'YOU ARE INSERTING AN ELEMENT *NOT* INTO THE TRACKING PART OF THE LATTICE!', &
                  'ELEMENT: ' // insert_ele_copy%name)
endif

if (branch%ix_to_ele >= ix_ele) branch%ix_to_ele = branch%ix_to_ele + 1

do ix = 0, ubound(lat%branch, 1)
  branch2 => lat%branch(ix)
  if (branch2%ix_from_ele >= ix_ele .and. branch2%ix_from_branch == ix_br) &
                                                  branch2%ix_from_ele = branch2%ix_from_ele + 1
enddo

! Update fork info

do ix = 0, ubound(lat%branch, 1)
  branch2 => lat%branch(ix)
  do i = 1, branch2%n_ele_track
    ele0 => branch2%ele(i)
    if (ele0%key /= fork$ .and. ele0%key /= photon_fork$) cycle
    ele2 => pointer_to_next_ele(ele0, 1, follow_fork = .true.)
    if (ele2%ix_branch /= ix_br) cycle
    if (ele2%ix_ele >= ix_ele) ele0%value(ix_to_element$) = ele2%ix_ele + 1
  enddo
enddo

call set_flags_for_changed_attribute (inserted_ele)
call nametable_add(lat%nametable, insert_ele%name, ele_nametable_index(inserted_ele))

end subroutine
