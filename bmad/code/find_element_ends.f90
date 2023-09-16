!+
! Subroutine find_element_ends (ele, ele1, ele2, ix_multipass)
!
! Subroutine to find the end points of an element in the tracking part of the 
! lattice.
!
! Note: This routine does not make sense in a number of cases. For example,
! for multipass_lord elements.
!
! Input:
!   ele          -- Ele_struct: Element to find the ends for.
!   ix_multipass -- integer, optional: Which multipass pass to follow. Default is 1.
!                     This is ignored if there is no multipass elements.
!
! Output:
!   ele1        -- Ele_struct, pointer:  Pointer to element just before ele. 
!   ele2        -- Ele_struct, pointer:  Pointer to ele itself or the last sub-element within ele.
!
! Note: ele1 and ele2 will be nullified if ele is in the lord 
!       part of the lattice and does not have any slaves.
!
! Note: For an element in the tracking part of the lattice:
!       ele1%ix_ele = ele%ix_ele - 1
!       ele2        => ele
!  Exception: For Beginning element (index 0), ele1 => ele
!-

subroutine find_element_ends (ele, ele1, ele2, ix_multipass)

use pointer_to_ele_mod, except_dummy => find_element_ends

implicit none
                                                       
type (ele_struct), target :: ele
type (ele_struct), pointer :: ele1, ele2

integer, optional :: ix_multipass
integer ix_start, ix_end, ix_start_branch, ix_end_branch, ix_mult

! 

nullify (ele1, ele2)
ix_mult = integer_option(1, ix_multipass)

if (ele%n_slave == 0) then

  if (ele%ix_ele > ele%branch%n_ele_track) return

  if (ele%ix_ele == 0) then
    ele1 => ele
    ele2 => ele

  else
    ele1 => pointer_to_next_ele (ele, -1)
    ele2 => ele
  endif

elseif (ele%lord_status == super_lord$ .or. ele%key == girder$ .or. &
        ele%lord_status == multipass_lord$) then
  ele1 => ele
  do
    if (ele1%n_slave == 0) exit
    if (ele1%lord_status == multipass_lord$) then
      if (ix_mult < 1 .or. ix_mult > ele1%n_slave) return
      ele1 => pointer_to_slave(ele1, ix_mult)
    else
      ele1 => pointer_to_slave(ele1, 1)
    endif
  enddo
  ele1 => pointer_to_next_ele (ele1, -1)

  ele2 => ele
  do
    if (ele2%n_slave == 0) exit
    if (ele2%lord_status == multipass_lord$) then
      if (ix_mult < 1 .or. ix_mult > ele2%n_slave) return
      ele2 => pointer_to_slave(ele2, ix_mult)
    else
      ele2 => pointer_to_slave(ele2, ele2%n_slave)
    endif
  enddo

! For overlays and groups: The idea is to look at all the slave elements in the tracking 
! part of the lattice and find the minimum and maximum element indexes.
! An element with a greater %ix_branch is always considered to be greater independent of %ix_ele.
! The small complication is that overlays or groups lords can control other overlays or 
! groups, etc.
! So we must "recursively" follow the slave tree.
! ix_slave_array/ix_branch_array holds the list of slaves we need to look at.

else  ! overlay$, group$

  ix_start = 1000000
  ix_start_branch = 100000

  ix_end = 0
  ix_end_branch = -1

  call find_these_ends (ele)

  ele1 => pointer_to_ele (ele%branch%lat, ix_start, ix_start_branch)
  ele2 => pointer_to_ele (ele%branch%lat, ix_end, ix_end_branch)

endif

!-------------------------------------------------------
contains

recursive subroutine find_these_ends (lord)

type (ele_struct), target :: lord

integer i

! If in tracking part of the lattice then see if it is at the edge.

if (lord%ix_ele <= lord%branch%n_ele_track) then

  if (lord%ix_branch < ix_start_branch .or. &
        (lord%ix_branch == ix_start_branch .and. lord%ix_ele - 1 < ix_start)) then
    ix_start = lord%ix_ele - 1
    ix_start_branch = lord%ix_branch
  endif
  if (lord%ix_branch > ix_end_branch .or. &
        (lord%ix_branch == ix_end_branch .and. lord%ix_ele > ix_end)) then
    ix_end = lord%ix_ele 
    ix_end_branch = lord%ix_branch
  endif

  return
endif

! If a lord then follow the slaves.

if (lord%lord_status == multipass_lord$) then
  if (ix_mult < 1 .or. ix_mult > lord%n_slave) return
  call find_these_ends (pointer_to_slave(lord, ix_mult))

else
  do i = 1, lord%n_slave
    call find_these_ends (pointer_to_slave(lord, i))
  enddo
endif

end subroutine find_these_ends

end subroutine find_element_ends
