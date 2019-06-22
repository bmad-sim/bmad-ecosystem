!+
! Subroutine tao_control_tree_list (ele, tree)
!
! Routine to create a list of slaves and lords of a given lattice element.
! Also included is all lords of any element that is in the list.
! The list is semi-ordered with slave elements appearing before lord elements.
!
! Input:
!   ele     -- ele_struct: Lattice element to start at.
!
! Ouput:
!   tree(:) -- ele_pointer_struct, allocatable: Array of elements.
!-

subroutine tao_control_tree_list (ele, tree)

use bmad

implicit none

type (ele_struct) ele
type (ele_pointer_struct), allocatable :: tree(:)

integer n_found, i

! First walk downwards and get a list of ultimate slaves

n_found = 0
call walk_downwards (ele, tree, n_found)

! And now walk upwards

do i = 1, n_found
  call walk_upwards (tree(i)%ele, tree, n_found)
enddo

call re_allocate_eles(tree, n_found, .true., .true.)

!------------------------------
contains

recursive subroutine walk_downwards (ele, tree, n_found)

type (ele_struct), target :: ele
type (ele_struct), pointer :: slave
type (ele_pointer_struct), allocatable :: tree(:)

integer n_found, i

!

if (ele%n_slave + ele%n_slave_field == 0) then  ! Ultimate slave
  call re_allocate_eles(tree, n_found+10, .true.)
  n_found = n_found + 1
  tree(n_found)%ele => ele

else  ! loop over slaves
  do i = 1, ele%n_slave + ele%n_slave_field
    slave => pointer_to_slave(ele, i)
    call walk_downwards(ele, tree, n_found)
  enddo
endif

end subroutine walk_downwards

!------------------------------
! contains

recursive subroutine walk_upwards (ele, tree, n_found)

type (ele_struct), target :: ele
type (ele_struct), pointer :: lord
type (ele_pointer_struct), allocatable :: tree(:)

integer n_found, i

! Put element in list if not already there.

do i = 1, n_found
  if (associated(tree(i)%ele, ele)) exit
  if (i == n_found) then
      call re_allocate_eles(tree, n_found+10, .true.)
      n_found = n_found + 1
      tree(n_found)%ele => ele
  endif
enddo

! Now go through lords

do i = 1, ele%n_lord
  lord => pointer_to_lord(ele, i)
  call walk_upwards (lord, tree, n_found)
enddo

end subroutine walk_upwards

end subroutine tao_control_tree_list
