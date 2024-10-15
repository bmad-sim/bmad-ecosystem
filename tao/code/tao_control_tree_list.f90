!+
! Subroutine tao_control_tree_list (ele, tree)
!
! Routine to create a list of slaves and lords of a given lattice element.
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
type (ele_pointer_struct), allocatable :: tree(:), tree2(:)

integer n_found, n_found2, i

! Walk downwards and upwards to get a list of slaves and lords

n_found = 0
call walk_downwards (ele, tree, n_found)

n_found2 = 0
call walk_upwards (ele, tree2, n_found2)

! Combine lists. ele is common to both lists so do not double count.

call re_allocate_eles(tree, n_found+n_found2-1, .true., .true.)
do i = 2, n_found2
  tree(n_found+n_found2+1-i) = tree2(i)
enddo

!------------------------------
contains

recursive subroutine walk_downwards (ele, tree, n_found)

type (ele_struct), target :: ele
type (ele_struct), pointer :: slave
type (ele_pointer_struct), allocatable :: tree(:)

integer n_found, i

! Put element in list if not already there.

do i = 1, n_found
  if (associated(tree(i)%ele, ele)) return
enddo

if (n_found + 1 > size(tree)) call re_allocate_eles(tree, n_found+10, .true.)
n_found = n_found + 1
tree(n_found)%ele => ele

! Loop over slaves

do i = 1, ele%n_slave + ele%n_slave_field
  slave => pointer_to_slave(ele, i)
  call walk_downwards(slave, tree, n_found)
enddo

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
  if (associated(tree(i)%ele, ele)) return
enddo

if (n_found + 1 > size(tree)) call re_allocate_eles(tree, n_found+10, .true.)
n_found = n_found + 1
tree(n_found)%ele => ele

! Now go through lords

do i = 1, ele%n_lord
  lord => pointer_to_lord(ele, i)
  call walk_upwards (lord, tree, n_found)
enddo

end subroutine walk_upwards

end subroutine tao_control_tree_list
