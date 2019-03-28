!+
! Subroutine slice_lattice (lat, ele_list, error)
!
! Routine to discard from the lattice all elements not in ele_list.
! Note controllers that control elements that remain will not be cut.
!
! Input:
!   lat       -- lat_struct: Lattice
!   ele_list  -- character(*): List of elements to retain
!
! Output:
!   lat       -- lat_struct: Lattice with unwanted elements sliced out.
!   error     -- logical: Set True if there is an error Set False if not.
!-

subroutine slice_lattice (lat, ele_list, error)

use bmad, dummy => slice_lattice

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_pointer_struct), allocatable :: eles(:)

integer ie, ib, n_loc
logical error, err

character(*) ele_list
character(*), parameter :: r_name = 'slice_lattice'

!

error = .true.

call lat_ele_locator (ele_list, lat, eles, n_loc, err)
if (err) return
if (n_loc == 0) then
  call out_io (s_error$, r_name, 'NO LATTICE ELEMENTS FOUND: ' // ele_list, 'LATTICE NOT SLICED.')
  return
endif

! Use ele%ixx to tag elements to be deleted which is everything not in eles list.

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  branch%ele%ixx = -1
enddo

do ie = 1, n_loc
  eles(ie)%ele%ixx = 0  ! Do not delete
enddo

! Now go through and save all controllers that can be saved

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  do ie = 1, branch%n_ele_max
    if (branch%ele(ie)%ixx == -1) cycle
    call add_back_controllers(branch%ele(ie))
  enddo
enddo


! And remove

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  do ie = 1, branch%n_ele_max
    if (branch%ele(ie)%ixx == -1) branch%ele(ie)%key = -1
  enddo
enddo

call remove_eles_from_lat (lat)
call lattice_bookkeeper (lat)

error = .false.

!-------------------------------------------------------------------------------------
contains

recursive subroutine add_back_controllers (ele)

type (ele_struct) ele
type (ele_struct), pointer :: ele2
integer i

do i = 1, ele%n_lord
  ele2 => pointer_to_lord (ele, i)
  ele2%ixx = 0
  call add_back_controllers(ele2)
enddo

end subroutine add_back_controllers

end subroutine slice_lattice
