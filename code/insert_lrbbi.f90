!+
! Subroutine insert_LRBBI (ring, ring_oppos, cross_positions, ix_LRBBI)
!            
! Uses a ring and a list of parasitic crossing points to create and insert
!   beambeam elements at each crossing point. Returns the new ring, and a list
!   of indices of the new elements.
!
! Modules needed:
!   use bmad
!
! Input:
!   Ring               -- Ring_struct: ring to be updated with lrbbi.
!   Ring_oppos         -- Ring_struct: ring with positions of opposite bunches
!   cross_positions(:) -- Real(rp): array of parasitic crossing positions
!
! Output:
!   Ring               -- Updated ring within beambeam elements inserted
!   ix_LRBBI(1:*)      -- Integer: array of indices of inserted elements
!-

#include "CESR_platform.inc"

subroutine insert_LRBBI (ring, oppos_ring, cross_positions, ix_LRBBI)

  use bmad_struct
  use bmad_interface, except => insert_LRBBI

  implicit none

  type (ring_struct)  ring
  type (ring_struct) oppos_ring
  type (ele_struct), dimension(:), allocatable ::  insert_ele

  real(rp), dimension(:), intent(inout) :: cross_positions
  integer, dimension(:), intent(inout) :: ix_LRBBI

  real(rp) :: s_split
  character*16 :: call_it
  integer :: ix_ele, ix_split, ix_split_oppos, ierr, i
  logical :: split_done

  integer :: loc_smallest, end
  real(rp) :: smallest
  integer, dimension(1) :: minloc_array

! Order the crossing points from earliest to latest.

  end = size(cross_positions)

  do i = 1, end
    smallest = minval(cross_positions(i:end))
    minloc_array = minloc(cross_positions(i:end))
    loc_smallest = i - 1 + minloc_array(1)
    if (loc_smallest .ne. i) then
      cross_positions(loc_smallest) = cross_positions(i)
      cross_positions(i) = smallest
    endif
  enddo

  allocate(insert_ele(1:end), stat=ierr)
  if (ierr .ne. 0) then
    print*, "INSERT_ELE: ALLOCATION REQUEST DENIED."
    call err_exit
  endif

  call_it = 'lrbbi'

! Add a beambeam element to the ring at each crossing point.

  do i=1, size(cross_positions)

    s_split = cross_positions(i) * ring%param%total_length 
            !puts cross_positions in meters

    if (s_split == 0) cycle

    ix_split=0

    call split_ring(oppos_ring, s_split, ix_split_oppos, split_done)
    call ring_make_mat6(oppos_ring, -1)

    call split_ring(ring, s_split, ix_split, split_done)
    call ring_make_mat6(ring, -1)

    ix_ele = ix_split + 1

    call init_LRBBI(ring, oppos_ring, insert_ele(i), ix_split, ix_split_oppos)

    insert_ele(i)%s = s_split

    ix_LRBBI(i) = ix_ele

    call insert_element(ring, insert_ele(i), ix_ele)

    call ring_make_mat6(ring, ix_ele)

  enddo

  deallocate(insert_ele, stat=ierr)
  if (ierr .ne.0) then
    print*, "INSERT_ELE: DEALLOCATION REQUEST DENIED."
    call err_exit
  endif

end subroutine
