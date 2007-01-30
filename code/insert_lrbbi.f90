!+
! Subroutine insert_LRBBI (lat, lat_oppos, cross_positions, ix_LRBBI)
!            
! Uses a lat and a list of parasitic crossing points to create and insert
!   beambeam elements at each crossing point. Returns the new lat, and a list
!   of indices of the new elements.
!
! Modules needed:
!   use bmad
!
! Input:
!   Lat               -- lat_struct: lat to be updated with lrbbi.
!   Lat_oppos         -- lat_struct: lat with positions of opposite bunches
!   cross_positions(:) -- Real(rp): array of parasitic crossing positions
!
! Output:
!   Lat               -- Updated lat within beambeam elements inserted
!   ix_LRBBI(1:*)      -- Integer: array of indices of inserted elements
!-

#include "CESR_platform.inc"

subroutine insert_LRBBI (lat, oppos_lat, cross_positions, ix_LRBBI)

  use bmad_struct
  use bmad_interface, except => insert_LRBBI

  implicit none

  type (lat_struct)  lat
  type (lat_struct) oppos_lat
  type (ele_struct), pointer ::  insert_ele(:)

  real(rp), dimension(:), intent(inout) :: cross_positions
  integer, dimension(:), intent(inout) :: ix_LRBBI

  real(rp) :: s_split
  character(40) :: call_it
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

! Add a beambeam element to the lat at each crossing point.

  do i=1, size(cross_positions)

    s_split = cross_positions(i) * lat%param%total_length 
            !puts cross_positions in meters

    if (s_split == 0) cycle

    ix_split=0

    call split_lat(oppos_lat, s_split, ix_split_oppos, split_done)
    call lat_make_mat6(oppos_lat, -1)

    call split_lat(lat, s_split, ix_split, split_done)
    call lat_make_mat6(lat, -1)

    ix_ele = ix_split + 1

    call init_LRBBI(lat, oppos_lat, insert_ele(i), ix_split, ix_split_oppos)

    insert_ele(i)%s = s_split

    ix_LRBBI(i) = ix_ele

    call insert_element(lat, insert_ele(i), ix_ele)

    call lat_make_mat6(lat, ix_ele)

  enddo

  deallocate(insert_ele, stat=ierr)
  if (ierr .ne.0) then
    print*, "INSERT_ELE: DEALLOCATION REQUEST DENIED."
    call err_exit
  endif

end subroutine
