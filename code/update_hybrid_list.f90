!+
! Subroutine update_hybrid_list (lat, n_in, keep_ele, keep_overlays_and_groups)
!
! Subroutine to add elements to the keep_ele list needed by the routine
! make_hybrid_lat.
!
! Keep_ele is a list of elements that should appear in the hybridized lattice. 
! This list is used by the subroutine make_hybrid_lat. If an element is to be
! kept (not hybridized) then the associated lord and slave elements need
! to be also added to the keep_ele list. This routine does that bookkeeping
! for a single element keep_ele(n_in).
!
! Modules needed:
!   use bmad
!
! Input:
!   lat  -- lat_struct: Input lat structure.
!   n_in -- Integer: keep_ele(n_in) is the element whose associated lord and
!             slave elements are to be added to keep_ele.
!   keep_overlays_and_groups
!        -- Logical, optional: If present and False then overlay and group elements
!             are kept. Default is True.
!
! Output:
!   keep_ele(:) -- Logical: list of lat elements to be not hybridized.
!                   This is used with make_hybrid_lat.
!-

recursive subroutine update_hybrid_list (lat, n_in, keep_ele, keep_overlays_and_groups)

use bmad_struct

implicit none

type (lat_struct)  lat

logical keep_ele(:)
logical, optional :: keep_overlays_and_groups

integer ix, n_in, i, j

! 

keep_ele(n_in) = .true.

! now go through and put controlled elements in the list and make sure
! all appropriate controllers are on the list

do i = lat%ele(n_in)%ix1_slave, lat%ele(n_in)%ix2_slave
  ix = lat%control(i)%ix_slave
  if (keep_ele(ix)) cycle
  call update_hybrid_list (lat, ix, keep_ele)
enddo

do i = lat%ele(n_in)%ic1_lord, lat%ele(n_in)%ic2_lord
  j = lat%ic(i)
  ix = lat%control(j)%ix_lord

  if (keep_ele(ix)) cycle  ! If already set then don't need to do it again
  if (.not. logic_option(.true., keep_overlays_and_groups) .and. &
                    (lat%ele(ix)%lord_status == overlay_lord$ .or. &
                    lat%ele(ix)%lord_status == group_lord$)) cycle

  call update_hybrid_list (lat, ix, keep_ele, keep_overlays_and_groups)
enddo

end subroutine
