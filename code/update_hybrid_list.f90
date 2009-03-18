!+
! Subroutine update_hybrid_list (lat, n_in, keep_ele)
!
! Subroutine to add elements to the keep_ele list needed by the routine
! make_hybrid_lat.
!
! Keep_ele is a list of elements that should appear in the hyberdized lattice. 
! This list is used by the subroutine make_hybrid_lat. If an element is to be
! kept (not hyberdized) then the associated lord and slave elements need
! to be also added to the keep_ele list. This routine does that bookkeeping
! for a single element keep_ele(n_in).
!
! Modules needed:
!   use bmad
!
! Input:
!   lat -- lat_struct: Input lat structure.
!   n_in -- Integer: keep_ele(n_in) is the element whose associated lord and
!             slave elements are to be added to keep_ele.
!
! Output:
!   keep_ele(:) -- Logical: list of lat elements to be not hyberdized.
!                   This is used with make_hybrid_lat.
!
! Note: If keep_ele(n_in) = .false. then no updating is done
!-

recursive subroutine update_hybrid_list (lat, n_in, keep_ele)

use bmad_struct

implicit none

type (lat_struct)  lat

logical keep_ele(:)

integer ix, n_in, i, j

! see if any work needs to be done

if (.not. keep_ele(n_in)) return

! now go through and put controlled elements in the list and make sure
! all appropriate controllers are on the list

do i = lat%ele(n_in)%ix1_slave, lat%ele(n_in)%ix2_slave
  ix = lat%control(i)%ix_slave
  if (.not. keep_ele(ix)) then
    keep_ele(ix) = .true.
    call update_hybrid_list (lat, ix, keep_ele)
  endif
enddo

do i = lat%ele(n_in)%ic1_lord, lat%ele(n_in)%ic2_lord
  j= lat%ic(i)
  ix = lat%control(j)%ix_lord
  if (.not. keep_ele(ix)) then
    keep_ele(ix) = .true.
    call update_hybrid_list (lat, ix, keep_ele)
  endif
enddo

end subroutine
