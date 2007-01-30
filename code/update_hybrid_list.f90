!+
! Subroutine update_hybrid_list (lat, n_in, use_ele)
!
! Subroutine to add elements to the use_ele list needed by the routine
! make_hybrid_lat.
!
! Use_ele is a list of elements that should not be hyberdized. This list
! is used by the subroutine make_hybrid_lat. If an element is to be
! used (not hyberdized) then the associated lord and slave elements need
! to be also added to the use_ele list. This routine does that bookkeeping.
! for a single element use_ele(n_in).
!
! Modules needed:
!   use bmad
!
! Input:
!   lat -- lat_struct: Input lat structure.
!   n_in -- Integer: use_ele(n_in) is the element whose associated lord and
!             slave elements are to be added to use_ele.
!
! Output:
!   USE_ELE(:) -- Logical: list of lat elements to be not hyberdized.
!                   This is used with make_hybrid_lat.
!
! Note: If use_ele(n_in) = .false. then no updating is done
!-

#include "CESR_platform.inc"

recursive subroutine update_hybrid_list (lat, n_in, use_ele)

  use bmad_struct

  implicit none

  type (lat_struct)  lat

  logical use_ele(:)

  integer ix, n_in, i, j

! see if any work needs to be done

  if (.not. use_ele(n_in)) return

! now go through and put controlled elements in the list and make sure
! all appropriate controllers are on the list

  do i = lat%ele(n_in)%ix1_slave, lat%ele(n_in)%ix2_slave
    ix = lat%control(i)%ix_slave
    if (.not. use_ele(ix)) then
      use_ele(ix) = .true.
      call update_hybrid_list (lat, ix, use_ele)
    endif
  enddo

  do i = lat%ele(n_in)%ic1_lord, lat%ele(n_in)%ic2_lord
    j= lat%ic(i)
    ix = lat%control(j)%ix_lord
    if (.not. use_ele(ix)) then
      use_ele(ix) = .true.
      call update_hybrid_list (lat, ix, use_ele)
    endif
  enddo

end subroutine
