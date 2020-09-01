!+
! Subroutine ptc_bookkeeper (lat)
!
! Routine to do bookkeeping on the associated ptc layout if it exists.
!
! Input:
!   lat -- lat_struct: Bmad lattice.
!
! Output:
!   lat -- lat_struct:
!-

subroutine ptc_bookkeeper (lat)

use ptc_layout_mod, dummy => ptc_bookkeeper

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele

integer i, j

!

do i = lbound(lat%branch, 1), ubound(lat%branch, 1)
  branch => lat%branch(i)
  if (.not. associated(branch%ptc%m_t_layout)) cycle
  do j = 0, branch%n_ele_track
    ele => branch%ele(j)
    if (.not. associated (ele%ptc_fibre)) cycle
    if (ele%bookkeeping_state%ptc /= stale$) cycle
    call update_fibre_from_ele (ele)
  enddo
enddo

end subroutine ptc_bookkeeper

