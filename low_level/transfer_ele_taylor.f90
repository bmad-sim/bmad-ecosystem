!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine transfer_ele_taylor (ele_in, ele_out, taylor_order)
!
! Subroutine to transfer a Taylor map from one element to another.
!
! Input:
!   ele_in       -- Ele_struct: Element with the Taylor map.
!   taylor_order -- Integer, optional: Order to truncate the Taylor map at.
!
! Output:
!   ele_out      -- Ele_struct: Element receiving the Taylor map truncated to
!                     order taylor_order.
!-

subroutine transfer_ele_taylor (ele_in, ele_out, taylor_order)

use bmad_struct

implicit none

type (ele_struct) ele_in, ele_out
integer, optional :: taylor_order
integer it, ix, k 

!

do it = 1, 6

  if (present(taylor_order)) then
    ix = 0
    do k = 1, size(ele_in%taylor(it)%term)
      if (sum(ele_in%taylor(it)%term(k)%expn(:)) > taylor_order) cycle
      ix = ix + 1
    enddo
  else
    ix = size(ele_in%taylor(it)%term)
  endif

  if (.not. associated(ele_out%taylor(it)%term)) allocate (ele_out%taylor(it)%term(ix))
  if (size(ele_out%taylor(it)%term) /= ix) allocate (ele_out%taylor(it)%term(ix))

  ix = 0
  do k = 1, size(ele_in%taylor(it)%term)
    if (present(taylor_order)) then
      if (sum(ele_in%taylor(it)%term(k)%expn(:)) > taylor_order) cycle
    endif
    ix = ix + 1
    ele_out%taylor(it)%term(ix) = ele_in%taylor(it)%term(k)
  enddo

enddo

ele_out%taylor(:)%ref = ele_in%taylor(:)%ref

end subroutine transfer_ele_taylor

