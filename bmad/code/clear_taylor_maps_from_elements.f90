!+
! Subroutine clear_taylor_maps_from_elements (lat)
!
! Routine to clear the taylor maps from all the elements in a lattice.
! Exception: Taylor elements are not cleared since their maps were set in the lattice file.
!
! Input:
!   lat     -- lat_struct: Lattice
!
! Output:
!   lat     -- lat_struct: Lattice with all maps cleared
!-

subroutine clear_taylor_maps_from_elements (lat)

use taylor_mod, dummy => clear_taylor_maps_from_elements

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
integer ib, ie

!

do ib = 0, ubound(lat%branch,1)
  do ie = 1, lat%branch(ib)%n_ele_max
    ele => lat%branch(ib)%ele(ie)
    if (ele%key == taylor$) cycle

    if (associated(ele%taylor(1)%term)) then
      call kill_taylor(ele%taylor)
      call kill_taylor(ele%spin_taylor)
      ele%spin_q(0,0) = real_garbage$
    endif
  enddo
enddo

end subroutine
