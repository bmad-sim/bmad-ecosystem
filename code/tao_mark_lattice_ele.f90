!+
! Subroutine tao_mark_lattice_ele (lat)
!
! Routine to mark every ele%ix_pointer in lat to denote it's place in the lattice.
! This is used so that if an element's position is changed (say via superposition) then Tao can 
! identify where the element has moved to.
!
! Input:
!   lat   -- lat_struct: Input lattice
!
! Output:
!   lat   -- lat_struct: Lattice with elements marked.
!-

subroutine tao_mark_lattice_ele (lat)

use tao_struct

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele

integer ib, ie

!

do ib = 0, ubound(lat%branch, 1)
  do ie = 0, ubound(lat%branch(ib)%ele, 1)
    ele => lat%branch(ib)%ele(ie)
    ele%ix_pointer = 1000000*ib + ie
  enddo
enddo

end subroutine
