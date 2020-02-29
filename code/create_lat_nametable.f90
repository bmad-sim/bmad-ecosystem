!+
! Subroutine create_lat_nametable (lat, nametable)
!
! Routine to create a sorted nametable of element names for a lattice.
!
! Input:
!   lat         -- lat_struct: Lattice.
!
! Ouput:
!   nametable   -- lat_nametable_struct: Nametable of the elment names
!-

subroutine create_lat_nametable (lat, nametable)

use bmad_routine_interface, dummy => create_lat_nametable

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (lat_nametable_struct) nametable

integer ie, ib, n, n_tot

! Allocate arrays

n_tot = 0
do ib = 0, ubound(lat%branch, 1)
  n_tot = n_tot + lat%branch(ib)%n_ele_max + 1
enddo

allocate (nametable%name(n_tot), nametable%indexx(n_tot))

! And sort

n_tot = 0
do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  n = branch%n_ele_max
  do ie = 0, n
    nametable%name(n_tot+ie+1) = branch%ele(ie)%name
  enddo
  n_tot = n_tot + n + 1
enddo

call indexx (nametable%name, nametable%indexx)

end subroutine
