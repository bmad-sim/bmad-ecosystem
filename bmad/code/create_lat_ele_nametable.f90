!+
! Subroutine create_lat_ele_nametable (lat, nametable)
!
! Routine to create a sorted nametable of element names for a lattice.
!
! Input:
!   lat         -- lat_struct: Lattice.
!
! Ouput:
!   nametable   -- nametable_struct: Nametable of the elment names
!-

subroutine create_lat_ele_nametable (lat, nametable)

use bmad_routine_interface, dummy => create_lat_ele_nametable

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (nametable_struct), target :: nametable
type (nametable_struct), pointer :: nt

integer ie, ib, n, n_tot

! Allocate arrays

nt => nametable
n_tot = 0
do ib = 0, ubound(lat%branch, 1)
  n_tot = n_tot + lat%branch(ib)%n_ele_max + 1
enddo

call re_allocate2 (nt%name, 0, n_tot-1, .false.)
call re_allocate2 (nt%index, 0, n_tot-1, .false.)
nt%n_min = 0
nt%n_max = n_tot - 1

! And sort

n_tot = 0
do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  do ie = 0, branch%n_ele_max
    nt%name(n_tot) = branch%ele(ie)%name
    n_tot = n_tot + 1
  enddo
enddo

call indexer (nt%name(0:nt%n_max), nt%index(0:nt%n_max))
nt%index(0:nt%n_max) = nt%index(0:nt%n_max) - 1

end subroutine
