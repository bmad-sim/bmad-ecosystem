!+
! Subroutine create_lat_ele_sorted_nametable (lat, nametable)
!
! Routine to create a sorted nametable of element names for a lattice.
!
! Input:
!   lat         -- lat_struct: Lattice.
!
! Ouput:
!   nametable   -- lat_nametable_struct: Nametable of the elment names
!-

subroutine create_lat_ele_sorted_nametable (lat, nametable)

use bmad_routine_interface, dummy => create_lat_ele_sorted_nametable

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (lat_nametable_struct) nametable

integer ie, ib, n, n_tot
character(40), allocatable :: ele_name(:)

!

if (allocated(nametable%branch)) deallocate(nametable%branch)
allocate (nametable%branch(0:ubound(lat%branch, 1)))

n_tot = 0
do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  n = branch%n_ele_max
  allocate (nametable%branch(ib)%indexx(n+1))
  call indexx (branch%ele(0:n)%name, nametable%branch(ib)%indexx)
  nametable%branch(ib)%indexx = nametable%branch(ib)%indexx - 1
  n_tot = n_tot + n + 1
enddo

allocate (ele_name(n_tot))
allocate (nametable%all(n_tot), nametable%all_indexx(n_tot))

n_tot = 0
do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  n = branch%n_ele_max
  do ie = 0, n
    ele_name(n_tot+ie+1) = branch%ele(ie)%name
    nametable%all(n_tot+ie+1)%ele => branch%ele(ie)
  enddo
  n_tot = n_tot + n + 1
enddo

call indexx (ele_name, nametable%all_indexx)

end subroutine
