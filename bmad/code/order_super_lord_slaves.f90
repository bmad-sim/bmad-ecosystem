!+
! Subroutine order_super_lord_slaves (lat, ix_lord)
!
! Subroutine to properly order the slave element array of a super_lord
!
! Input:
!   lat    -- lat_struct: Lat.
!   ix_lord -- Integer: Index of lord element.
!
! Output
!   lat -- lat_struct: Lat with fixed controls.
!-

subroutine order_super_lord_slaves (lat, ix_lord)

use bmad_interface, except_dummy => order_super_lord_slaves

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: lord, slave
type (control_struct), allocatable :: cs(:)
type (branch_struct), pointer :: branch

integer i, ix_lord, ix1, ix2, ns
integer, allocatable :: ixx(:), iyy(:)

real(rp) ds, s1_lord, l_branch, ds_small
real(rp), allocatable :: s_rel(:)

character(32), parameter :: r_name = 'order_super_lord_slaves'

! Init setup.

lord => lat%ele(ix_lord)

if (lord%lord_status /= super_lord$) then
  call out_io (s_fatal$, r_name, 'ELEMENT NOT A SUPER_LORD')
  if (global_com%exit_on_error) call err_exit
endif

! Make an array of distances between the slave elements and the lord element.
! If the lord wraps around the ends of the lattice, we need to correct the distance.


ns = lord%n_slave
allocate (s_rel(ns), ixx(ns), iyy(ns), cs(ns))
ds_small = bmad_com%significant_length / 2

do i = 1, lord%n_slave
  slave => pointer_to_slave(lord, i)
  branch => lat%branch(slave%ix_branch)
  ds = slave%s - lord%s
  s1_lord = lord%s - lord%value(l$)
  l_branch = branch%param%total_length
  if (s1_lord < branch%ele(0)%s - ds_small .and. s1_lord + l_branch < slave%s) ds = ds - l_branch
  s_rel(i) = ds
enddo

! Sort slaves by distance.

call indexer (s_rel, ixx)

ix1 = lord%ix1_slave; ix2 = lord%ix1_slave+lord%n_slave-1
cs = lat%control(ix1:ix2) 

do i = 1, ns
  lat%control(i+ix1-1) = cs(ixx(i))
  iyy(ixx(i)) = i
enddo

do i = 1, lat%n_ic_max
  if (lat%ic(i) >= ix1 .and. lat%ic(i) <= ix2) then
    lat%ic(i) = iyy(lat%ic(i)+1-ix1) + ix1 - 1
  endif
enddo

deallocate (s_rel, ixx, iyy, cs)

end subroutine
