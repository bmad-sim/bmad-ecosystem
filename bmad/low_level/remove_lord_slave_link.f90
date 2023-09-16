!+
! Subroutine remove_lord_slave_link (lord, slave)
!
! Routine to remove the pointers between a lord and a slave.
! Note: This routine will not modify lord%lord_status and slave%slave_status.
!
! Input:
!   lord  -- ele_struct: Lord element
!   slave -- ele_struct: Slave element
!
! Output:
!   lord  -- ele_struct: Lord element with link info removed
!   slave -- ele_struct: Slave element with link info removed
!   lat   -- Lattice_struct: Lattice obtaind from lord%lat and slave%lat pointers.
!     %control(:)  -- Array modified to remove lord/slave link.
!     %ic(:)       -- Array modified to remove lord/slave link.
!-

subroutine remove_lord_slave_link (lord, slave)

use bmad_struct
implicit none

type (ele_struct), target :: lord, slave
type (ele_struct), pointer :: ele
type (lat_struct), pointer :: lat
type (branch_struct), pointer :: branch

integer ic_out, icon_out, ib, i, n, ic2_lord

!

if (.not. associated(lord%branch, slave%branch)) call err_exit  ! Should not be

lat => lord%branch%lat

! Find lat%control(:) and lat%ic(:) elements associated with this link

ic2_lord = slave%ic1_lord + slave%n_lord - 1
do ic_out = slave%ic1_lord, ic2_lord
  icon_out = lat%ic(ic_out)
  if (lat%control(icon_out)%lord%ix_ele == lord%ix_ele) exit
enddo

if (icon_out == ic2_lord+1) call err_exit ! Should not be

! Compress lat%control and lat%ic arrays.

n = lat%n_control_max
lat%control(icon_out:n-1) = lat%control(icon_out+1:n)
lat%n_control_max = n - 1

n = lat%n_ic_max
lat%ic(ic_out:n-1) = lat%ic(ic_out+1:n)
lat%n_ic_max = n - 1

do i = 1, n - 1
  if (lat%ic(i) > icon_out) lat%ic(i) = lat%ic(i) - 1
enddo

! Correct info in all elements

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)

  do i = 1, branch%n_ele_max
    ele => branch%ele(i)

    if (ele%ix1_slave >  icon_out) ele%ix1_slave = ele%ix1_slave - 1
    if (ele%ic1_lord >  ic_out) ele%ic1_lord = ele%ic1_lord - 1
  enddo
enddo

slave%n_lord = slave%n_lord - 1
if (slave%n_lord == 0) then
  slave%ic1_lord = 0
  slave%slave_status = free$
endif

lord%n_slave = lord%n_slave - 1
if (lord%n_slave == 0) then
  lord%ix1_slave = 0
  lord%lord_status = not_a_lord$
endif

end subroutine remove_lord_slave_link 
