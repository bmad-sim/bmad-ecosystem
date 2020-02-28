!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine get_slave_list (lord, slave_list, n_slave)
!
! Subroutine to get the list of slaves for a lord element.
!
! This is a list of ultimate slaves. That is, slaves in the tracking part 
! of the lattice. Thus if the element lord controls an
! element lord1 which controlls an element lord2, then lord2 will
! show up in the slave_list but lord1 will not.
!
! If the "lord" element does not have any slaves, 
! then the slave_list will just be the lord element.
!
! This routine will increase the size of slave_list if needed but will
! not decrease it.
!
! Input:
!   lord  -- Ele_struct: The lord element.
!
! Output:
!   slaves(:) -- Ele_pointer_struct, allocatable :: Array of slaves.
!   n_slave   -- Integer: Number of slaves.
!-

subroutine get_slave_list (lord, slaves, n_slave)

use equal_mod, dummy => get_slave_list

implicit none

type (ele_struct), target :: lord
type (ele_pointer_struct), allocatable :: slaves(:)

integer n_slave

!

if (lord%n_slave == 0) then
  n_slave = 1
  call re_allocate_eles (slaves, 1)
  slaves(1)%ele => lord
  return
endif

!

n_slave = 0
if (.not. allocated(slaves)) call re_allocate_eles (slaves, lord%n_slave)

call get_slaves (lord)

!--------------------------------------------------------------------------
contains

recursive subroutine get_slaves (lord)

type (ele_struct) lord
type (ele_struct), pointer :: slave
integer i, ix

!

do i = 1, lord%n_slave
  slave => pointer_to_slave(lord, i)
  if (slave%n_slave > 0) then
    call get_slaves (slave)
  else
    n_slave = n_slave + 1
    if (n_slave > size(slaves)) call re_allocate_eles(slaves, n_slave + 4, .true.)
    slaves(n_slave)%ele => slave
  endif
enddo

end subroutine get_slaves

end subroutine get_slave_list
