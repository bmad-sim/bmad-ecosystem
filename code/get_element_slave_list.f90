!+
! Subroutine get_element_slave_list (lat, ix_ele, slave_list, n_slave)
!
! Subroutine to get the list of slaves for an element.
!
! This is a list of ultimate slaves. That is, slaves in the regular part 
! of the lattice. Thus if the element lat%ele_(ix_ele) controls an
! element ix_ele2 which controlls an element ix_ele3 then ix_ele3 will
! show up in the slave_list but ix_ele2 will not.
!
! If the ix_ele element is in the regular part of the lattice 
! then slave_list will just be that element.
!
! This routine will increase the size of slave_list if needed but will
! not decrease it.
!
! Modules needed:
!   use bmad
!
! Input:
!   lat     -- Ring_struct: Lattice
!   ix_ele  -- Integer: Index element in the lattice.
!
! Output:
!   slave_list(:) -- Integer, allocatable: Array of slave elements.
!   n_slave       -- Integer: Number of slaves.
!-

subroutine get_element_slave_list (lat, ix_ele, slave_list, n_slave)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ring_struct) lat

  integer ix_ele, n_slave
  integer, allocatable :: slave_list(:)

!

  n_slave = 0
  if (.not. allocated(slave_list)) call re_allocate (slave_list, 100)

  call get_slaves (lat%ele_(ix_ele))

!--------------------------------------------------------------------------
contains

recursive subroutine get_slaves (lord)

  type (ele_struct) lord
  integer i, ix

!

  do i = lord%ix1_slave, lord%ix2_slave
    ix = lat%control_(i)%ix_slave
    if (ix > lat%n_ele_use) then
      call get_slaves (lat%ele_(ix))
    else
      n_slave = n_slave + 1
      if (n_slave > size(slave_list)) &
                              call re_allocate(slave_list, n_slave + 100)
      slave_list(n_slave) = ix
    endif
  enddo

end subroutine

end subroutine
