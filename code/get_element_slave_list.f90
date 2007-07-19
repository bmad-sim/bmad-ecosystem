!+
! Subroutine get_element_slave_list (lat, ix_lord, slave_list, n_slave)
!
! Subroutine to get the list of slaves for a lord element.
!
! This is a list of ultimate slaves. That is, slaves in the tracking part 
! of the lattice. Thus if the element lat%ele(ix_lord) controls an
! element ix_lord2 which controlls an element ix_lord3 then ix_lord3 will
! show up in the slave_list but ix_lord2 will not.
!
! If the ix_lord element is in the tracking part of the lattice 
! then the slave_list will just be that element.
!
! This routine will increase the size of slave_list if needed but will
! not decrease it.
!
! Modules needed:
!   use bmad
!
! Input:
!   lat     -- lat_struct: Lattice
!   ix_lord  -- Integer: Index of lord element in the lattice.
!
! Output:
!   slave_list(:) -- Integer, allocatable: Array of slave elements.
!   n_slave       -- Integer: Number of slaves.
!-

subroutine get_element_slave_list (lat, ix_lord, slave_list, n_slave)

  use bmad_struct
  use bmad_interface, except_dummy => get_element_slave_list

  implicit none

  type (lat_struct) lat

  integer ix_lord, n_slave
  integer, allocatable :: slave_list(:)

!

  n_slave = 0
  if (.not. allocated(slave_list)) call re_allocate (slave_list, 10)

  call get_slaves (lat%ele(ix_lord))

!--------------------------------------------------------------------------
contains

recursive subroutine get_slaves (lord)

  type (ele_struct) lord
  integer i, ix

!

  do i = lord%ix1_slave, lord%ix2_slave
    ix = lat%control(i)%ix_slave
    if (ix > lat%n_ele_track) then
      call get_slaves (lat%ele(ix))
    else
      n_slave = n_slave + 1
      if (n_slave > size(slave_list)) &
                              call re_allocate(slave_list, n_slave + 10)
      slave_list(n_slave) = ix
    endif
  enddo

end subroutine

end subroutine
