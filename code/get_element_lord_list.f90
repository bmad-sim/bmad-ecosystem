!+
! Subroutine get_element_lord_list (lat, ix_slave, lord_list, n_lord)
!
! Subroutine to get the list of lords for a slave element.
!
! Note: This routine will increase the size of lord_list if needed 
! but will not decrease it.
!
! Modules needed:
!   use bmad
!
! Input:
!   lat       -- lat_struct: Lattice
!   ix_slave  -- Integer: Index of slave element in the lattice.
!
! Output:
!   lord_list(:) -- Integer, allocatable: Array of lord elements.
!   n_lord       -- Integer: Number of lords.
!-

subroutine get_element_lord_list (lat, ix_slave, lord_list, n_lord)

  use bmad_struct
  use bmad_interface, except_dummy => get_element_lord_list

  implicit none

  type (lat_struct) lat

  integer ix_slave, n_lord
  integer i, i2, ix
  integer, allocatable :: lord_list(:)

!

  n_lord = lat%ele(ix_slave)%n_lord
  call re_allocate (lord_list, n_lord)

  do i = 1, n_lord
    i2 = lat%ic(i + lat%ele(ix_slave)%ic1_lord - 1)
    ix = lat%control(i2)%ix_lord
    lord_list(i) = ix
  enddo

end subroutine
