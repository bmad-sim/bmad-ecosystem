!+
! Subroutine get_element_lord_list (lat, ix_ele, lord_list, n_lord)
!
! Subroutine to get the list of lords for an element.
!
! Note: This routine will increase the size of lord_list if needed 
! but will not decrease it.
!
! Modules needed:
!   use bmad
!
! Input:
!   lat     -- lat_struct: Lattice
!   ix_ele  -- Integer: Index element in the lattice.
!
! Output:
!   lord_list(:) -- Integer, allocatable: Array of lord elements.
!   n_lord       -- Integer: Number of lords.
!-

subroutine get_element_lord_list (lat, ix_ele, lord_list, n_lord)

  use bmad_struct
  use bmad_interface

  implicit none

  type (lat_struct) lat

  integer ix_ele, n_lord
  integer i, i2, ix
  integer, allocatable :: lord_list(:)

!

  n_lord = lat%ele(ix_ele)%n_lord
  call re_allocate (lord_list, n_lord)

  do i = 1, n_lord
    i2 = lat%ic(i + lat%ele(ix_ele)%ic1_lord - 1)
    ix = lat%control(i2)%ix_lord
    lord_list(i) = ix
  enddo

end subroutine
