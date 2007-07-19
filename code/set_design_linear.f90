!+
! Subroutine set_design_linear (lat)
!
! Subroutine to set only those elements on that constitute the "design"
! lattice. That is, only quadrupoles, bends and wigglers will be set on.
!
! Modules needed:
!   use bmad
!
! Output:
!   lat -- lat_struct: Lat structure with quads, bends and wigglers only 
!           set on.
!-

#include "CESR_platform.inc"


subroutine set_design_linear (lat)

  use bmad_struct
  use bmad_interface, except_dummy => set_design_linear

  implicit none

  type (lat_struct) lat
  type (coord_struct), allocatable :: orb(:)

  integer i, key

  allocate (orb(0:ubound(lat%ele, 1)))

!

  do i = 1, lat%n_ele_track
    key = lat%ele(i)%key
    if (key == quadrupole$ .or. key == sbend$ .or. key == wiggler$) then
      lat%ele(i)%is_on = .true.
    else
      lat%ele(i)%is_on = .false.
    endif
  enddo

  call lat_make_mat6 (lat, -1, orb)

  deallocate(orb)

end subroutine
