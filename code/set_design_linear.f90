!+
! Subroutine set_design_linear (ring)
!
! Subroutine to set only those elements on that constitute the "design"
! lattice. That is, only quadrupoles, bends and wigglers will be set on.
!
! Modules needed:
!   use bmad
!
! Output:
!   ring -- Ring_struct: Ring structure with quads, bends and wigglers only 
!           set on.
!-

#include "CESR_platform.inc"


subroutine set_design_linear (ring)

  use bmad_struct
  use bmad_interface, except => set_design_linear

  implicit none

  type (ring_struct) ring
  type (coord_struct), allocatable :: orb_(:)

  integer i, key

  allocate (orb_(0:ring%n_ele_maxx))

!

  do i = 1, ring%n_ele_use
    key = ring%ele_(i)%key
    if (key == quadrupole$ .or. key == sbend$ .or. key == wiggler$) then
      ring%ele_(i)%is_on = .true.
    else
      ring%ele_(i)%is_on = .false.
    endif
  enddo

  call ring_make_mat6 (ring, -1, orb_)

  deallocate(orb_)

end subroutine
