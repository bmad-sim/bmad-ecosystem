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

!$Id$
!$Log$
!Revision 1.3  2002/02/23 20:32:24  dcs
!Double/Single Real toggle added
!
!Revision 1.2  2001/09/27 18:31:57  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine set_design_linear (ring)

  use bmad

  implicit none

  type (ring_struct) ring
  type (coord_struct) orb_(0:n_ele_maxx)

  integer i, key

!

  do i = 1, ring%n_ele_ring
    key = ring%ele_(i)%key
    if (key == quadrupole$ .or. key == sbend$ .or. key == wiggler$) then
      ring%ele_(i)%is_on = .true.
    else
      ring%ele_(i)%is_on = .false.
    endif
  enddo

  call ring_make_mat6 (ring, -1, orb_)

end subroutine
