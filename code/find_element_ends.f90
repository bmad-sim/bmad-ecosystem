!+
! Subroutine find_element_ends (ring, ix_ele, ix_start, ix_end)
!
! Subroutine to find the end points of an element in the regular part of the 
! ring.
!
! Modules Needed:
!   use bmad
!
! Input:
!   ring   -- Ring_struct: Ring holding the lattice
!   ix_ele -- Integer: Index of element to find the ends for.
!
! Output:
!   ix_start -- Integer: Index of the element just before the element.
!   ix_end   -- Integer: Index of element itself or index of the
!                   last sub-element within the element.
!
! Note: For an element in the regular part of the ring:
!       ix_start = ix_ele - 1
!       ix_end = ix_ele
!-

!$Id$
!$Log$
!Revision 1.4  2002/02/23 20:32:15  dcs
!Double/Single Real toggle added
!
!Revision 1.3  2002/01/08 21:44:39  dcs
!Aligned with VMS version  -- DCS
!
!Revision 1.2  2001/09/27 18:31:51  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine find_element_ends (ring, ix_ele, ix_start, ix_end)

  use bmad
  use nr, only: indexx

  implicit none
                                                         
  type (ring_struct) ring

  integer ix_ele, ix_start, ix_end
  integer ix1, ix2, n, ix_(100), indx(100)

!

  ix1 = ring%ele_(ix_ele)%ix1_slave
  ix2 = ring%ele_(ix_ele)%ix2_slave

  if (ring%ele_(ix_ele)%n_slave == 0) then
    ix_start = ix_ele - 1
    ix_end = ix_ele
  elseif (ring%ele_(ix_ele)%control_type == super_lord$) then
    ix_start = ring%control_(ix1)%ix_slave - 1
    ix_end = ring%control_(ix2)%ix_slave
  else  ! overlay_lord$ or group_lord$
    n = ring%ele_(ix_ele)%n_slave
    ix_(1:n) = ring%control_(ix1:ix2)%ix_slave
    ix_(2:n) = ix_(2:n) + ring%n_ele_ring * &
                             nint(float(ix_(1) - ix_(2:n)) / ring%n_ele_ring)
    call indexx (ix_(1:n), indx(1:n))
    ix_start = ix_(indx(1)) - 1
    ix_end = ix_(indx(n))
  endif

end subroutine
