!+
! Subroutine s_calc (ring)
!
! Subroutine to calculate the longitudinal distance S for the elements
! in a ring.
!
! Modules Needed:
!   use bmad
!
! Input:
!   ring -- Ring_struct:
!
! Output:
!   ring -- Ring_struct:
!-

!$Id$
!$Log$
!Revision 1.7  2003/03/18 20:33:04  dcs
!%total_length calc subtracts off %ele_(0)%s
!
!Revision 1.6  2003/03/08 01:31:22  dcs
!Allow fixing initial s position
!
!Revision 1.4  2002/02/23 20:32:24  dcs
!Double/Single Real toggle added
!
!Revision 1.3  2002/01/08 21:44:43  dcs
!Aligned with VMS version  -- DCS
!
!Revision 1.2  2001/09/27 18:31:57  rwh24
!UNIX compatibility updates

#include "CESR_platform.inc"

subroutine s_calc (ring)

  use bmad_struct
  use bmad_interface
  
  implicit none

  type (ring_struct)  ring

  integer n, ix2
  real*8 ss

! Just go through all the elements and add up the lengths.

  ss = ring%ele_(0)%s

  do n = 1, ring%n_ele_ring
    ss = ring%ele_(n-1)%s + ring%ele_(n)%value(l$)
    ring%ele_(n)%s = ss
  enddo

  ring%param%total_length = ss - ring%ele_(0)%s

! now get fill in the positions of the super_lords

  do n = ring%n_ele_ring+1, ring%n_ele_max
    if (ring%ele_(n)%control_type == super_lord$) then
      ix2 = ring%control_(ring%ele_(n)%ix2_slave)%ix_slave
      ring%ele_(n)%s = ring%ele_(ix2)%s
    else
      ring%ele_(n)%s = 0
    endif
  enddo

end subroutine
