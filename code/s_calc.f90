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

#include "CESR_platform.inc"

subroutine s_calc (ring)

  use bmad_struct
  use bmad_interface, except => s_calc
  
  implicit none

  type (ring_struct)  ring

  integer n, ix2
  real*8 ss

! Just go through all the elements and add up the lengths.

  ss = ring%ele_(0)%s

  do n = 1, ring%n_ele_use
    ss = ring%ele_(n-1)%s + ring%ele_(n)%value(l$)
    ring%ele_(n)%s = ss
  enddo

  ring%param%total_length = ss - ring%ele_(0)%s

! now get fill in the positions of the super_lords

  do n = ring%n_ele_use+1, ring%n_ele_max
    if (ring%ele_(n)%control_type == super_lord$) then
      ix2 = ring%control_(ring%ele_(n)%ix2_slave)%ix_slave
      ring%ele_(n)%s = ring%ele_(ix2)%s
    else
      ring%ele_(n)%s = 0
    endif
  enddo

end subroutine
