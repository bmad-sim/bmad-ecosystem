!+
! Subroutine NEW_CONTROL (RING, IX_ELE)
!
! Subroutine to create a new control element.
!
! Modules Needed:
!   use bmad
!
! Input:
!     RING -- Ring_struct: Ring used
!
! Output
!     IX_ELE -- Integer: Index of the new control element
!-

#include "CESR_platform.inc"

subroutine new_control (ring, ix_ele)

  use bmad_struct
  use bmad_interface, except => new_control

  implicit none

  type (ring_struct)  ring
  integer ix_ele

!

  ring%n_ele_max = ring%n_ele_max + 1
  ix_ele = ring%n_ele_max

  if (ix_ele > ubound(ring%ele_, 1))  call allocate_ring_ele_(ring)
  call init_ele (ring%ele_(ix_ele))

end subroutine
