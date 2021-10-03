!+
! Subroutine ele_to_spin_sprint (ele, orb0)
!
! Routine to calculate the spin map for a lattice element using the sprint expansion.
!
! Input:
!   ele       -- ele_struct: Element to form map for.
!   orb0      -- coord_struct, optional: Coordinates about which map is made.
!
! Output:
!   ele       -- ele_struct: Element with map.
!     %spin_taylor(:)   -- Taylor map.
!     %spin_q           -- Map to 1st order.
!-

subroutine ele_to_spin_sprint (ele, orb0)

use bmad, dummy => ele_to_spin_sprint

implicit none

type (ele_struct) ele
type (coord_struct), optional :: orb0


end subroutine
