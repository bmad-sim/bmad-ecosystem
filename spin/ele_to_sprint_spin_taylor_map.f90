!+
! Subroutine ele_to_sprint_spin_taylor_map (ele)
!
! Routine to calculate the spin Taylor map for a lattice element using the sprint formalism.
!
! Input:
!   ele       -- ele_struct: Element to form map for.
!
! Output:
!   ele       -- ele_struct: Element with map.
!     %spin_taylor(:)   -- Taylor map.
!-

subroutine ele_to_sprint_spin_taylor_map (ele)

use bmad, dummy => ele_to_sprint_spin_taylor_map

implicit none

type (ele_struct) ele

!

ele%spin_taylor_ref_orb_in = 0  ! Sprint ref is always the zero orbit

select case (ele%key)
case (drift$)

case (quadrupole$)

end select

end subroutine
