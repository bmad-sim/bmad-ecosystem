!+
! Subroutine insert_phase_trombone (branch)
!
! Routine to insert a match element in phase trombone mode as the first element in a lattice branch.
!
! Input:
!   branch      -- branch_struct: Lattice branch.
!
! Output:
!   branch      -- branch_struct: Lattice branch with trumbone at branch%ele(1).
!-

subroutine insert_phase_trombone (branch)

use bmad

implicit none

type (branch_struct), target :: branch
type (ele_struct) ele

!

call init_ele (ele, match$)
ele%value(phase_trombone$) = true$
call insert_element (branch%lat, ele, 1, branch%ix_branch)

end subroutine
