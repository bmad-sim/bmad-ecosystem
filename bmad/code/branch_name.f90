!+
! Function branch_name(branch) result (name)
!
! Routine to return a string with the lattice branch name encoded.
! This routine is useful for error messages.
!
! Input:
!   branch    -- branch_struct: Lattice branch
!
! Output:
!   name      -- character(40): Encoded name
!-

function branch_name(branch) result (name)

use bmad_struct

implicit none

type (branch_struct), target :: branch
character(40) name

!

write (name, '(i0, 2a)') branch%ix_branch, ':', trim(branch%name)

end function branch_name
