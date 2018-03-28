!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine transfer_branches (branch1, branch2)
!
! Subroutine to set branch2 = branch1. 
! This is a plain transfer of information not using the overloaded equal sign.
! Thus, at the end, branch2's pointers point to the same memory as branch1's.
!
! NOTE: Do not use this routine unless you know what you are doing!
!
! Input:
!   branch1(:) -- Branch_struct:
!
! Output:
!   branch2(:) -- Branch_struct:
!-

subroutine transfer_branches (branch1, branch2)

use bmad_struct

implicit none

type (branch_struct) :: branch1(:)
type (branch_struct) :: branch2(:)

branch2 = branch1

end subroutine transfer_branches

