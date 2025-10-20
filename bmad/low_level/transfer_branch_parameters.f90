!+
! Subroutine transfer_branch_parameters (branch_in, branch_out)
!
! Subroutine to transfer the branch parameters (such as branch%name, branch%param, etc.) from 
! one branch to another. Not touched are arrays and PTC related. Specifically not touched is:
!   branch%ele(:)
!   branch%wall3d(:)
!   branch%ptc
!
! Input:
!   branch_in   -- branch_struct: Input branch.
!
! Output:
!   branch_out  -- branch_struct: Output branch with parameters set.
!-

subroutine transfer_branch_parameters (branch_in, branch_out)

use bmad_struct

implicit none

type (branch_struct), intent(in) :: branch_in
type (branch_struct) :: branch_out

! Note: branch_out inherits %ix_branch from branch_in.
! In most cases this is fine but there are exceptions.

branch_out%name           = branch_in%name
branch_out%ix_branch      = branch_in%ix_branch
branch_out%ix_from_branch = branch_in%ix_from_branch
branch_out%ix_from_ele    = branch_in%ix_from_ele
branch_out%ix_to_ele      = branch_in%ix_to_ele
branch_out%ix_fixer       = branch_in%ix_fixer
branch_out%n_ele_track    = branch_in%n_ele_track
branch_out%n_ele_max      = branch_in%n_ele_max
branch_out%param          = branch_in%param
branch_out%a              = branch_in%a
branch_out%b              = branch_in%b
branch_out%z              = branch_in%z
branch_out%ele%ix_branch  = branch_in%ix_branch

end subroutine transfer_branch_parameters

