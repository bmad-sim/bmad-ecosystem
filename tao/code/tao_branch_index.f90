!+
! Function tao_branch_index (ix_branch) result (ix_this)
!
! Fnction to return the branch number.
! Generally ix_this = ix_branch except:
!   ix_this = -1  -> ix_this = s%global%default_branch.
!
! Input:
!   ix_branch           -- integer: Nominal branch number.
!
! Output:
!   ix_this      -- integer: Branch number. 
!-

function tao_branch_index (ix_branch) result (ix_this)

use tao_struct

implicit none

integer ix_branch, ix_this

!

select case (ix_branch)
case (-1)
  ix_this = s%global%default_branch

case default
  ix_this = ix_branch
end select

end function tao_branch_index

