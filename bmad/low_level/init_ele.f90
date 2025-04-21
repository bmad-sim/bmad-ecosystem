!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine init_ele (ele, key, sub_key, ix_ele, branch)
!
! Subroutine to initialize a Bmad element.
!
! Input:
!   key       -- Integer, optional: Key to initialize to. EG: quadrupole$, etc.
!   sub_key   -- Integer, optional: Sub-key to initialize to.
!   ix_ele    -- Integer, optional: ix_ele index to initalize to. Default = -1.
!   branch    -- branch_struct, optional: Branch to point ele%branch and ele%ix_branch to.
!                  Otherwise ele%branch is nullified and ele%ix_branch = 0.
!
! Output:
!   ele -- Ele_struct: Initialized element.
!-

subroutine init_ele (ele, key, sub_key, ix_ele, branch)

use bmad_routine_interface, dummy => init_ele

implicit none

type (ele_struct)  ele
type (branch_struct), optional, target :: branch
integer, optional :: key, sub_key
integer, optional :: ix_ele

!

call deallocate_ele_pointers (ele)

ele = ele_struct()

if (present(branch)) then
  ele%branch => branch
  ele%ix_branch = branch%ix_branch
endif

if (present(ix_ele)) ele%ix_ele = ix_ele

if (present(key)) then
  ele%key = key
  call set_ele_defaults(ele)
endif

if (present(sub_key)) ele%sub_key = sub_key

end subroutine init_ele

