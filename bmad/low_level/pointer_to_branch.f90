!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function pointer_to_branch_given_name (branch_name, lat, parameter_is_branch0, blank_branch) result (branch_ptr)
!
! Function to point to the named lattice branch.
! This routine is overloaded by the routine: pointer_to_branch.
! See pointer_to_branch for more details.
! 
! Input:
!   branch_name          -- character(*): May be a branch name or a branch index.
!   lat                  -- lat_struct: Lattice to search.
!   parameter_is_branch0 -- logical, optional: If True, 'PARAMETER' is taken to be an alternative
!                             name for branch(0). Default is False. Used by Bmad parser.
!   blank_branch         -- integer, optional: Branch index if branch_name = ''. Default is blank is an error.
!
! Output:
!   branch_ptr  -- branch_struct, pointer: Pointer to the nameed branch.
!                   Nullified if there is no such branch.
!-

function pointer_to_branch_given_name (branch_name, lat, parameter_is_branch0, blank_branch) result (branch_ptr)

use bmad_struct

implicit none

type (branch_struct), pointer :: branch_ptr
type (lat_struct), target :: lat

integer ib, ios
integer, optional :: blank_branch
logical, optional :: parameter_is_branch0

character(*) branch_name
character(40) b_name
character(*), parameter :: r_name = 'pointer_to_branch_given_name'

! parameter

call str_upcase (b_name, branch_name)

if (logic_option(.false., parameter_is_branch0) .and. b_name == 'PARAMETER') then
  branch_ptr => lat%branch(0)
  return
endif

! Init in case of error

nullify(branch_ptr)

if (b_name == '') then
  if (.not. present(blank_branch)) return
  if (blank_branch < 0 .or. blank_branch > ubound(lat%branch, 1)) return
  branch_ptr => lat%branch(blank_branch)
  return
endif

! Is index.

if (is_integer(trim(b_name), ib)) then
  if (ib < 0 .or. ib > ubound(lat%branch, 1)) return
  branch_ptr => lat%branch(ib)
  return
endif

! Is name.

do ib = lbound(lat%branch, 1), ubound(lat%branch, 1)
  if (lat%branch(ib)%name == b_name) then
    branch_ptr => lat%branch(ib)
    return
  endif
enddo

end function pointer_to_branch_given_name

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function pointer_to_branch_given_ele (ele) result (branch_ptr)
!
! Function to point to the lattice branch associated with an element.
!
! This routine is overloaded by the routine: pointer_to_branch.
! See pointer_to_branch for more details.
!
! Note: Result is ambiguous if the element is associated with multiple branches which
! can happen, for example, with overlay lord elements.
!
! Input:
!   ele      -- Ele_struct: Element.
!
! Output:
!   branch_ptr  -- branch_struct, pointer: Pointer to the branch associated with the element.
!                   Nullified if the element is not associated with a lattice.
!-

recursive function pointer_to_branch_given_ele (ele) result (branch_ptr)

use bmad_struct

implicit none

type (ele_struct), target :: ele
type (ele_struct), pointer :: ele2
type (branch_struct), pointer :: branch_ptr
integer ixb

!

if (ele%ix_ele == ix_slice_slave$) then
  ele2 => ele%lord
else
  ele2 => ele
endif

! Not associated with a lattice case

if (.not. associated(ele2%branch)) then
  nullify(branch_ptr)
  return
endif

! A null_ele that is a lord is the result of superposition. This null_ele has no slaves and the
! original branch info is in %value(ix_branch$).

if (ele2%key == null_ele$ .and. ele2%ix_branch == 0 .and. ele2%ix_ele > ele2%branch%n_ele_track) then
  ixb = nint(ele2%value(ix_branch$))
  branch_ptr => ele2%branch%lat%branch(ixb)
  return
endif

! Not a lord case

if (ele2%n_slave == 0) then
  branch_ptr => ele2%branch
  return
endif

! If a lord then look to the first slave for the associated branch

branch_ptr => pointer_to_branch_given_ele(pointer_to_slave(ele2, 1))

end function pointer_to_branch_given_ele
