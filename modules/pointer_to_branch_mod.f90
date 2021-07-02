module pointer_to_branch_mod

use bmad_struct

implicit none

private pointer_to_branch_given_name, pointer_to_branch_given_ele

!+
! Function pointer_to_branch
!
! Routine to return a pointer to the lattice branch associated with a given name
! or a given element.
!
! This routine is an overloaded name for:
!   pointer_to_branch_given_ele (ele) result (branch_ptr, parameter_is_branch0, blank_is_branch0)
!   pointer_to_branch_given_name (branch_name, lat) result (branch_ptr)
!
! The lattice branch *associated* with a given element is not necessarily the
! branch where the element is *located*. For example, all lords live in branch #0.
! But the branch associated with a super_lord element is the branch of its slaves.
!
! To get the branch where the element is located, simply use ele%ix_branch.
! 
! Note: Result is ambiguous if ele argument is associated with multiple branches 
! which can happen, for example, with overlay elements.
!
! Input:
!   ele                  -- Ele_struct: Element contained in the branch.
!   branch_name          -- Character(*): May be a branch name or a branch index.
!   lat                  -- Lat_struct: Lattice to search.
!   parameter_is_branch0 -- logical, optional: If True, 'PARAMETER' is taken to be
!                             an alternative name for branch(0). Default is False.
!   blank_is_branch0     -- logical, optional: If True, branch_name = '' is taken to be
!                             an alternative name for branch(0). Default is False.
!
! Output:
!   branch_ptr  -- branch_struct, pointer: Pointer to the branch.
!                   Nullified if there is no associated branch.
!-

interface pointer_to_branch
  module procedure pointer_to_branch_given_ele
  module procedure pointer_to_branch_given_name
end interface

contains

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function pointer_to_branch_given_name (branch_name, lat, parameter_is_branch0, blank_is_branch0) result (branch_ptr)
!
! Function to point to the named lattice branch.
! This routine is overloaded by the routine: pointer_to_branch.
! See pointer_to_branch for more details.
! 
! Input:
!   branch_name -- Character(*): May be a branch name or a branch index.
!   lat         -- Lat_struct: Lattice to search.
!   parameter_is_branch0 -- logical, optional: If True, 'PARAMETER' is taken to be
!                     an alternative name for branch(0). Default is False.
!
! Output:
!   branch_ptr  -- branch_struct, pointer: Pointer to the nameed branch.
!                   Nullified if there is no such branch.
!-

function pointer_to_branch_given_name (branch_name, lat, parameter_is_branch0, blank_is_branch0) result (branch_ptr)

type (branch_struct), pointer :: branch_ptr
type (lat_struct), target :: lat

integer ib, ios
logical, optional :: parameter_is_branch0, blank_is_branch0

character(*) branch_name
character(40) b_name
character(*), parameter :: r_name = 'pointer_to_branch_given_name'

! parameter

call str_upcase (b_name, branch_name)

if (logic_option(.false., parameter_is_branch0) .and. b_name == 'PARAMETER' .or. &
    logic_option(.false., blank_is_branch0) .and. b_name == '') then
  branch_ptr => lat%branch(0)
  return
endif

! Init in case of error

nullify(branch_ptr)

! Is index.

if (is_integer(trim(b_name))) then
  read (b_name, *, iostat = ios) ib
  if (ib < 0 .or. ib > ubound(lat%branch, 1)) then
    !! call out_io (s_error$, r_name, 'BRANCH INDEX OUT OF RANGE: ' // b_name)
    return
  endif
  branch_ptr => lat%branch(ib)

! Is name.

else
  do ib = lbound(lat%branch, 1), ubound(lat%branch, 1)
    if (lat%branch(ib)%name == b_name) then
      branch_ptr => lat%branch(ib)
      return
    endif
  enddo
  !! call out_io (s_error$, r_name, 'BRANCH NAME NOT FOUND: ' // branch_name)
  return
endif

end function pointer_to_branch_given_name

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function pointer_to_branch_given_ele (ele) result (branch_ptr)
!
! Function to point to the lattice branch associated with an element.
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

type (ele_struct), target :: ele
type (branch_struct), pointer :: branch_ptr
integer ixb

! Not associated with a lattice case

if (.not. associated(ele%branch)) then
  nullify(branch_ptr)
  return
endif

! A null_ele that is a lord is the result of superposition. This null_ele has no slaves and the
! original branch info is in %value(ix_branch$).

if (ele%key == null_ele$ .and. ele%ix_branch == 0 .and. ele%ix_ele > ele%branch%n_ele_track) then
  ixb = nint(ele%value(ix_branch$))
  branch_ptr => ele%branch%lat%branch(ixb)
  return
endif

! Not a lord case

if (ele%n_slave == 0) then
  branch_ptr => ele%branch
  return
endif

! If a lord then look to the first slave for the associated branch

branch_ptr => pointer_to_branch_given_ele(pointer_to_slave(ele, 1))

end function pointer_to_branch_given_ele

end module
