!+
! Subroutine tao_pointer_to_branches (branch_str, branches, unis, err)
!
! Routine to point branches to use.
!
! Examples:
!   ""          -- Default branch in default universe.
!   "[1,3]"     -- Default branch in universes 1 and 3.
!   "2"         -- Branch #2 in default universe.
!   "*@1"       -- #1 branch in all universes.
!
! Input:
!   branch_str      -- character(*): String specifying what branches to use.
!
! Output:
!   branches(:)     -- branch_pointer_struct, allocatable: Array of pointers to branches.
!   unis(:)         -- tao_universe_pointer_struct, allocatable: Array of corresponding universes.
!   err             -- logical: Set True if there is an error.
!-

subroutine tao_pointer_to_branches (branch_str, branches, unis, err)

use tao_interface, dummy => tao_pointer_to_branches

implicit none

type (branch_pointer_struct), allocatable :: branches(:)
type (tao_universe_pointer_struct), allocatable, target :: unis(:)
type (tao_universe_struct), pointer :: u

integer i
logical err
character(40) b_str
character(*) branch_str
character(*), parameter :: r_name = 'tao_pointer_to_branches'

!

if (allocated(branches)) deallocate(branches)
call tao_pointer_to_universes (branch_str, unis, err, b_str); if (err) return
allocate(branches(size(unis)))

do i = 1, size(unis)
  u => unis(i)%u
  if (b_str == '') then
    branches(i)%branch => u%model%lat%branch(s%global%default_branch)
  else
    branches(i)%branch => pointer_to_branch(b_str, u%model%lat)
  endif
enddo

end subroutine tao_pointer_to_branches
