!+
! Subroutine nametable_add (nametable, name, ix_name)
!
! Routine to add a name to the nametable at index ix_name.
!
! The existing names with indexes in the range [ix_name:] will have their indexes increased by 1.
!
! Input:
!   nametable   -- nametable_struct: Nametable.
!   name        -- character(*): Name to add.
!   ix_name     -- integer: Index to add at.
!
! Output:
!   nametable   -- nametable_struct: Nametable with name added.
!-

subroutine nametable_add (nametable, name, ix_name)

use sim_utils, dummy => nametable_add

implicit none

type (nametable_struct), target :: nametable
integer ix_name, i, n, ix_id
integer, pointer :: n_min, n_max
character(*) name
character(*), parameter :: r_name = 'nametable_add'

!

n_min => nametable%n_min
n_max => nametable%n_max

if (ix_name < n_min .or. ix_name > n_max + 1) then
  call out_io (s_error$, r_name, 'INDEX OUT OF RANGE: \i0\ ', i_array = [ix_name])
  return
endif

ix_id = nametable_bracket_indexx(nametable, name)

n_max = n_max + 1
if (n_max > ubound(nametable%name, 1)) then
  n = n_max - n_min + 10
  call re_allocate2(nametable%index, n_min, n_max+n, .false.)
  call re_allocate2(nametable%name, n_min, n_max+n, .false.)
endif

!

if (n_max > ix_name) nametable%name(ix_name+1:n_max) = nametable%name(ix_name:n_max-1)
nametable%name(ix_name) = name


do i = n_min, n_max-1
  if (nametable%index(i) >= ix_name) nametable%index(i) = nametable%index(i) + 1
enddo

if (n_max > ix_id+1) nametable%index(ix_id+2:n_max) = nametable%index(ix_id+1:n_max-1)
nametable%index(ix_id+1) = ix_name

end subroutine nametable_add

