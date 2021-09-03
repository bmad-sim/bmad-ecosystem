!+
! Subroutine nametable_remove (nametable, ix_name)
!
! Routine to remove a name from the nametable at index ix_name.
!
! The existing names with indexes in the range [ix_name:] will have their indexes decreased by 1.
!
! Input:
!   nametable   -- nametable_struct: Nametable.
!   ix_name     -- integer: Index of name to remove.
!
! Output:
!   nametable   -- nametable_struct: Nametable with name removeed.
!-

subroutine nametable_remove (nametable, ix_name)

use output_mod, dummy => nametable_remove

implicit none

type (nametable_struct), target :: nametable
integer ix_name, i0, i, ix_id
integer, pointer :: n_min, n_max
character(*), parameter :: r_name = 'nametable_remove'

!

n_min => nametable%n_min
n_max => nametable%n_max

if (ix_name < n_min .or. ix_name > n_max) then
  call out_io (s_error$, r_name, 'INDEX OUT OF RANGE: \i0\ ', i_array = [ix_name])
  return
endif

ix_id = nametable_bracket_indexx(nametable, nametable%name(ix_name))

!

n_max = n_max - 1
nametable%name(ix_name:n_max) = [nametable%name(ix_name+1:n_max+1)]

do i0 = ix_id, n_min, -1
  if (nametable%index(i0) == ix_name) exit
enddo

do i = n_min, n_max+1
  if (nametable%index(i) >= ix_name) nametable%index(i) = nametable%index(i) - 1
enddo

nametable%index(i0:n_max) = nametable%index(i0+1:n_max+1)

end subroutine nametable_remove
