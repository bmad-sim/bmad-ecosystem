!+
! Subroutine nametable_change1 (nametable, name, ix_name)
!
! Routine to change one entry in a nametable.
!
! The existing names with indexes in the range [ix_name:] will have their indexes increased by 1.
!
! Input:
!   nametable   -- nametable_struct: Nametable.
!   name        -- character(*): New name.
!   ix_name     -- integer: Name index where change is made.
!
! Output:
!   nametable   -- nametable_struct: Nametable with name changed
!-

subroutine nametable_change1 (nametable, name, ix_name)

use sim_utils, dummy => nametable_change1

implicit none

type (nametable_struct), target :: nametable
integer ix_name, i, n, ix_old, ix_new
integer, pointer :: n_min, n_max
character(*) name
character(*), parameter :: r_name = 'nametable_change1'

!

n_min => nametable%n_min
n_max => nametable%n_max

if (ix_name < n_min .or. ix_name > n_max) then
  call out_io (s_error$, r_name, 'INDEX OUT OF RANGE: \i0\ ', i_array = [ix_name])
  return
endif

ix_old = nametable_bracket_indexx(nametable, nametable%name(ix_name))
do while (nametable%index(ix_old) /= ix_name)
  ix_old = ix_old - 1
enddo

ix_new = nametable_bracket_indexx(nametable, name) + 1

!

nametable%name(ix_name) = name

if (ix_new > ix_old) then
  ix_new = ix_new - 1
  nametable%index(ix_old:ix_new-1) = nametable%index(ix_old+1:ix_new)
elseif (ix_new < ix_old) then
  nametable%index(ix_new+1:ix_old) = nametable%index(ix_new:ix_old-1)
endif

nametable%index(ix_new) = ix_name

end subroutine nametable_change1

