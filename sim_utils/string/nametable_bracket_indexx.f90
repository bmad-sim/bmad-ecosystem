!+
! Function nametable_bracket_indexx (nametable, name, n_match) result (ix_max)
!
! Routine to find the index ix_max such that:
!   nametable%name(nametable%index(ix_max)) <= name < nametable%name(nametable%index(ix_max+1))
!
! Input:
!   nametable   -- nametable_struct: Nametable.
!   name        -- character(*): Name to match against. Names are limited to 40 characters.
!
! Output:
!   n_match     -- integer, optional: Number of matches.
!   ix_max      -- integer: Upper bound of match region
!-

function nametable_bracket_indexx (nametable, name, n_match) result (ix_max)

use sim_utils, dummy => nametable_bracket_indexx

implicit none

type (nametable_struct) nametable
character(*) name
integer, optional :: n_match
integer ix_max
integer ix1, ix3

!

ix1 = nametable%n_min
ix3 = nametable%n_max

if (ix1 > ix3) then
  ix_max = ix1 - 1
  if (present(n_match)) n_match = 0
  return
endif

!

do
  if (ix1 == ix3) then
    if (nametable%name(nametable%index(ix1)) < name) then
      ix_max = ix1
    elseif (nametable%name(nametable%index(ix1)) == name) then
      ix_max = ix1
    else
      ix_max = ix1 - 1
    endif
    if (present(n_match)) n_match = calc_n_match(nametable, name, ix_max)
    return
  endif

  if (ix1 + 1 == ix3) then
    if (name < nametable%name(nametable%index(ix1))) then
      ix_max = ix1 - 1
    elseif (name >= nametable%name(nametable%index(ix3))) then
      ix_max = ix3
    else
      ix_max = ix1
    endif
    if (present(n_match)) n_match = calc_n_match(nametable, name, ix_max)
    return
  endif

  ix_max = (ix1 + ix3) / 2 

  if (nametable%name(nametable%index(ix_max)) == name) then
    do ! if there are duplicate nametable%name in the list choose the last one
      if (ix_max == nametable%n_max) exit
      if (nametable%name(nametable%index(ix_max+1)) /= nametable%name(nametable%index(ix_max))) exit
      ix_max = ix_max + 1
    enddo

    if (present(n_match)) n_match = calc_n_match(nametable, name, ix_max)
    return

  elseif (nametable%name(nametable%index(ix_max)) < name) then
    ix1 = ix_max + 1

  else
    ix3 = ix_max - 1
  endif
enddo

!------------------------------------------
contains

function calc_n_match(nametable, name, ix_max) result (n_match)

type (nametable_struct) nametable
character(40) name
integer :: n_match
integer i, ix_max

!

do i = ix_max, nametable%n_min, -1
  if (nametable%name(nametable%index(i)) /= name) exit
enddo

n_match = ix_max - i

end function calc_n_match

end function nametable_bracket_indexx
