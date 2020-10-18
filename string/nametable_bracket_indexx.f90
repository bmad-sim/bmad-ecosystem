!+
! Function nametable_bracket_indexx (nametable, name) result (ix_indexx)
!
! Routine to find the index ix_indexx such that:
!   nametable%name(nametable%indexx(ix_indexx)) <= name < nametable%name(nametable%indexx(ix_indexx+1))
!
! Input:
!   nametable   -- nametable_struct: Nametable.
!   name        -- character(*): Name to match against.
!
! Output:
!   ix_indexx   -- integer: Index.
!-

function nametable_bracket_indexx (nametable, name) result (ix_indexx)

use sim_utils, dummy => nametable_bracket_indexx

implicit none

type (nametable_struct) nametable
character(40) name
integer ix_indexx
integer ix1, ix3

!

ix1 = nametable%n_min
ix3 = nametable%n_max

if (ix1 > ix3) then
  ix_indexx = ix1 - 1
  return
endif

!

do
  if (ix1 == ix3) then
    if (nametable%name(nametable%indexx(ix1)) < name) then
      ix_indexx = ix1
    elseif (nametable%name(nametable%indexx(ix1)) == name) then
      ix_indexx = ix1
    else
      ix_indexx = ix1 - 1
    endif
    return
  endif

  if (ix1 + 1 == ix3) then
    if (name < nametable%name(nametable%indexx(ix1))) then
      ix_indexx = ix1 - 1
    elseif (name >= nametable%name(nametable%indexx(ix3))) then
      ix_indexx = ix3
    else
      ix_indexx = ix1
    endif
    return
  endif

  ix_indexx = (ix1 + ix3) / 2 

  if (nametable%name(nametable%indexx(ix_indexx)) == name) then
    do ! if there are duplicate nametable%name in the list choose the last one
      if (ix_indexx == nametable%n_max) exit
      if (nametable%name(nametable%indexx(ix_indexx+1)) /= nametable%name(nametable%indexx(ix_indexx))) exit
      ix_indexx = ix_indexx + 1
    enddo
    return

  elseif (nametable%name(nametable%indexx(ix_indexx)) < name) then
    ix1 = ix_indexx + 1

  else
    ix3 = ix_indexx - 1
  endif
enddo

end function nametable_bracket_indexx
