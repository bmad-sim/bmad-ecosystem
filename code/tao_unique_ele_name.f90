!+
! Function tao_unique_ele_name (ele, nametable) result (unique_name)
!
! Routine to create a unique name for a lattice element.
! For example "Q3##2" would be a unique name if there were multiple Q3 elements.
!
! The idea is to return a name that has the best chance of being "correct" even if the
! lattice is modified. For example, just using the element index as the name will break if
! elements are added to the lattice.
!
! Input:
!   ele           -- ele_struct: Element under consideration.
!   nametable     -- lat_nametable_struct: Lattice nametable.
!
! Output:
!   unique_name   -- character(40): Unique name to identify element. 
!                     Set to "???" if ele_name is not in the nametable.
!-

function tao_unique_ele_name (ele, nametable) result (unique_name)

use bmad_struct

implicit none

type (ele_struct), target :: ele
type (ele_struct), pointer :: ele0
type (lat_nametable_struct) nametable

integer i, j, ix, ix2, n_match, n, nc, n_in

character(40) unique_name
character(*), parameter :: r_name = 'tao_unique_ele_name'

!

unique_name = '???'
call find_indexx (ele%name, nametable%name, nametable%indexx, size(nametable%name), ix, ix2, n_match = n_match)
if (n_match == 0) then 
  call out_io (s_error$, r_name, 'CANNOT MATCH LATTICE ELEMENT NAME TO NAMELIST!' // ele%name)
  return
endif

! Unique name?

if (n_match == 1) then
  unique_name = ele%name
  return
endif

! If element name is unique for the branch then use "ix_branch>>" prefix.

n_in = 0
do n = ix2, ix2+n_match-1
  j = nametable%indexx(n)
  ele0 => nametable%ele(j)%ele
  if (ele0%ix_branch == ele%ix_branch) n_in = n_in + 1
enddo

if (n_in == 1) then
  write (unique_name, '(i0, 2a)') ele%ix_branch, '>>', ele%name
  return
endif

! When all else fails, use the "##" construct.

nc = 0
do n = ix2, ix2+n_match-1
  j = nametable%indexx(n)
  ele0 => nametable%ele(j)%ele
  if (ele0%ix_branch == ele%ix_branch) nc = nc + 1
  if (.not. associated (ele0, ele)) cycle
  write (unique_name, '(2a, i0)') trim(ele%name), '##', nc
  exit
enddo

! If there are elements of the same name in other branches then use a branch prefix.

if (n_in /= n_match .and. size(ele%branch%lat%branch) /= 1) then
  write (unique_name, '(i0, 2a)') ele%ix_branch, '>>', unique_name
endif

end function tao_unique_ele_name
