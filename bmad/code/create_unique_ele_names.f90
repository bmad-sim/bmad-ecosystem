!+
! Subroutine create_unique_ele_names (lat, key, suffix)
!
! Routine to give elements in a lattice unique names.
!
! This is done by appending the "suffix" argument to all elements that share 
! a common name. The "suffix" argument must have a single "?" charater in it. 
! When the suffix is applied, to the n^th element having a common name, 
! the number "n" is substituted for "?".
! 
! For example: If suffix is "_?" and the following elements are in the lattice:
!				QA    QB    QX    QA    QB     QB
! then after the suffix is applied the elements will have names:
!				QA_1  QB_1  QX    QA_2  QB_2   QB_3
! Notice that the suffix is not applied to any element with a unique name.
!
! The key argument is for restricting what elements can get their names modified. 
! For example, if key = quadrupole$ then only quadrupole elements will be looked at. 
! key = 0 means that all elements will be considered.
!
! Also see the routine: lat_ele_unique_name.
!
! Input:
!	  lat    -- Lat_struct: Lattice holding the elements.
!	  key    -- Integer: Class key of elements to consider.
!	  suffix -- Character(*): Suffix string. Must have a single "?" character.
!
! Output:
!		lat    -- Lat_struct: Lattice with names made unique.
!-

subroutine create_unique_ele_names (lat, key, suffix)

use bmad_interface, except => create_unique_ele_names

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele

integer key
integer i, j, ib, ix_p, ix, ixx, n_max
integer, allocatable :: name_indexx(:)

character(*) suffix
character(40) suff
character(40), allocatable :: original_name(:)
character(40) :: r_name = 'create_unique_ele_names'

! Find '?' character

ix_p = index(suffix, '?')
if (ix_p == 0) then
	call out_io (s_error$, r_name, 'SUFFIX DOES NOT HAVE A "?" CHARACTER: ' // suffix)
	return
endif

suff = suffix
call str_upcase (suff, suff)

!

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  branch%ele%ixx = 0

  ! Sort element names

  n_max = branch%n_ele_max
  call re_allocate (name_indexx, n_max)
  call re_allocate (original_name, n_max)
  original_name = branch%ele(1:n_max)%name
  call indexer (original_name, name_indexx)

  ! Find repeated names

  do i = 1, n_max
  	ele => branch%ele(i)

  	if (key /= 0 .and. ele%key /= key) cycle
  	call find_index (ele%name, original_name, name_indexx, n_max, ix, ixx)
  	if (ixx == n_max) cycle  ! Name is unique
  	j = name_indexx(ixx+1)
  	if (original_name(j) /= ele%name) cycle ! Name is unique

  	! Now add the suffix

  	branch%ele(ix)%ixx = branch%ele(ix)%ixx + 1
  	write (ele%name, '(2a, i0, a)') trim(ele%name), suff(1:ix_p-1), branch%ele(ix)%ixx, trim(suff(ix_p+1:))
  enddo
enddo

call create_lat_ele_nametable(lat, lat%nametable)

end subroutine
