module indexer_mod

use re_allocate_mod
use output_mod
use swap_mod

implicit none

! Note: To initialize use the routine init_str_index_struct

interface indexer
  module procedure indexer_str
  module procedure indexer_var_str
  module procedure indexer_real
  module procedure indexer_int
end interface

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine find_index(...)
!
! This routine is an overloaded name for: 
!   Subroutine find_index0 (name, var_str_index, ix_match, ix2_match, add_to_list, n_match, has_been_added)
!   Subroutine find_index1 (name, nametable, ix_match, ix2_match, add_to_list, n_match, has_been_added)
!   Subroutine find_index2 (name, names, an_index, n_max, ix_match, ix2_match, add_to_list, n_match, had_been_added)
!   Subroutine find_index3 (name, names, n_min, an_index, n_max, ix_match, ix2_match, add_to_list, n_match, has_been_added)
!
! Subroutine to find a matching name in a list of names.
! The routine indexer should be used to create an_index.
! Also see: find_index2
!
! Note: The difference between the nametable_struct and str_index_struct is that the nametable_struct stores 40 char strings
! and the str_index_struct stores var length strings.
!
! Note: names and an_index arguments are not allocatable so it is up to the calling routine to make sure that their size is large
! enough if add_to_list = True.
!
! If add_to_list = True and if name is not in the list: Name is added to the names_list and
! an_index is updated using the prescription:
!   First: Find ix2_match.
!   an_index(ix2_match+1:n_max+1) = an_index(ix2_match:n_max)
!   an_index(ix2_match) = n_max + 1
!   names(n_max+1) = name
!   ix_match = n_max+1
!
! Input:
!   name           -- Character(*): Name to match to.
!   nametable      -- nametable_struct: Structure containing array of names and index array.
!   var_str_index  -- str_index_struct: Structure containing array of names and index array.
!   n_min          -- integer: Min index of names(:). For find_index0 and find_index1 effectively n_min = 1.
!   names(:)       -- Character(*): Array of names.
!   an_index(:)    -- Integer: Sorted index for names(:) array. names(an_index(i)) is in alphabetical order.
!   n_max          -- Integer: Use only names(n_min:n_max) part of array.
!   add_to_list    -- Logical, optional: If present and True and name does not appear in the names(:) array, 
!                         add name to the end of the names(:) array and update the str_indexer%index array.
!                         Default is False.
!  
! Output:
!   n_max          -- Integer: Increased by 1 if name is added to names(:) array.
!   ix_match       -- Integer: If a match is found then: names(ix_match) = name
!                              If no match is found then: ix_match = n_min - 1.
!                              Note: There will always be a match if add_to_list = T
!   ix2_match      -- Integer, optional: 
!                       If a match is found then
!                                an_index(ix2_match) = ix_match
!                                names(an_index(ix2_match-1)) /= name
!                       If no match is found then 
!                       For j = an_index(ix2_match): names(j) > name
!                       And if ix2_match > 1 then for j = an_index(ix2_match-1): names(j) < name
!   names(:)       -- Character(*): Updated if add_to_list = True.
!   var_str_index  -- str_index_struct: Updated if add_to_list = True.
!   nametable      -- nametable_struct: Structure containing array of names and index array.
!   an_index(:)    -- Integer: Updated if add_to_list = True.
!   n_match        -- integer, optional: Number items in the updated names list that match .
!   has_been_added -- Logical, optional: Set True if name has been added to the list. False otherwise.
!-

interface find_index
  module procedure find_index0
  module procedure find_index1
  module procedure find_index2
  module procedure find_index3
end interface

contains

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine init_str_index_struct(var_str_indexer, n_min)
! 
! Routine to initialize a str_index_struct instance.
!
! Input:
!   n_min   -- integer, optional: Lower bound of array to be indexed. Default is 1.
!
! Output:
!   var_str_indexer -- str_index_struct: Variable to be initialized.
!-

subroutine init_str_index_struct (var_str_indexer, n_min)

type (str_index_struct) var_str_indexer
integer, optional :: n_min

!

var_str_indexer%n_min = integer_option(1, n_min)
var_str_indexer%n_max = var_str_indexer%n_min - 1

end subroutine init_str_index_struct

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine find_index0 (name, var_str_index, ix_match, ix2_match, add_to_list, n_match, has_been_added)
!
! Subroutine to find a matching name in a list of names.
! This name is overloaded by find_index.
! See the find_index documentation for more details.
!-

subroutine find_index0 (name, var_str_index, ix_match, ix2_match, add_to_list, n_match, has_been_added)

character(*) name
type (str_index_struct), target :: var_str_index
integer  ix_match
integer, optional :: ix2_match
integer, optional :: n_match
integer ix1, ix2, ix3, n
integer, pointer :: n_min, n_max
logical, optional :: add_to_list, has_been_added

! simple case where there are no elements

if (present(has_been_added)) has_been_added = .false.

n_min => var_str_index%n_min
n_max => var_str_index%n_max

if (n_max < n_min) then
  if (present(ix2_match)) ix2_match = n_min
  ix_match = n_min - 1
  if (logic_option(.false., add_to_list)) then
    call re_allocate2(var_str_index%name, n_min, n_min+10)
    call re_allocate2(var_str_index%index, n_min, n_min+10)
    ix_match = n_min
    var_str_index%index(n_min) = n_min
    var_str_index%name(n_min)%str = name
    n_max = n_min
    if (present(n_match)) n_match = 1
    if (present(has_been_added)) has_been_added = .true.
  else
    if (present(n_match)) n_match = 0
  endif
  return
endif

!

ix1 = n_min
ix3 = n_max

do
  ix2 = (ix1 + ix3) / 2 

  if (var_str_index%name(var_str_index%index(ix2))%str == name) then
    do ! if there are duplicate name in the list choose the first one
      if (ix2 == n_min) exit
      if (var_str_index%name(var_str_index%index(ix2-1))%str /= var_str_index%name(var_str_index%index(ix2))%str) exit
      ix2 = ix2 - 1
    enddo
    ix_match = var_str_index%index(ix2)
    exit
  elseif (var_str_index%name(var_str_index%index(ix2))%str < name) then
    ix1 = ix2 + 1
  else
    ix3 = ix2 - 1
  endif
                     
  if (ix1 > ix3) then
    ix_match = n_min - 1
    if (var_str_index%name(var_str_index%index(ix2))%str < name) ix2 = ix2 + 1
    exit
  endif
enddo

if (present(ix2_match)) ix2_match = ix2

if (logic_option(.false., add_to_list) .and. ix_match < n_min) then
  if (n_max + 1 > ubound(var_str_index%name, 1)) then
    n = n_max - n_min
    call re_allocate2(var_str_index%name, n_min, n_max+n)
    call re_allocate2(var_str_index%index, n_min, n_max+n)
  endif
  if (present(has_been_added)) has_been_added = .true.
  var_str_index%index(ix2+1:n_max+1) = var_str_index%index(ix2:n_max)
  var_str_index%index(ix2) = n_max + 1
  var_str_index%name(n_max+1)%str = name
  ix_match = n_max + 1
  n_max = n_max+1
endif

if (present(n_match)) then
  n_match = 0
  if (ix_match == n_min - 1) return
  do
    if (var_str_index%name(var_str_index%index(ix2))%str /= name) return
    n_match = n_match +1
    ix2 = ix2 + 1
    if (ix2 > n_max) return
  enddo
endif

end subroutine find_index0

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine find_index1 (name, nametable, ix_match, ix2_match, add_to_list, n_match, has_been_added)
!
! Subroutine to find a matching name in a list of names.
! This name is overloaded by find_index.
! See the find_index documentation for more details.
!-

subroutine find_index1 (name, nametable, ix_match, ix2_match, add_to_list, n_match, has_been_added)

type (nametable_struct), target :: nametable

integer n, ix1, ix2, ix3, ix_match
integer, optional :: ix2_match
integer, pointer :: n_min, n_max
integer, optional :: n_match

character(*) name
logical, optional :: add_to_list, has_been_added

! simple case where there are no elements

n_min => nametable%n_min
n_max => nametable%n_max

if (present(has_been_added)) has_been_added = .false.

if (n_max < n_min) then
  if (present(ix2_match)) ix2_match = n_min
  ix_match = n_min - 1
  if (logic_option(.false., add_to_list)) then
    call re_allocate2(nametable%name, n_min, n_min+10)
    call re_allocate2(nametable%index, n_min, n_min+10)
    ix_match = n_min
    nametable%index(n_min) = n_min
    nametable%name(n_min) = name
    n_max = n_min
    if (present(n_match)) n_match = 1
    if (present(has_been_added)) has_been_added = .true.
  else
    if (present(n_match)) n_match = 0
  endif
  return
endif

!

ix1 = n_min
ix3 = n_max

do
  ix2 = (ix1 + ix3) / 2 

  if (nametable%name(nametable%index(ix2)) == name) then
    do ! if there are duplicate nametable%name in the list choose the first one
      if (ix2 == n_min) exit
      if (nametable%name(nametable%index(ix2-1)) /= nametable%name(nametable%index(ix2))) exit
      ix2 = ix2 - 1
    enddo
    ix_match = nametable%index(ix2)
    exit
  elseif (nametable%name(nametable%index(ix2)) < name) then
    ix1 = ix2 + 1
  else
    ix3 = ix2 - 1
  endif
                     
  if (ix1 > ix3) then
    ix_match = n_min - 1
    if (nametable%name(nametable%index(ix2)) < name) ix2 = ix2 + 1
    exit
  endif
enddo

if (present(ix2_match)) ix2_match = ix2

if (logic_option(.false., add_to_list) .and. ix_match < n_min) then
  if (n_max + 1 > ubound(nametable%name, 1)) then
    n = n_max - n_min
    call re_allocate2(nametable%name, n_min, n_max+n)
    call re_allocate2(nametable%index, n_min, n_max+n)
  endif
  if (present(has_been_added)) has_been_added = .true.
  nametable%index(ix2+1:n_max+1) = nametable%index(ix2:n_max)
  nametable%index(ix2) = n_max + 1
  nametable%name(n_max+1) = name
  ix_match = n_max + 1
  n_max = n_max+1
endif

if (present(n_match)) then
  n_match = 0
  if (ix_match == n_min - 1) return
  do
    if (nametable%name(nametable%index(ix2)) /= name) return
    n_match = n_match +1
    ix2 = ix2 + 1
    if (ix2 > n_max) return
  enddo
endif

end subroutine find_index1

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine find_index2 (name, names, an_index, n_max, ix_match, ix2_match, add_to_list, n_match, has_been_added)
!
! Subroutine to find a matching name in a list of names.
! This name is overloaded by find_index.
! See the find_index documentation for more details.
!-

subroutine find_index2 (name, names, an_index, n_max, ix_match, ix2_match, add_to_list, n_match, has_been_added)

integer ix1, ix2, ix3, n_max, ix_match
integer, optional :: ix2_match
integer an_index(:)
integer, optional :: n_match

character(*) name, names(:)

logical, optional :: add_to_list, has_been_added

! 

call find_index3 (name, names, 1, an_index, n_max, ix_match, ix2_match, add_to_list, n_match, has_been_added)

end subroutine find_index2

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine find_index3 (name, names, n_min, an_index, n_max, ix_match, ix2_match, add_to_list, n_match, has_been_added)
!
! Subroutine to find a matching name in a list of names.
! This name is overloaded by find_index.
! See the find_index documentation for more details.
!-

subroutine find_index3 (name, names, n_min, an_index, n_max, ix_match, ix2_match, add_to_list, n_match, has_been_added)

integer ix1, ix2, ix3, n_min, n_max, ix_match
integer, optional :: ix2_match
integer an_index(n_min:)

character(*) name, names(n_min:)

logical, optional :: add_to_list, has_been_added
integer, optional :: n_match

! simple case where there are no elements

if (present(has_been_added)) has_been_added = .false.

if (n_max < n_min) then
  if (present(ix2_match)) ix2_match = n_min
  ix_match = n_min - 1
  if (logic_option(.false., add_to_list)) then
    ix_match = n_min
    an_index(n_min) = n_min
    names(n_min) = name
    n_max = n_min
    if (present(n_match)) n_match = 1
    if (present(has_been_added)) has_been_added = .true.
  else
    if (present(n_match)) n_match = 0
  endif
  return
endif

!

ix1 = n_min
ix3 = n_max

do

  ix2 = (ix1 + ix3) / 2 

  if (names(an_index(ix2)) == name) then
    do ! if there are duplicate names in the list choose the first one
      if (ix2 == n_min) exit
      if (names(an_index(ix2-1)) /= names(an_index(ix2))) exit
      ix2 = ix2 - 1
    enddo
    ix_match = an_index(ix2)
    exit
  elseif (names(an_index(ix2)) < name) then
    ix1 = ix2 + 1
  else
    ix3 = ix2 - 1
  endif
                     
  if (ix1 > ix3) then
    ix_match = n_min - 1
    if (names(an_index(ix2)) < name) ix2 = ix2 + 1
    exit
  endif

enddo

if (present(ix2_match)) ix2_match = ix2

if (logic_option(.false., add_to_list) .and. ix_match < n_min) then
  if (present(has_been_added)) has_been_added = .true.
  an_index(ix2+1:n_max+1) = an_index(ix2:n_max)
  an_index(ix2) = n_max + 1
  names(n_max+1) = name
  ix_match = n_max + 1
  n_max = n_max+1
endif

if (present(n_match)) then
  n_match = 0
  if (ix_match == n_min - 1) return
  do
    if (names(an_index(ix2)) /= name) return
    n_match = n_match +1
    ix2 = ix2 + 1
    if (ix2 > n_max) return
  enddo
endif

end subroutine find_index3

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine indexer_str (str_arr, indx_arr)
!
! Subroutine to sort a character array.
! This subroutine is used to overload the generic name indexer.
! This is an extension of the Numerical Recipes routine indexer.
!
! See also find_index
!
! Input:
!   str_arr(:)  -- character(*): Array of strings
!
! Output:
!   indx_arr(:) -- integer: sorted list such that For all j:
!                    str_arr(i1) < str_arr(i2) where i1 = indx_arr(j), i2 = indx_arr(j+1).
!-

subroutine indexer_str(str_arr, indx_arr)

character(*) :: str_arr(:)
integer :: indx_arr(:)
integer, parameter :: nn = 15, nstack = 50
character(len(str_arr)) :: a
integer :: n, k, i, j, indext, jstack, l, r
integer, dimension(nstack) :: istack
n = assert_equal([size(indx_arr), size(str_arr)], 'indexer_str')
indx_arr = [(k, k = 1, n)]
jstack = 0
l = 1
r = n
do
  if (r-l < NN) then
    do j = l+1, r
      indext = indx_arr(j)
      a = str_arr(indext)
      do i = j-1, l, -1
        if (str_arr(indx_arr(i)) <= a) exit
        indx_arr(i+1) = indx_arr(i)
      end do
      indx_arr(i+1) = indext
    end do
    if (jstack == 0) RETURN
    r = istack(jstack)
    l = istack(jstack-1)
    jstack = jstack-2
  else
    k = (l+r)/2
    call swap(indx_arr(k), indx_arr(l+1))
    call icomp_xchg(indx_arr(l), indx_arr(r))
    call icomp_xchg(indx_arr(l+1), indx_arr(r))
    call icomp_xchg(indx_arr(l), indx_arr(l+1))
    i = l+1
    j = r
    indext = indx_arr(l+1)
    a = str_arr(indext)
    do
      do
        i = i+1
        if (str_arr(indx_arr(i)) >= a) exit
      end do
      do
        j = j-1
        if (str_arr(indx_arr(j)) <= a) exit
      end do
      if (j < i) exit
      call swap(indx_arr(i), indx_arr(j))
    end do
    indx_arr(l+1) = indx_arr(j)
    indx_arr(j) = indext
    jstack = jstack+2
    if (jstack > NSTACK) call err_exit('indexer: NSTACK too small')
    if (r-i+1 >= j-l) then
      istack(jstack) = r
      istack(jstack-1) = i
      r = j-1
    else
      istack(jstack) = j-1
      istack(jstack-1) = l
      l = i
    end if
  end if
end do

!-----------------
contains
subroutine icomp_xchg(i, j)
integer, intent(inout) :: i, j
integer :: swp
if (str_arr(j) < str_arr(i)) then
  swp = i
  i = j
  j = swp
end if
end subroutine icomp_xchg
end subroutine indexer_str

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine indexer_var_str (str_arr, indx_arr)
!
! Subroutine to sort a character array.
! This subroutine is used to overload the generic name indexer.
! This is an extension of the Numerical Recipes routine indexer.
!
! See also find_index
!
! Input:
!   str_arr(:)  -- var_length_string_struct: Array of var_strings
!
! Output:
!   indx_arr(:) -- integer: sorted list such that:
!                 For all j, str_arr(i1) < str_arr(i2) where i1 = indx_arr(j), i2 = indx_arr(j+1).
!-

subroutine indexer_var_str(str_arr, indx_arr)

type(var_length_string_struct) :: str_arr(:)
integer :: indx_arr(:)
integer, parameter :: nn = 15, nstack = 50
type(var_length_string_struct) :: a
integer :: n, k, i, j, indext, jstack, l, r
integer, dimension(nstack) :: istack
n = assert_equal([size(indx_arr), size(str_arr)], 'indexer_var_str')
indx_arr = [(k, k = 1, n)]
jstack = 0
l = 1
r = n
do
  if (r-l < NN) then
    do j = l+1, r
      indext = indx_arr(j)
      a = str_arr(indext)
      do i = j-1, l, -1
        if (str_arr(indx_arr(i))%str <= a%str) exit
        indx_arr(i+1) = indx_arr(i)
      end do
      indx_arr(i+1) = indext
    end do
    if (jstack == 0) RETURN
    r = istack(jstack)
    l = istack(jstack-1)
    jstack = jstack-2
  else
    k = (l+r)/2
    call swap(indx_arr(k), indx_arr(l+1))
    call icomp_xchg(indx_arr(l), indx_arr(r))
    call icomp_xchg(indx_arr(l+1), indx_arr(r))
    call icomp_xchg(indx_arr(l), indx_arr(l+1))
    i = l+1
    j = r
    indext = indx_arr(l+1)
    a = str_arr(indext)
    do
      do
        i = i+1
        if (str_arr(indx_arr(i))%str >= a%str) exit
      end do
      do
        j = j-1
        if (str_arr(indx_arr(j))%str <= a%str) exit
      end do
      if (j < i) exit
      call swap(indx_arr(i), indx_arr(j))
    end do
    indx_arr(l+1) = indx_arr(j)
    indx_arr(j) = indext
    jstack = jstack+2
    if (jstack > NSTACK) call err_exit('indx_arrx: NSTACK too small')
    if (r-i+1 >= j-l) then
      istack(jstack) = r
      istack(jstack-1) = i
      r = j-1
    else
      istack(jstack) = j-1
      istack(jstack-1) = l
      l = i
    end if
  end if
end do

!-----------------
contains

subroutine icomp_xchg(i, j)
integer, intent(inout) :: i, j
integer :: swp
if (str_arr(j)%str < str_arr(i)%str) then
  swp = i
  i = j
  j = swp
end if
end subroutine icomp_xchg
end subroutine indexer_var_str

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine indexer_real(arr, index)

implicit none
real(rp), dimension(:), intent(in) :: arr
integer, dimension(:), intent(out) :: index
integer, parameter :: nn = 15, nstack = 50
real(rp) :: a
integer :: n, k, i, j, indext, jstack, l, r
integer, dimension(nstack) :: istack
n = assert_equal([size(index), size(arr)], 'indexer_real')
index = [(k, k = 1, n)]
jstack = 0
l = 1
r = n
do
  if (r-l < NN) then
    do j = l+1, r
      indext = index(j)
      a = arr(indext)
      do i = j-1, l, -1
        if (arr(index(i)) <= a) exit
        index(i+1) = index(i)
      end do
      index(i+1) = indext
    end do
    if (jstack == 0) RETURN
    r = istack(jstack)
    l = istack(jstack-1)
    jstack = jstack-2
  else
    k = (l+r)/2
    call swap(index(k), index(l+1))
    call icomp_xchg(index(l), index(r))
    call icomp_xchg(index(l+1), index(r))
    call icomp_xchg(index(l), index(l+1))
    i = l+1
    j = r
    indext = index(l+1)
    a = arr(indext)
    do
      do
        i = i+1
        if (arr(index(i)) >= a) exit
      end do
      do
        j = j-1
        if (arr(index(j)) <= a) exit
      end do
      if (j < i) exit
      call swap(index(i), index(j))
    end do
    index(l+1) = index(j)
    index(j) = indext
    jstack = jstack+2
    if (jstack > NSTACK) call err_exit('indexer: NSTACK too small')
    if (r-i+1 >= j-l) then
      istack(jstack) = r
      istack(jstack-1) = i
      r = j-1
    else
      istack(jstack) = j-1
      istack(jstack-1) = l
      l = i
    end if
  end if
end do

!-------------
contains

subroutine icomp_xchg(i, j)
integer, intent(inout) :: i, j
integer :: swp
if (arr(j) < arr(i)) then
  swp = i
  i = j
  j = swp
end if
end subroutine icomp_xchg

end subroutine indexer_real

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine indexer_int(iarr, index)

implicit none
integer, intent(in) :: iarr(:)
integer, intent(out) :: index(:)
integer, parameter :: nn = 15, nstack = 50
integer :: a
integer :: n, k, i, j, indext, jstack, l, r
integer, dimension(nstack) :: istack
n = assert_equal([size(index), size(iarr)], 'indexer_int')
index = [(k, k = 1, n)]
jstack = 0
l = 1
r = n
do
  if (r-l < nn) then
    do j = l+1, r
      indext = index(j)
      a = iarr(indext)
      do i = j-1, l, -1
        if (iarr(index(i)) <= a) exit
        index(i+1) = index(i)
      end do
      index(i+1) = indext
    end do
    if (jstack == 0) return
    r = istack(jstack)
    l = istack(jstack-1)
    jstack = jstack-2
  else
    k = (l+r)/2
    call swap(index(k), index(l+1))
    call icomp_xchg(index(l), index(r))
    call icomp_xchg(index(l+1), index(r))
    call icomp_xchg(index(l), index(l+1))
    i = l+1
    j = r
    indext = index(l+1)
    a = iarr(indext)
    do
      do
        i = i+1
        if (iarr(index(i)) >= a) exit
      end do
      do
        j = j-1
        if (iarr(index(j)) <=  a) exit
      end do
      if (j < i) exit
      call swap(index(i), index(j))
    end do
    index(l+1) = index(j)
    index(j) = indext
    jstack = jstack+2
    if (jstack > NSTACK) call err_exit('indexer: NSTACK too small')
    if (r-i+1 >= j-l) then
      istack(jstack) = r
      istack(jstack-1) = i
      r = j-1
    else
      istack(jstack) = j-1
      istack(jstack-1) = l
      l = i
    end if
  end if
end do

!------------------------------
contains

subroutine icomp_xchg(i, j)
integer, intent(inout) :: i, j
integer :: swp
if (iarr(j) < iarr(i)) then
  swp = i
  i = j
  j = swp
end if
end subroutine icomp_xchg

end subroutine indexer_int

end module
