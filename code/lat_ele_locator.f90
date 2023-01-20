!+
! Subroutine lat_ele_locator (loc_str, lat, eles, n_loc, err, above_ubound_is_err, ix_dflt_branch, order_by_index)
!
! Routine to locate all the elements in a lattice that corresponds to loc_str.
! If there are multiple elements of the same name, pointers to all such elements are returned.
!
! loc_str is a list of element names or element ranges.
! A space or a comma delimits the elements.
! An ampersand "&" can be used to form the intersection of two groups.
!
! An element "name" (which can match to multiple elements) can be of the form
!   {~}{branch>>}{key::}ele_id{##N}{+/-offset}
! Where
!   ~         -- Negation character. See below.
!   key       -- Optional key name ("quadrupole", "sbend", etc.)
!   branch    -- Name or index of branch. May contain the wild cards "*" and "%".
!   ele_id    -- Name or index of element. May contain the wild cards "*" and "%".
!                 If a name and no branch is given, all branches are searched.
!                 If an index and no branch is given, branch 0 is assumed.
!   ##N       -- N = integer. N^th instance of ele_id in the branch.
!   +/-offset -- Element offset. For example, "Q1+1" is the element after "Q1" and 
!                 "Q1-2" is the second element before "Q1".
! Note: An old syntax that is still supported is:
!   {key::}{branch>>}ele_id{##N}
!
! An element range is of the form:
!   {key::}ele1:ele2
! Where:
!   key      -- Optional key name ("quadrupole", "sbend", etc.). 
!               Also key may be "type", "alias", or "descrip" in which case the %type, %alias, or 
!               %descrip field is matched to instead of the element name.
!   ele1     -- Starting element of the range.
!   ele2     -- Ending element of the range. 
! The range includes ele1 and ele2.
! If ele2%ix_ele < ele1%ix_ele, the range wraps around the ends of the lattice. That is, the range is
!   [ele1%ix_ele, branch%n_ele_track] + [0, ele2%ix_ele]
! Note: ele1 and ele2 must be in the same branch. Branch 0 is the default.
! If ele1 or ele2 is a super_lord, the elements in the range are determined by the position of the super_slave elements.
! For example, if loc_str is "Q1:Q1" and Q1 is *not* a super_lord, the eles list will simply be Q1.
! If Q1 is a super_lord, the eles list will be the super_slaves of Q1.
! It is an error if ele1 or ele2 is a multipass_lord. Also it is is an error for ele1 or ele2
! to be an overlay, group, or girder if the number of slaves is not one. 
!
! Note: Use quotation marks for matching to type, alias, or descrip strings that have a blank.
! Note: matching is always case insensitive.
!
! Examples:
!   x_br>>quad::q*     All quadrupoles of branch "x_br" whose name begins with "q"
!   3,15:17            Elements with index 3, 15, 16, and 17 in branch 0.
!   2>>45:51           Elements 45 through 51 of branch 2.
!   q1:q5              Elements between "q1" and "q5".
!   sbend::q1:q5       All sbend elements between "q1" and "q5".
!   marker::a*##2      2^nd marker element in each branch whose name begins with "a".
!   type::bpm*         All elements whose %type field starts with bpm.
!   alias::'my duck'   Match to all elements whose %alias matches "my duck"
!
! The negation character "~" signifies that elements are to be removed from the list. For example:
!   "quadrupole::*, ~q3"
! The element list would be all quadrupoles except quadrupoles named "Q3". Example:
!   "*::*, ~octupole::*, oct23"
! Here the list would be all elements except the octupole elements that are not named "oct23".
! This shows that order is important when negation is used since adding/subtracting elements from
! the list is done left to right.
!
! An ampersand "&" can be used to form the intersection of two groups.
! Example:
!   1:10 & BPM*   ! All element with index in range [1,10] and with name starting with "BPM".
! 
! Note: For something like loc_str = "quad::*", if order_by_index = True, the eles(:) array will
! be ordered by element index. If order_by_index = False, the eles(:) array will be ordered by
! s-position. This is the same as order by index except in the case where where there are super_lord
! elements. Since super_lord elements always have a greater index (at least in branch 0), order by index
! will place any super_lord elements at the end of the list.
!
! Note: When there are multiple element names in loc_str (which will be separated by a comma or blank), 
! the elements in the eles(:) array will be in the same order as they appear loc_str. For example,
! with loc_str = "quad::*,sbend::*", all the quadrupoles will appear in eles(:) before all of the sbends.
! This is independent of the setting of order_by_index but the setting of order_by_index can affect
! the ordering of the quads among themselves and the ordering of the sbends among themselves.
! 
! Input:
!   loc_str        -- character(*): Element names or indexes. May be lower case.
!   lat            -- lat_struct: Lattice to search through.
!   above_ubound_is_err
!                  -- logical, optional: If the upper bound "e2" on an "e1:e2" range construct 
!                       is above the maximum element index then treat this as an error? 
!                       Default is True. If False, treat e2 as the maximum element index. 
!   ix_dflt_branch -- integer, optional: If present and not -1 then restrict search to specified branch.
!                       If not present or -1: Search all branches. Exception: For elements specified using 
!                       an integer index (EG: "43"), if ix_dflt_branch is not present or -1 use branch 0.
!   order_by_index -- logical, optional: False is default. If True, order a component of loc_str like "quad::*" 
!                       by element index instead of longitudinal s-position. Index ordering and s-position ordering
!                       are different when there are super lords and super slaves.
!
! Output:
!   eles(:)       -- ele_pointer_struct, allocatable: Array of matching elements.
!                    Note: This routine does not try to deallocate eles.
!                     It is up to you to deallocate eles if needed.
!   n_loc         -- integer: Number of locations found.
!                      Set to zero if no elements are found.
!   err           -- logical, optional: Set True if there is a decode error.
!                      Note: Not finding any matching element is not an error.
!-

subroutine lat_ele_locator (loc_str, lat, eles, n_loc, err, above_ubound_is_err, ix_dflt_branch, order_by_index)

use bmad_interface, dummy => lat_ele_locator

implicit none

type (lat_struct), target :: lat
type (ele_pointer_struct), allocatable :: eles(:)
type (ele_pointer_struct), allocatable, target :: eles2(:)
type (ele_struct), pointer :: ele_start, ele_end
type (branch_struct), pointer :: branch

character(*) loc_str
character(200) str
character(60) name, name2
character(8) branch_str
character(*), parameter :: r_name = 'lat_ele_locator'
character(1) delim

integer, optional :: ix_dflt_branch
integer i, j, ib, ios, ix, n_loc, n_loc2, match_name_to
integer in_range, ix_word, key
integer, parameter :: ele_name$ = -1

logical, optional :: above_ubound_is_err, err, order_by_index
logical err2, delim_found, above_ub_is_err, s_ordered, negate, names_are_integers, intersection

! init

if (present(err)) err = .true.
n_loc = 0
intersection = .false.
branch_str = ''
delim = ''
str = loc_str
call str_upcase (str, str)
! In_range:
!   0 -> not in range construct
!   1 -> processing start
!   2 -> processing stop
in_range = 0   

! Loop over all items in the list

do

  ! Get next item. 
  ! If the split is in between "::" then need to piece the two haves together.

  if (delim == '&') intersection = .true.

  call word_read (str, ':, &', name, ix_word, delim, delim_found, str, ignore_interior = .true.)

  ! Need to handle old style "class::branch>>name" and new style "branch>>class::name" and not be
  ! confused by "name1:name2" range syntax

  if (in_range == 0) then
    negate = .false.
    key = 0
    match_name_to = ele_name$
  endif

  if (name(1:1) == '~') then
    negate = .true.
    name = name(2:)
  endif

  ix = index(name, '>>')
  if (ix /= 0) then
    if (in_range /= 0) then
      call out_io (s_error$, r_name, 'BAD BRANCH/RANGE CONSTRUCT')
      return
    endif

    branch_str = name(:ix-1)
    name = name(ix+2:)
  endif

  ! Have key?
  if (str(1:1) == ':') then
    if (in_range /= 0) then
      call out_io (s_error$, r_name, 'BAD CLASS/RANGE CONSTRUCT')
      return
    endif

    key = extended_key_name_to_key_index (name, .true., match_name_to)
    if (key < 0) then
      call out_io (s_error$, r_name, 'BAD ELEMENT CLASS: ' // name)
      return
    endif

    call word_read (str(2:), ':, &', name, ix_word, delim, delim_found, str, ignore_interior = .true.)
    if (name == '') then
      call out_io (s_error$, r_name, '"CLASS::NAME" CONSTRUCT NOT VALID WITH BLANK NAME. IN: ' // loc_str, &
                                     'Note: USE "*" TO MATCH TO ALL NAMES')
      return
    endif
  endif

  !

  if (name(1:1) == "'" .or. name(1:1) == '"') then
    name = unquote(name)

  else
    ix = index(name, '>>')
    if (ix /= 0) then
      branch_str = name(:ix-1)
      name = name(ix+2:)
    endif
  endif

  if (delim == ':' .or. in_range /= 0) in_range = in_range + 1 
  if (name == '') exit

  ! Get list of elements for this item

  ! In a range the key is applied to the list of elements in the range and not
  ! applied to picking the end elements.
  above_ub_is_err = (logic_option(.true., above_ubound_is_err) .or. in_range /= 2)
  s_ordered = (.not. logic_option(.false., order_by_index))
  if (in_range == 0) then
    call lat_ele1_locator (branch_str, key, name, match_name_to, lat, eles2, n_loc2, err2, &
                                                     above_ub_is_err, ix_dflt_branch, s_ordered)
  else
    call lat_ele1_locator (branch_str, 0, name, match_name_to, lat, eles2, n_loc2, err2, &
                                                     above_ub_is_err, ix_dflt_branch, s_ordered)
    if (in_range == 1) then
      names_are_integers = (is_integer(name))
    else
      names_are_integers = (is_integer(name) .and. names_are_integers)
    endif
  endif
  if (err2) return

  ! If not a range construct then just put this on the list

  if (in_range == 0) then
    if (negate) then
      call negate_eles(eles, n_loc, eles2(1:n_loc2))
    elseif (intersection) then
      call intersection_eles (eles, n_loc, eles2(1:n_loc2))
      intersection = .false.
    else
      if (.not. allocated(eles)) allocate(eles(n_loc+n_loc2))
      call re_allocate_eles(eles, n_loc+n_loc2, .true.)
      eles(n_loc+1:n_loc+n_loc2) = eles2(1:n_loc2)
      n_loc = n_loc + n_loc2
    endif
    cycle
  endif

  ! Handle range construct...
  ! Each item in the construct must represent only one element.

  if (n_loc2 == 0) then
    if (key_name_to_key_index(name, .false.) > 0) then
      call out_io (s_error$, r_name, 'NO ELEMENT ASSOCIATED WITH: ' // name, &
                                     'THIS IS NOT PERMITTED IN A RANGE CONSTRUCT: ' // loc_str, &
                                     'IF THIS IS AN ELEMENT CLASS YOU NEED TWO COLONS "::" AND NOT ONE.')
    else
      call out_io (s_error$, r_name, 'NO ELEMENT ASSOCIATED WITH: ' // name, &
                                     'THIS IS NOT PERMITTED IN A RANGE CONSTRUCT: ' // loc_str)
    endif
    return
  elseif (n_loc2 > 1) then
    call out_io (s_error$, r_name, 'MULTIPLE ELEMENTS ASSOCIATED WITH: ' // name, &
                                   'THIS IS NOT PERMITTED IN A RANGE CONSTRUCT: ' // loc_str)
    return
  endif

  if (in_range == 1) ele_start => eles2(1)%ele
  if (in_range == 2) ele_end   => eles2(1)%ele

  if (delim == ':') then
    if (in_range >= 3) then
      call out_io (s_error$, r_name, 'TOO MANY ":" IN RANGE CONSTRUCT: ' // loc_str)
      return
    endif
    cycle
  endif

  ! if we have the range then add it to the eles list.
  ! Normally a range "A:B" (for example, "Q1:Q2") is interpreted so that if "A" or "B" is a lord, 
  ! use the corresponding slave To construct a range.
  ! Exception: Something like "1:100" where integers are used is interpreted as elements 1 through 100.

  ele_end => eles2(1)%ele
  if (.not. names_are_integers) then
    ele_start => find_this_end(ele_start, entrance_end$, err2); if (err2) return
    ele_end => find_this_end(ele_end, exit_end$, err2);  if (err2) return
  endif

  if (ele_start%ix_branch /= ele_end%ix_branch) then
    call out_io (s_error$, r_name, 'ELEMENTS NOT OF THE SAME BRANCH IN RANGE: ' // loc_str)
    return
  endif
  branch => lat%branch(ele_end%ix_branch)
  ib = branch%ix_branch

  call re_allocate_eles(eles2, branch%n_ele_track, .false.)
  n_loc2 = 0
  if (ele_start%ix_ele > ele_end%ix_ele) then
    do i = ele_start%ix_ele, branch%n_ele_track
      if (key > 0 .and. branch%ele(i)%key /= key .and. branch%ele(i)%sub_key /= key) cycle
      n_loc2 = n_loc2 + 1; eles2(n_loc2)%ele => branch%ele(i)
    enddo

    do i = 0, ele_end%ix_ele
      if (key > 0 .and. branch%ele(i)%key /= key .and. branch%ele(i)%sub_key /= key) cycle
      n_loc2 = n_loc2 + 1; eles2(n_loc2)%ele => branch%ele(i)
    enddo

  else
    do i = ele_start%ix_ele, ele_end%ix_ele
      if (key > 0 .and. branch%ele(i)%key /= key .and. branch%ele(i)%sub_key /= key) cycle
      n_loc2 = n_loc2 + 1; eles2(n_loc2)%ele => branch%ele(i)
    enddo
  endif

  if (negate) then
    call negate_eles(eles, n_loc, eles2(1:n_loc2))

  elseif (intersection) then
    call intersection_eles(eles, n_loc, eles2(1:n_loc2))
    intersection = .false.

  else
    call re_allocate_eles(eles, n_loc+n_loc2, .true.)
    do i = 1, n_loc2
      eles(n_loc+i)%ele => eles2(i)%ele
      eles(n_loc+i)%loc = lat_ele_loc_struct(eles2(i)%ele%ix_ele, ib)
    enddo
    n_loc = n_loc + n_loc2
  endif

  in_range = 0
enddo

! Check that no ranges are open.

if (in_range == 1) then
  call out_io (s_error$, r_name, 'MISSING ELEMENT AFTER ":" IN: ' // loc_str)
  return
endif

if (present(err)) err = .false.

!---------------------------------------------------------------------------
contains

!+
! Subroutine lat_ele1_locator (branch_str, key, name_in, match_name_to, lat, eles, n_loc, err, &
!                                                                      above_ub_is_err, ix_dflt_branch, s_ordered)
!
! Routine to locate all the elements in a lattice that corresponds to loc_str. 
! Note: This routine is private and meant to be used only by lat_ele_locator
!-

subroutine lat_ele1_locator (branch_str, key, name_in, match_name_to, lat, eles, n_loc, err, &
                                                                      above_ub_is_err, ix_dflt_branch, s_ordered)

use super_recipes_mod

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch, branch2
type (ele_pointer_struct), allocatable :: eles(:)
type (ele_struct), pointer :: ele, slave
type (nametable_struct), pointer :: nt

character(*) branch_str, name_in
character(40) name
character(*), parameter :: r_name = 'lat_ele1_locator'

integer, optional :: ix_dflt_branch
integer i, n, ib, ix, ix_branch, ix_branch_old, ixp, ios, ix_ele, n_loc, n_max, ix1, ix2
integer key, target_instance, n_instance_found, match_name_to, ix_slave, offset
integer, allocatable :: ix_nt(:)

logical above_ub_is_err, s_ordered
logical err, do_match_wild, Nth_instance_found

! init

err = .true.
n_loc = 0
name = name_in
offset = 0

! Look for "+N" or "-N" offset suffix.

ix = index(name, '+')
if (ix /= 0) then
  if (.not. is_integer(name(ix+1:), offset)) then
    call out_io (s_error$, r_name, 'INVALID OFFSET SYNTAX: ' // name_in)
    return
  endif
  name = name(:ix-1)
endif

ix = index(name, '-')
if (ix /= 0) then
  if (.not. is_integer(name(ix+1:), offset)) then
    call out_io (s_error$, r_name, 'INVALID OFFSET SYNTAX: ' // name_in)
    return
  endif
  offset = -offset
  name = name(:ix-1)
endif

! Look for "##N" suffix to indicate which instance to choose from when there
! are multiple elements of a given name.

ix = index(name, '##')
target_instance = 0
if (ix /= 0 .and. match_name_to == ele_name$ .and. is_integer(name(ix+2:))) then
  read (name(ix+2:), *) target_instance
  name = name(:ix-1)
endif

! Read branch name or index which is everything before an '>>'.

ix_branch = integer_option (-1, ix_dflt_branch)
if (ix_branch < -1 .or. ix_branch > ubound(lat%branch, 1)) then
  call out_io (s_error$, r_name, 'BRANCH INDEX OUT OF RANGE: \i0\ ', i_array = [ix_branch])
  return
endif

! Note: Branch name not found is treated the same as element name not found.
! That is, no match and is not an error

if (branch_str /= '') then
  if (index(branch_str, '*') == 0 .and. index(branch_str, '%') == 0) then
    branch => pointer_to_branch (branch_str, lat)
    if (.not. associated(branch)) then
      err = .false.
      return
    endif
    ix_branch = branch%ix_branch
  else
    ix_branch = -2    ! Match to branch_str
  endif
endif

! Read integer if present

if (is_integer(name) .and. match_name_to == ele_name$) then
  read (name, *, iostat = ios) ix_ele
  if (ios /= 0) then
    call out_io (s_error$, r_name, 'BAD ELEMENT LOCATION INDEX: ' // name)
    return
  endif

  if (ix_branch == -2) then  ! Match to branch_str
    do ib = 0, ubound(lat%branch, 1)
      branch => lat%branch(ib)
      if (.not. match_wild(branch%name, branch_str)) cycle
      if (.not. above_ub_is_err .and. ix_ele > branch%n_ele_max) then
        ix_ele = branch%n_ele_max
      elseif (ix_ele < 0 .or. ix_ele > branch%n_ele_max) then
        call out_io (s_error$, r_name, 'ELEMENT LOCATION INDEX OUT OF RANGE: ' // name)
        return
      endif
      n_loc = n_loc + 1
      if (.not. allocated(eles) .or. size(eles) < n_loc) call re_allocate_eles (eles, 2*n_loc, .true.)
      eles(n_loc)%ele => branch%ele(ix_ele)
    enddo

  else
    if (ix_branch == -1) ix_branch = 0
    branch => lat%branch(ix_branch)
    if (.not. above_ub_is_err .and. ix_ele > branch%n_ele_max) then
      ix_ele = branch%n_ele_max
    elseif (ix_ele < 0 .or. ix_ele > branch%n_ele_max) then
      call out_io (s_error$, r_name, 'ELEMENT LOCATION INDEX OUT OF RANGE: ' // name)
      return
    endif
    n_loc = 1
    call re_allocate_eles (eles, 1)
    eles(1)%ele => branch%ele(ix_ele)
  endif

  call add_offset (eles, n_loc, offset, err)
  return
endif

! search for matches.

do_match_wild = (index(name, "*") /= 0 .or. index(name, "%") /= 0)
Nth_instance_found = .false.

! If a element name search then can use nametable to quickly find the matches.

if (.not. do_match_wild .and. match_name_to == ele_name$) then
  nt => lat%nametable
  call find_index(name, nt, ix, ix2_match = ix1)
  if (ix < 0) then
    err = .false.
    return
  endif
  ix2 = ix1

  do
    if (ix2 == nt%n_max) exit
    if (nt%name(nt%index(ix2+1)) /= name) exit
    ix2 = ix2 + 1
  enddo

  allocate (ix_nt(ix1:ix2))
  do i = ix1, ix2
    ix_nt(i) = nt%index(i)
  enddo

  ! Due to the way check_this_match works, put first slaves in the list
  if (s_ordered) then
    i = ix1 - 1
    do
      i = i + 1
      if (i > ix2) exit
      ele => pointer_to_ele(lat, ix_nt(i))
      if (ele%n_slave == 0) cycle
      slave => pointer_to_slave (ele, 1)
      do while (slave%n_slave /= 0)
        slave => pointer_to_slave(slave, 1)
      enddo
      j = ele_nametable_index(slave)
      if (any(j == ix_nt)) cycle  ! No duplicates
      ix2 = ix2 + 1
      call re_allocate2(ix_nt, ix1, ix2)
      ix_nt(ix2) = j
    enddo
  endif

  call super_sort(ix_nt)

  ix_branch_old = -1
  do i = ix1, ix2
    ele => pointer_to_ele(lat, ix_nt(i))
    branch => lat%branch(ele%ix_branch)
    if (ix_branch == -2 .and. .not. match_wild(branch%name, branch_str)) cycle
    if (ix_branch > -1 .and. ele%ix_branch /= ix_branch) cycle
    if (ele%ix_branch /= ix_branch_old) then
      n_instance_found = 0
      ix_branch_old = ele%ix_branch
      Nth_instance_found = .false.
    endif
    if (s_ordered .and. ele%ix_branch == 0 .and. ele%ix_ele > branch%n_ele_track) cycle
    call check_this_match (name, ele, key, do_match_wild, match_name_to, n_instance_found, n_loc, eles, &
                                                            target_instance, Nth_instance_found, s_ordered)
  enddo

  ! Elements passed over in above loop

  if (s_ordered .and. .not. Nth_instance_found .and. (ix_branch == -1 .or. ix_branch == 0)) then
    do i = ix1, ix2
      ele => pointer_to_ele(lat, ix_nt(i))
      if (ele%n_slave /= 0 .or. ele%ix_branch /= 0 .or. ele%ix_ele <= lat%n_ele_track) cycle
      call check_this_match(name, ele, key, do_match_wild, match_name_to, n_instance_found, n_loc, eles, &
                                                              target_instance, Nth_instance_found, s_ordered)
    enddo
  endif

  call add_offset (eles, n_loc, offset, err)
  return
endif

! Not using the nametable.

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  if (ix_branch == -2 .and. .not. match_wild(branch%name, branch_str)) cycle
  if (ix_branch > -1 .and. ib /= ix_branch) cycle

  n_instance_found = 0
  Nth_instance_found = .false.

  n_max = branch%n_ele_max
  if (s_ordered) n_max = branch%n_ele_track

  do i = 0, n_max
    ele => branch%ele(i)
    call check_this_match (name, ele, key, do_match_wild, match_name_to, n_instance_found, n_loc, eles, &
                                                            target_instance, Nth_instance_found, s_ordered)
    if (Nth_instance_found) exit
  enddo
enddo 

! During lattice parsing, elements may get added to the lord section for temporary storage and
! these elements will not have any slaves.
! If S-ordered, these elements have been so far overlooked and this has to be rectified.

if (s_ordered .and. .not. Nth_instance_found .and. (ix_branch == -1 .or. ix_branch == 0)) then
  do i = lat%n_ele_track+1, lat%n_ele_max
    ele => lat%ele(i)
    if (ele%n_slave /= 0) cycle
    call check_this_match(name, ele, key, do_match_wild, match_name_to, n_instance_found, n_loc, eles, &
                                                            target_instance, Nth_instance_found, s_ordered)
  enddo
endif

call add_offset (eles, n_loc, offset, err)

end subroutine lat_ele1_locator

!---------------------------------------------------------------------------
! contains

subroutine add_offset (eles, n_loc, offset, err)

type (ele_pointer_struct), allocatable, target :: eles(:)
type (ele_struct), pointer :: ele
type (branch_struct), pointer :: branch
integer n_loc, offset, ie, nl
logical err

!

if (offset == 0) then
  err = .false.
  return
endif

!

do ie = 1, n_loc
  ele => eles(ie)%ele

  select case (ele%lord_status)
  case (multipass_lord$, overlay_lord$, group_lord$, ramper_lord$)
    call out_io (s_error$, r_name, 'OFFSET CANNOT BE APPLIED TO MULTIPASS, OVERLAY, OR GROUP LORD: ' // ele%name)
    return

  case (super_lord$, girder_lord$)
    if (offset > 0) then
      ele => pointer_to_slave(ele, ele%n_slave)
    else
      ele => pointer_to_slave(ele, 1)
    endif
  end select

  branch => ele%branch
  nl = ele%ix_ele + offset
  if (nl < 0) nl = ele%branch%n_ele_track + nl + 1
  if (nl > branch%n_ele_track) nl = nl - branch%n_ele_track - 1
  eles(ie)%ele => branch%ele(nl)
enddo

err = .false.

end subroutine add_offset

!---------------------------------------------------------------------------
! contains

recursive subroutine check_this_match (name, ele, key, do_match_wild, match_name_to, n_instance_found, n_loc, eles, &
                                                            target_instance, Nth_instance_found, s_ordered)

type (ele_struct), target :: ele
type (ele_struct), pointer :: lord
type (ele_pointer_struct), allocatable :: eles(:)

integer key, match_name_to, target_instance, ix_slave, n_instance_found, n_loc, il
character(*) name
logical do_match_wild, Nth_instance_found, s_ordered

! S-ordered means check lords before slaves

if (s_ordered) then
  do il = 1, ele%n_lord
    lord => pointer_to_lord(ele, il, ix_slave_back = ix_slave)
    if (ix_slave /= 1) cycle
    call check_this_match(name, lord, key, do_match_wild, match_name_to, n_instance_found, n_loc, eles, &
                                                             target_instance, Nth_instance_found, s_ordered)
    if (Nth_instance_found) return
  enddo
endif

!

if ((key /= 0 .and. ele%key /= key) .and. (ele%key /= sbend$ .or. ele%sub_key /= key)) return

if (do_match_wild) then
  select case (match_name_to)
  case (alias$)
    if (ele%alias == '') return
    if (.not. match_wild(upcase(ele%alias), name)) return
  case (descrip$)
    if (.not. associated(ele%descrip)) return
    if (.not. match_wild(upcase(ele%descrip), name)) return
  case (type$)
    if (ele%type == '') return
    if (.not. match_wild(upcase(ele%type), name)) return
  case (ele_name$)
    if (.not. match_wild(ele%name, name)) return
  end select

else
  select case (match_name_to)
  case (alias$)
    if (upcase(ele%alias) /= name) return
  case (descrip$)
    if (.not. associated(ele%descrip)) return
    if (upcase(ele%descrip) /= name) return
  case (type$)
    if (upcase(ele%type) /= name) return
  case (ele_name$)
    if (ele%name /= name) return
  end select
endif

! Is matched

if (target_instance > 0) then
  n_instance_found = n_instance_found + 1
  if (n_instance_found /= target_instance) return
  Nth_instance_found = .true.
endif

n_loc = n_loc + 1
if (.not. allocated(eles) .or. size(eles) < n_loc) call re_allocate_eles (eles, 2*n_loc, .true.)
eles(n_loc)%ele => ele
eles(n_loc)%loc = lat_ele_loc_struct(ele%ix_ele, ele%ix_branch)

end subroutine check_this_match

!---------------------------------------------------------------------------
! contains

function find_this_end(ele, which_end, err) result (end_ele)

type (ele_struct), target :: ele
type (ele_struct), pointer :: end_ele, orig_ele
integer which_end
logical err

!

err = .true.
orig_ele => ele
end_ele => ele

do
  select case (end_ele%lord_status)
  case (overlay_lord$, group_lord$, girder_lord$)
    if (end_ele%n_slave /= 1) then
      call out_io (s_error$, r_name, &
          'FOR ELEMENT: ' // orig_ele%name, 'WHICH IS PART OF A RANGE CONSTRUCT IN: ' // loc_str, &
          'THE PROBLEM IS THAT THIS ELEMENT IS AN OVERLAY, GROUP OR GIRDER LORD AND DOES ', &
          'NOT HAVE A UNIQUE SLAVE ELEMENT TO COMPUTE THE RANGE END POSITION.')
      return
    endif
    end_ele => pointer_to_slave(end_ele, 1)

  case (ramper_lord$)
    call out_io (s_error$, r_name, &
        'FOR ELEMENT: ' // orig_ele%name, 'WHICH IS PART OF A RANGE CONSTRUCT IN: ' // loc_str, &
        'THE PROBLEM IS THAT THIS ELEMENT IS A RAMPER LORD AND DOES ', &
        'NOT HAVE A UNIQUE SLAVE ELEMENT TO COMPUTE THE RANGE END POSITION.')
    return

  case (multipass_lord$)
    if (orig_ele%lord_status == multipass_lord$) then
      call out_io (s_error$, r_name, &
          'FOR ELEMENT: ' // orig_ele%name, 'WHICH IS PART OF A RANGE CONSTRUCT IN: ' // loc_str, &
          'THE PROBLEM IS THAT THIS ELEMENT IS A MULTIPASS_LORD DOES ', &
          'NOT HAVE A UNIQUE MULTIPASS_SLAVE ELEMENT TO COMPUTE THE RANGE END POSITION.')
    else
      call out_io (s_error$, r_name, &
          'FOR ELEMENT: ' // orig_ele%name, 'WHICH IS PART OF A RANGE CONSTRUCT IN: ' // loc_str, &
          'THE PROBLEM IS THAT THIS ELEMENT IS AN OVERLAY, GROUP OR GIRDER LORD ELEMENT ', &
          'WITH A MULTIPASS_LORD AS A SLAVE SO THERE IS NO UNIQUE RANGE END POSITION.')
    endif
    return

  case (super_lord$)
    if (which_end == entrance_end$) then
      end_ele => pointer_to_slave(end_ele, 1)
    else
      end_ele => pointer_to_slave(end_ele, end_ele%n_slave)
    endif
    exit

  case (not_a_lord$)
    exit

  case default
    call out_io (s_fatal$, r_name, 'INTERNAL ERR. PLEASE REPORT!')
    stop
  end select
enddo

err = .false.

end function find_this_end

!---------------------------------------------------------------------------
! contains

function extended_key_name_to_key_index (name, abbrev_allowed, match_name_to) result (key_index)

character(*) name
integer match_name_to, key_index
logical abbrev_allowed

!

select case (name)
case ('TYPE')
  match_name_to = type$
  key_index = 0
case ('DESCRIP')
  match_name_to = descrip$
  key_index = 0
case ('ALIAS')
  match_name_to = alias$
  key_index = 0
case default
  match_name_to = ele_name$
  key_index = key_name_to_key_index(name, abbrev_allowed)
end select

end function extended_key_name_to_key_index

!---------------------------------------------------------------------------
! contains

subroutine negate_eles(eles, n_ele, eles2)

type (ele_pointer_struct) eles(:), eles2(:)
integer n_ele
integer i, n

!

do i = 1, n_ele
  eles(i)%ele%ixx = 1   ! Mark to keep
enddo

do i = 1, size(eles2)
  eles2(i)%ele%ixx = 0  ! Mark to remove
enddo

n = 0
do i = 1, n_ele
  if (eles(i)%ele%ixx == 0) cycle
  n = n + 1
  eles(n)%ele => eles(i)%ele
enddo

n_ele = n

end subroutine negate_eles

!---------------------------------------------------------------------------
! contains

subroutine intersection_eles(eles, n_ele, eles2)

type (ele_pointer_struct) eles(:), eles2(:)
integer n_ele
integer i, n

!

do i = 1, n_ele
  eles(i)%ele%ixx = 0
enddo

do i = 1, size(eles2)
  eles2(i)%ele%ixx = 1  ! Mark to keep
enddo

n = 0
do i = 1, n_ele
  if (eles(i)%ele%ixx == 0) cycle
  n = n + 1
  eles(n)%ele => eles(i)%ele
enddo

n_ele = n

end subroutine intersection_eles

end subroutine lat_ele_locator
