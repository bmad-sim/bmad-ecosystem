module lat_ele_loc_mod

use bmad_utils_mod

private lat_ele1_locator

contains

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine lat_ele_locator (loc_str, lat, eles, n_loc, err, above_ubound_is_err, ix_dflt_branch)
!
! Routine to locate all the elements in a lattice that corresponds to loc_str. 
!
! If there are multiple elements of the same name, pointers to all such elements
! are returned.
!
! loc_str is a list of element names.
! A space or a comma delimits the elements.
!
! An element name can be of the form
!   {key::}{branch>>}ele_id{##N}
! Where
!   key     = Optional key name ("quadrupole", "sbend", etc.)
!   branch  = Name or index of branch. May contain the wild cards "*" and "%".
!   ele_id  = Name or index of element. May contain the wild cards "*" and "%".
!               If a name and no branch is given, all branches are searched.
!               If an index and no branch is given, branch 0 is assumed.
!   ##N     = N^th instance of ele_id in the branch.
!
! An element range is of the form:
!   {key::}ele1:ele2{:step}
! Where:
!   key      = Optional key name ("quadrupole", "sbend", etc.)
!   ele1     = Starting element of the range.
!   ele2     = Ending element of the range. 
!   step     = Optional step increment Default is 1. 
! Note: ele1 and ele2 must be in the same branch.
!
! Examples:
!   "quad::x_br>>q*"   All quadrupoles of branch "x_br" whose name begins with "q"
!   "3,5:7"            Elements with index 3, 5, 6, and 7 in branch 0.
!   "2>>45:51"         Elements 45 through 51 of branch 2.
!   "q1:q5"            Elements between "q1" and "q5"
!   "marker::a*##2"    2^nd marker element in each branch whose name begins with "a".
! 
! Modules Needed:
!   use lat_ele_loc_mod
!
! Input:
!   loc_str  -- Character(*): Element name.
!   lat      -- lat_struct: Lattice to search through.
!   above_ubound_is_err
!            -- Logical, optional: If the upper bound "e2" on an "e1:e2" range construct 
!                 is above the maximum element index then treat this as an error? 
!                 Default is True. If False then set e2 to the maximum element index. 
!   ix_dflt_branch
!            -- Integer, optional: If present and not -1 then restrict search to specified branch.
!                 If not present or -1: Search all branches. Exception: For elements specified using 
!                 an integer index (EG: "43"), if ix_dflt_branch is not present or -1 use branch 0.
!
! Output:
!   eles(:) -- Ele_pointer_struct, allocatable: Array of matching elements.
!              Note: This routine does not try to deallocate eles.
!               It is up to you to deallocate eles if needed.
!   n_loc   -- Integer: Number of locations found.
!                Set to zero if no elements are found.
!   err     -- Logical: Set True if there is a decode error.
!                Note: not finding any element is not an error.
!-

subroutine lat_ele_locator (loc_str, lat, eles, n_loc, err, above_ubound_is_err, ix_dflt_branch)

implicit none

type (lat_struct) lat
type (ele_pointer_struct), allocatable :: eles(:)
type (ele_pointer_struct), allocatable, target :: eles2(:)
type (ele_struct), pointer :: ele_start, ele_end

character(*) loc_str
character(200) str
character(60) name, name2
character(1) delim
character(20) :: r_name = 'lat_ele_locator'

integer, optional :: ix_dflt_branch
integer i, j, ib, ios, ix, n_loc, n_loc2
integer in_range, step, ix_word, key

logical, optional :: above_ubound_is_err
logical err, err2, delim_found, above_ub_is_err

! init

err = .true.
n_loc = 0
str = loc_str
call str_upcase (str, str)
in_range = 0   ! 0 -> not in range construct, 1 -> read start, 2 -> read stop
step = 1
key = 0

! Loop over all items in the list

do

  ! Get next item. 
  ! If the split is in between "::" then need to piece the two haves together.

  call word_read (str, ':, ', name, ix_word, delim, delim_found, str)

  if (str(1:1) == ':') then
    call word_read (str(2:), ', ', name2, ix_word, delim, delim_found, str)
    ! If name2 contains a ":" then we have a "<key>::<ele1>:<ele2>" range construct.
    ix = index(name2, ':') 
    if (ix == 0) then  ! No range
      name = trim(name) // '::' // trim(name2)
    else
      key = key_name_to_key_index (name, .true.)
      if (key < 0) then
        call out_io (s_error$, r_name, 'ERROR: BAD ELEMENT TYPE: ' // name)
        return           ! Return on error
      endif
      name = name2(1:ix-1)
      str = trim(name2(ix+1:)) // delim // trim(str)
      delim = ':'
    endif
  endif

  if (name == '') exit

  ! Get list of elements for this item

  if (in_range == 2) then  ! Must be step
    read (name, *, iostat = ios) step
    if (ios /= 0) then
      call out_io (s_error$, r_name, 'ERROR: BAD STEP: ' // loc_str)
      return
    endif
    in_range = in_range + 1
  else
    above_ub_is_err = (logic_option(.true., above_ubound_is_err) .or. in_range /= 1)
    call lat_ele1_locator (name, lat, eles2, n_loc2, err2, above_ub_is_err, ix_dflt_branch)
    if (err2) return
  endif

  ! If not a range construct then just put this on the list

  if (delim /= ":" .and. in_range == 0) then
    if (.not. allocated(eles)) allocate(eles(n_loc+n_loc2))
    call re_allocate_eles(eles, n_loc+n_loc2, .true.)
    eles(n_loc+1:n_loc+n_loc2) = eles2(1:n_loc2)
    n_loc = n_loc + n_loc2
    cycle
  endif

  ! Handle range construct...
  ! Each item in the construct must represent only one element.

  if (n_loc2 /= 1) then
    call out_io (s_error$, r_name, 'MULTIPLE ELEMENTS ASSOCIATED WITH: ' // name, &
                                   'THIS IS NOT PERMITTED IN A RANGE CONSTRUCT: ' // loc_str)
    return
  endif

  if (in_range == 0) ele_start => eles2(1)%ele
  if (in_range == 1) ele_end   => eles2(1)%ele
  in_range = in_range + 1

  if (delim == ':') then
    if (in_range >= 3) then
      call out_io (s_error$, r_name, 'TOO MANY ":" IN RANGE CONSTRUCT: ' // loc_str)
      return
    endif
    cycle
  endif

  ! if we have the range then add it to the eles list.

  if (in_range > 1) then
    ele_end => eles2(1)%ele
    if (ele_start%ix_branch /= ele_end%ix_branch) then
      call out_io (s_error$, r_name, 'ERROR: ELEMENTS NOT OF THE SAME BRANCH IN RANGE: ' // loc_str)
      return
    endif
    ib = ele_start%ix_branch

    if (key > 0) then
      n_loc2 = 0
      do i = ele_start%ix_ele, ele_end%ix_ele, step
        if (lat%branch(ib)%ele(i)%key == key) n_loc2 = n_loc2 + 1
      enddo
    else
      n_loc2 = (ele_end%ix_ele - ele_start%ix_ele) / step + 1
    endif

    call re_allocate_eles(eles, n_loc+n_loc2, .true.)

    n_loc2 = 0
    do i = ele_start%ix_ele, ele_end%ix_ele, step
      if (key > 0 .and. lat%branch(ib)%ele(i)%key /= key) cycle
      n_loc2 = n_loc2 + 1
      eles(n_loc+n_loc2)%ele => lat%branch(ib)%ele(i)
    enddo

    n_loc = n_loc + n_loc2
    in_range = 0
    key = 0
    step = 1
  endif

enddo

! Check that no ranges are open.

if (in_range == 1) then
  call out_io (s_error$, r_name, 'MISSING ELEMENT AFTER ":" IN: ' // loc_str)
  return
endif

err = .false.

end subroutine lat_ele_locator

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine lat_ele1_locator (name, lat, eles, n_loc, err, above_ub_is_err, ix_dflt_branch)
!
! Routine to locate all the elements in a lattice that corresponds to loc_str. 
! Note: This routine is private and meant to be used only by lat_ele_locator
!-

subroutine lat_ele1_locator (name, lat, eles, n_loc, err, above_ub_is_err, ix_dflt_branch)

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_pointer_struct), allocatable :: eles(:)

character(*) name
character(20) :: r_name = 'lat_ele1_locator'

integer, optional :: ix_dflt_branch
integer i, k, ix, ix_branch, ixp, ios, ix_ele, n_loc
integer key, ix_dup, n_dup

logical err, do_match_wild
logical :: above_ub_is_err

! init

err = .true.
n_loc = 0

! key::name construct

ix = index(name, '::')
if (ix == 0) then
  key = 0
else
  if (name(:ix-1) == "*") then
    key = 0
  else
    key = key_name_to_key_index (name(:ix-1), .true.)
    if (key < 1) return
  endif
  name = name(ix+2:)
endif

! Look for "##N" suffix to indicate which instance to choose from when there
! are multiple elements of a given name.

ix = index(name, '##')
ix_dup = 0
if (ix /= 0) then
  if (.not. is_integer(name(ix+2:))) return
  read (name(ix+2:), *) ix_dup
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

ixp = index(name, '>>')
if (ixp > 0) then
  branch => pointer_to_branch (name(1:ixp-1), lat)
  if (.not. associated(branch)) then
    err = .false.
    return
  endif
  ix_branch = branch%ix_branch
  name = name(ixp+2:)
endif

! Read integer if present

if (is_integer(name)) then
  read (name, *, iostat = ios) ix_ele
  if (ios /= 0) then
    call out_io (s_error$, r_name, 'BAD ELEMENT LOCATION INDEX: ' // name)
    return
  endif
  if (.not. allocated(eles)) allocate (eles(1))
  if (ix_branch == -1) ix_branch = 0
  if (.not. above_ub_is_err .and. ix_ele > lat%branch(ix_branch)%n_ele_max) then
    ix_ele = lat%branch(ix_branch)%n_ele_max
  elseif (ix_ele < 0 .or. ix_ele > lat%branch(ix_branch)%n_ele_max) then
    call out_io (s_error$, r_name, 'ELEMENT LOCATION INDEX OUT OF RANGE: ' // name)
    return
  endif
  eles(1)%ele => lat%branch(ix_branch)%ele(ix_ele)
  n_loc = 1
  err = .false.
  return
endif

! We have a "class::name" construct if there is a "::" in the string
! or a wild card "*" or "%".

do_match_wild = .false.  
if (index(name, "*") /= 0 .or. index(name, "%") /= 0) do_match_wild = .true.

! search for matches

do k = lbound(lat%branch, 1), ubound(lat%branch, 1)
  n_dup = 0
  if (ix_branch /= -1 .and. k /= ix_branch) cycle
  do i = 0, lat%branch(k)%n_ele_max
    if (key /= 0 .and. lat%branch(k)%ele(i)%key /= key) cycle
    if (do_match_wild) then
      if (.not. match_wild(lat%branch(k)%ele(i)%name, name)) cycle
    else
      if (lat%branch(k)%ele(i)%name /= name) cycle
    endif
    if (ix_dup > 0) then
      n_dup = n_dup + 1
      if (n_dup > ix_dup) exit 
      if (n_dup /= ix_dup) cycle
    endif
    n_loc = n_loc + 1
    if (.not. allocated(eles) .or. size(eles) < n_loc) call re_allocate_eles (eles, 2*n_loc, .true.)
    eles(n_loc)%ele => lat%branch(k)%ele(i)
  enddo
enddo 

err = .false.

end subroutine lat_ele1_locator

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine re_allocate_eles (eles, n, save_old, exact)
!
! Routine to allocate an array of ele_pointer_structs.
!
! Modules Needed:
!   use lat_ele_loc_mod
!
! Input:
!   eles(:) -- ele_pointer_struct, allocatable: Array of element pointers with possible old data.
!   n        -- Integer: Array size to set.
!   save_old -- Logical, optional: If present and True then save the old data.
!   exact    -- Logical, optional: If present and True then eles will have size = n
!                 If False (default), reallcation will not be done if eles is already large enough
!
! Output:
!   eles(:) -- ele_pointer_struct, allocatable: Array of element pointers.
!-

subroutine re_allocate_eles (eles, n, save_old, exact)

implicit none

type (ele_pointer_struct), allocatable :: eles(:)
type (ele_pointer_struct), allocatable :: l_temp(:)
integer n, n_old, nn
logical, optional :: save_old, exact

!

if (.not. allocated(eles)) then
  allocate (eles(n))
  return
endif

if (.not. logic_option(.false., exact) .and. size(eles) >= n) return
if  (size(eles) == n) return

if (logic_option (.false., save_old)) then
  call move_alloc (eles, l_temp)
  allocate (eles(n))
  nn = min(size(l_temp), n)
  eles(1:nn) = l_temp(1:nn)
  deallocate (l_temp)
else
  deallocate (eles)
  allocate (eles(n))
endif

end subroutine re_allocate_eles

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function ele_to_lat_loc (ele) result (ele_loc)
!
! Function to return an lat_ele_loc_struct identifying where an element is in the lattice.
!
! Input:
!   ele -- Ele_struct: Element to be identified
!
! Output:
!   ele_loc -- Lat_ele_loc_struct: Element identifier.
!-

function ele_to_lat_loc (ele) result (ele_loc)

type (ele_struct) ele
type (lat_ele_loc_struct) ele_loc

!

ele_loc%ix_ele = ele%ix_ele
ele_loc%ix_branch = ele%ix_branch

end function ele_to_lat_loc

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine get_slave_list (lord, slave_list, n_slave)
!
! Subroutine to get the list of slaves for a lord element.
!
! This is a list of ultimate slaves. That is, slaves in the tracking part 
! of the lattice. Thus if the element lord controls an
! element lord1 which controlls an element lord2, then lord2 will
! show up in the slave_list but lord1 will not.
!
! If the "lord" element does not have any slaves, 
! then the slave_list will just be the lord element.
!
! This routine will increase the size of slave_list if needed but will
! not decrease it.
!
! Modules needed:
!   use lat_ele_loc_mod
!
! Input:
!   lord  -- Ele_struct: The lord element.
!
! Output:
!   slaves(:) -- Ele_pointer_struct, allocatable :: Array of slaves.
!   n_slave   -- Integer: Number of slaves.
!-

subroutine get_slave_list (lord, slaves, n_slave)

implicit none

type (ele_struct) :: lord
type (ele_pointer_struct), allocatable :: slaves(:)

integer n_slave

!

n_slave = 0
if (.not. allocated(slaves)) call re_allocate_eles (slaves, lord%n_slave)

call get_slaves (lord)

!--------------------------------------------------------------------------
contains

recursive subroutine get_slaves (lord)

type (ele_struct) lord
type (ele_struct), pointer :: slave
integer i, ix

!

do i = 1, lord%n_slave
  slave => pointer_to_slave(lord, i)
  if (slave%n_slave > 0) then
    call get_slaves (slave)
  else
    n_slave = n_slave + 1
    if (n_slave > size(slaves)) call re_allocate_eles(slaves, n_slave + 4, .true.)
    slaves(n_slave)%ele => slave
  endif
enddo

end subroutine get_slaves

end subroutine get_slave_list

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine get_field_ele_list (ele, field_eles, dz_offset, n_field_ele)
!
! Subroutine to get the list of elements that specify the field for the given element.
!
! This is a list of elements that possibly have field info for a given element. 
! For example: a slice_slave under a super_slave under a multipass_slave. In this case
! the multipass_lord will store the field info.
! Note: groups, overlays, and girders will never have field info.
!
! If the given element does not have any lords, 
! the element list will just consist of the element itself.
!
! This routine will increase the size of field_ele_list if needed but will
! not decrease it.
!
! Modules needed:
!   use lat_ele_loc_mod
!
! Input:
!   lat   -- lat_struct: Lattice
!   ele   -- Ele_struct: Element whose fields are to be evaluated.
!
! Output:
!   field_eles(:) -- Ele_pointer_struct, allocatable :: Array of field_eles.
!   dz_offset(:)  -- Offsets of ele from field elements.
!   n_field_ele   -- Integer: Number of field_eles.
!-

subroutine get_field_ele_list (ele, field_eles, dz_offset, n_field_ele)

implicit none

type (ele_struct), target :: ele
type (ele_pointer_struct), allocatable :: field_eles(:)

real(rp), allocatable :: dz_offset(:)
real(rp) offset
integer n_field_ele

!

n_field_ele = 0
if (.not. allocated(field_eles)) call re_allocate_eles (field_eles, 3)
allocate (dz_offset(3))

offset = 0
call get_field_eles (ele, offset)

!--------------------------------------------------------------------------
contains

recursive subroutine get_field_eles (this_ele, thiz_offset)

type (ele_struct), target :: this_ele
type (ele_struct), pointer :: field_ele
real(rp) thiz_offset, thiz_offset2
integer i, ix

!

if (this_ele%field_calc == refer_to_lords$) then
  do i = 1, this_ele%n_lord
    field_ele => pointer_to_lord(this_ele, i)

    select case (field_ele%key)
    case (overlay$, group$, girder$); cycle
    end select

    if (field_ele%lord_status == multipass_lord$) cycle

    thiz_offset2 = thiz_offset + (this_ele%s - this_ele%value(l$)) - (field_ele%s - field_ele%value(l$))

    call get_field_eles (field_ele, thiz_offset2)

  enddo

else
  n_field_ele = n_field_ele + 1
  if (n_field_ele > size(field_eles) .or. n_field_ele > size(dz_offset)) then
    call re_allocate_eles(field_eles, n_field_ele + 10, .true.)
    call re_allocate (dz_offset, n_field_ele + 10, .false.)
  endif
  field_eles(n_field_ele)%ele => this_ele
  dz_offset(n_field_ele) = thiz_offset
endif

end subroutine get_field_eles

end subroutine get_field_ele_list

end module
