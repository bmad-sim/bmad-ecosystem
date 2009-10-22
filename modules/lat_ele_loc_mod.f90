module lat_ele_loc_mod

use bmad_struct
use bmad_interface

contains

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine lat_ele_locator (loc_str, lat, eles, n_loc, err)
!
! Routine to locate all the elements of a certain key 
! and certain name. loc_str is of the form:
!   <key>:<name> 
! or 
!   <name>
! or
!   <indexes>
! Where
!   <key>     = key name ("quadrupole", "sbend", etc.)
!   <name>    = Name of element. May contain the wild cards "*" and "%".
!   <indexes> = List of indexes in the lat%ele list.
! Example:
!   "quad:q*"    All quadrupoles whose name begins with "q"
!   "*:*"        All elements.
!   "3,5:7"      Elements with index 3, 5, 6, and 7 in branch 0.
!   "1.45:51"    Elements 45 through 51 of branch 1.
! 
! Modules Needed:
!   use lat_ele_loc_mod
!
! Input:
!   loc_str  -- Character(*): Element name.
!   lat      -- lat_struct: Lattice to search through.
!
! Output:
!   eles(:) -- Ele_pointer_struct, allocatable: Array of matching elements.
!              Note: This routine does not try to deallocate eles.
!               It is up to you to deallocate eles if needed.
!   n_loc   -- Integer: Number of locations found.
!                Set to zero if no elements are found.
!   err     -- Logical: Set True if there is a decode error.
!-

subroutine lat_ele_locator (loc_str, lat, eles, n_loc, err)

implicit none

type (lat_struct) lat

character(*) loc_str
character(40) name
type (ele_pointer_struct), allocatable :: eles(:)
integer i, j, k, ix, key, n_loc
logical err, do_match_wild

! index array

err = .true.
n_loc = 0

if (is_integer(loc_str(1:1))) then
  call lat_location_decode (loc_str, lat, eles, err)
  n_loc = size(eles)
  return
endif

! key:name construct

ix = index(loc_str, ':')
if (ix == 0) then
  key = 0
  call str_upcase (name, loc_str)
else
  if (loc_str(:ix-1) == "*") then
    key = 0
  else
    key = key_name_to_key_index (loc_str(:ix-1), .true.)
    if (key < 1) return
  endif
  call str_upcase (name, loc_str(ix+1:))
endif

! Save time by deciding if we need to call match_wild or not.

do_match_wild = .false.  
if (index(name, "*") /= 0 .or. index(name, "%") /= 0) do_match_wild = .true.

! search for matches

n_loc = 0
do k = lbound(lat%branch, 1), ubound(lat%branch, 1)
  do i = 0, lat%branch(k)%n_ele_max
    if (key /= 0 .and. lat%branch(k)%ele(i)%key /= key) cycle
    if (do_match_wild) then
      if (.not. match_wild(lat%branch(k)%ele(i)%name, name)) cycle
    else
      if (lat%branch(k)%ele(i)%name /= name) cycle
    endif
    n_loc = n_loc + 1
    if (.not. allocated(eles)) call re_allocate_eles (eles, 1, .true.)
    if (size(eles) < n_loc) call re_allocate_eles (eles, 2*n_loc, .true.)
    eles(n_loc)%ele => lat%branch(k)%ele(i)
  enddo
enddo

err = .false.

end subroutine lat_ele_locator

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine lat_location_decode (loc_str, lat, eles, err)
!
! Routine to parse a numeric list of element locations.
!
! Input:
!   loc_str  -- Character*(*): Array of locations.
!                     "," or a space delimits location numbers 
!                     A ":" is used for a range of locations. 
!                     A second ":" specifies a step increment.
!   lat      -- lat_struct: Lattice to search through.
!
! Modules Needed:
!   use lat_ele_loc_mod
!
! Output:
!   eles(:) -- Ele_pointer_struct: Array of matching locations.
!   err     -- Logical: Set True if location does not correspond to a physical element.
!
! Example:
!     string = '3:37 98, 101:115:2'
! Then:
!     call location_decode (string, eles, err)
! Gives:
!     eles(:) points to 3 to 37, 98, and odd numbers between 101 to 115
!-

subroutine lat_location_decode(loc_str, lat, eles, err)

implicit none

type (lat_struct), target :: lat
type (ele_pointer_struct), allocatable :: eles(:)
type (branch_struct), pointer :: branch

integer i, j, k
integer n_loc, ios, ix_next, ix_branch, step,start_loc, end_loc, ix_word
integer ix_ele, ixp

character(*) loc_str
character(1) delim
character(len(loc_str)) str
character :: r_name = 'lat_location_decode'

logical err
logical range_found, step_found

! initialize eles

str = loc_str

err = .true.
range_found = .false.
step_found = .false.
ix_next = 1
ix_branch = 0
branch => lat%branch(0)
n_loc = 0
step = 1

do k = lbound(lat%branch, 1), ubound(lat%branch, 1)
  lat%branch(k)%ele(:)%bmad_logic = .false.
enddo

!

do

  call string_trim2 (str(ix_next:), ':,', str, ix_word, delim, ix_next)

  if (ix_word == 0) exit

  ! Look for a name match

  ix_ele = 0

  ! If there is no name match then assume it is a number

  if (ix_ele == 0) then
    ixp = index(str(:ix_word), '.')
    if (ixp /= 0 .and. ixp <= ix_word) then
      read (str(1:ixp), *, iostat = ios) ix_branch
      if (ios /= 0) then
        call out_io (s_error$, r_name, 'ERROR: BAD LOCATION: ' // str(:ix_word))
        return
      endif
      if (ix_branch < 0 .or. ix_branch > ubound(lat%branch, 1)) then
        call out_io (s_error$, r_name, 'ERROR: BRANCH INDEX OUT OF RANGE: ' // str(:ix_word))
        return
      endif
      branch => lat%branch(ix_branch)
      str = str(ixp+1:)
      ix_word = ix_word - ixp
    endif
    read (str(:ix_word), *, iostat = ios) ix_ele
    if (ios /= 0) then
      call out_io (s_error$, r_name, 'ERROR: BAD LOCATION: ' // str(:ix_word))
      return
    endif
  endif

  if (step_found) step = ix_ele
    
  if (delim == ':') then
    if (range_found) then
      if (step_found) then
        call out_io (s_error$, r_name, 'ERROR: BAD RANGE ' // str(1:20))
        return
      else
        step_found = .true.
        end_loc = ix_ele
      endif
    else
      range_found = .true.
      start_loc = ix_ele
    endif
  else
    if (range_found) then
      branch%ele(start_loc:end_loc:step)%bmad_logic = .true.
      n_loc = n_loc + (end_loc - start_loc) / step
      range_found = .false.
      step_found = .false.
      step = 1
    else
      if (ix_ele < 0 .or. ix_ele > ubound(branch%ele, 1)) then
        call out_io (s_error$, r_name, 'ERROR: ELEMENT INDEX OUT OF RANGE: ' // str(:ix_word))
        return
      endif
      branch%ele(ix_ele)%bmad_logic = .true.
      n_loc = n_loc + 1
    endif
  endif

  if (ix_next == 0) exit

enddo

!--------

if (range_found) then
  call out_io (s_error$, r_name, 'ERROR IN LOCATION_DECODE: OPEN RANGE')
  return
endif

! count number of elements in arrray

call re_allocate_eles (eles, n_loc)
j = 0
do k = lbound(lat%branch, 1), ubound(lat%branch, 1)
  do i = 0, lat%branch(k)%n_ele_max
    if (.not. lat%branch(k)%ele(i)%bmad_logic) cycle
    j = j + 1
    eles(j)%ele => lat%branch(k)%ele(i)
  enddo
enddo

err = .false.

end subroutine lat_location_decode

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine re_allocate_eles (eles, n, save)
!
! Routine to allocate an array of ele_pointer_structs.
!
! Modules Needed:
!   use lat_ele_loc_mod
!
! Input:
!   n -- Integer: Array size to set.
!   save -- Logical, optional: If present and True then save the old data.
!
! Output:
!   eles(:) -- ele_pointer_struct, allocatable: Array of element pointers.
!-

subroutine re_allocate_eles (eles, n, save)

implicit none

type (ele_pointer_struct), allocatable :: eles(:)
type (ele_pointer_struct), allocatable :: l_temp(:)
integer n, n_old
logical, optional :: save

!

if (.not. allocated(eles)) then
  allocate (eles(n))
  return
endif

if  (size(eles) == n) return

if (logic_option (.false., save)) then
  n_old = min(size(eles), n)
  allocate (l_temp(n_old))
  l_temp = eles(1:n_old)
endif

deallocate (eles)
allocate (eles(n))

if (logic_option (.false., save)) then
  eles(1:n_old) = l_temp
  deallocate (l_temp)
endif

end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function pointer_to_ele (lat, ix_ele, ix_branch) result (ele_ptr)
!
! Function to return a pointer to an element in a lattice.
! See pointer_to_ele for more details.
!-

function pointer_to_ele (lat, ix_ele, ix_branch) result (ele_ptr)

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele_ptr

integer ix_branch, ix_ele

!

ele_ptr => null()

if (ix_branch < 0 .or. ix_branch > ubound(lat%branch, 1)) return
if (ix_ele < 0 .or. ix_ele > lat%branch(ix_branch)%n_ele_max) return

ele_ptr => lat%branch(ix_branch)%ele(ix_ele)

end function

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function pointer_to_slave (lat, lord, ix_slave, ix_contrl) result (slave_ptr)
!
! Function to point to a slave of a lord.
!
! Modules Needed:
!   use lat_ele_loc_mod
!
! Input:
!   lat      -- lat_struct: Lattice containing the lord
!   lord     -- Ele_struct: Pointer to the lord element
!   ix_slave -- Integer: Index of the slave. ix_slave goes from 1 to lord%n_lord
!
! Output:
!   slave_ptr  -- Ele_struct, pointer: Pointer to the slave.
!   ix_control -- Integer, optional :: index of appropriate lat%control(:) element.
!-

function pointer_to_slave (lat, lord, ix_slave, ix_control) result (slave_ptr)

implicit none

type (lat_struct), target :: lat
type (ele_struct) lord
type (ele_struct), pointer :: slave_ptr
type (control_struct), pointer :: con

integer, optional :: ix_control
integer ix_slave, icon

!

icon = lord%ix1_slave + ix_slave - 1
con => lat%control(icon)
slave_ptr => lat%branch(con%ix_branch)%ele(con%ix_slave)
if (present(ix_control)) ix_control = icon

end function

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function pointer_to_lord (lat, slave, ix_lord, ix_control) result (lord_ptr)
!
! Function to point to a lord of a slave.
!
! Modules Needed:
!   use lat_ele_loc_mod
!
! Input:
!   lat        -- lat_struct: Lattice containing the lord
!   slave      -- Ele_struct: Slave element.
!   ix_lord    -- Integer: Index of the lord. ix_lord goes from 1 to slave%n_lord
!
! Output:
!   lord_ptr   -- Ele_struct, pointer: Pointer to the lord.
!   ix_control -- Integer, optional :: index of appropriate lat%control(:) element.
!-

function pointer_to_lord (lat, slave, ix_lord, ix_control) result (lord_ptr)

implicit none

type (lat_struct), target :: lat
type (ele_struct) slave
type (ele_struct), pointer :: lord_ptr
type (control_struct), pointer :: con

integer, optional :: ix_control
integer ix_lord, icon

!

icon = lat%ic(slave%ic1_lord + ix_lord - 1)
lord_ptr => lat%ele(lat%control(icon)%ix_lord)
if (present(ix_control)) ix_control = icon

end function

end module
