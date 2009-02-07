!+
! Subroutine elements_locator (ele_name, lat, indx, err)
!
! Subroutine to locate all the elements of a certain key 
! and certain name. ele_name is of the form:
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
!   "3,5:7"      Elements with index 3, 5, 6, and 7.
! 
! Modules Needed:
!   use bmad
!
! Input:
!   ele_name  -- Character(*): Element name.
!   lat       -- lat_struct: Lattice to search through.
!
! Output:
!   indx(:) -- Integer, allocatable: Array of indexes of the elements found.
!              The size of indx will be zero if no elements are found.
!              Note: This routine does not try to deallocate indx.
!               It is up to you to deallocate indx if needed.
!   err     -- Logical: Set True if there is a decode error.
!-

subroutine elements_locator (ele_name, lat, indx, err)

use bmad_struct
use bmad_interface, except_dummy => elements_locator

implicit none

type (lat_struct) lat

character(*) ele_name
character(40) name
integer, allocatable :: indx(:)
integer i, j, ix, num, key
logical err

! index array

err = .true.

if (is_integer(ele_name(1:1))) then
  call location_decode (ele_name, lat%ele(:)%bmad_logic, 0, num)
  if (num == -1) return

! key:name construct

else

  ix = index(ele_name, ':')
  if (ix == 0) then
    key = 0
    call str_upcase (name, ele_name)
  else
    key = key_name_to_key_index (ele_name(:ix-1), .true.)
    if (key < 1) return
    call str_upcase (name, ele_name(ix+1:))
  endif

  num = 0
  lat%ele(:)%bmad_logic = .false.
  do i = 0, lat%n_ele_max
    if (key /= 0 .and. lat%ele(i)%key /= key) cycle
    if (.not. match_wild(lat%ele(i)%name, name)) cycle
    lat%ele(i)%bmad_logic = .true.
    num = num + 1
  enddo

endif

! Fill in indx array

call re_allocate (indx, num)
j = 0
do i = 0, lat%n_ele_max
  if (.not. lat%ele(i)%bmad_logic) cycle
  j = j + 1
  indx(j) = i
enddo

err = .false.

end subroutine
