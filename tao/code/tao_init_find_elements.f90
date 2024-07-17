!+
! subroutine tao_init_find_elements (u, search_string, eles, attribute, found_one)
!
! Routine to search the lattice for the specified element and flags eles(:)
! This routine is meant to be used by initialization routines.
! Otherwise use tao_locate_elements instead.
!
! Input:
!   u             -- tao_universe_struct: Universe to search
!   search_string -- character(*): What to search for
!   attribute     -- character(*), optional: Check that attribute of element is free to vary.
!
! Output:
!   eles(:)       -- ele_pointer_struct, allocatable: List of matching elements.
!                     Size is zero if no elements found.
!   found_one     -- logical, optional: Set True if a matching element is found. 
!                     However: Not set if no matching element found.
!-


subroutine tao_init_find_elements (u, search_string, eles, attribute, found_one)

use tao_interface, dummy => tao_init_find_elements
use attribute_mod, only: attribute_index, attribute_free

implicit none

type (tao_universe_struct), target :: u
type (ele_struct), pointer :: ele
type (ele_pointer_struct), allocatable :: eles(:)

character(*) search_string
character(*), optional :: attribute
character(len(search_string)) string
character(40) ele_name, key_name_in
character(*), parameter :: r_name = 'tao_init_find_elements'

integer key, found_key, ix_attrib, t
integer i, k, ix, ii, j, ix0, ix1, ix2, n_ele

logical, optional :: found_one
logical no_slaves, no_lords, err, warn_given

! Sort switches

no_slaves = .false.
no_lords = .false.

call string_trim(search_string, string, ix)

do
  if (string(1:1) /= '-') exit

  select case (string(1:ix))
  case ('-no_lords') 
    no_lords = .true.
  case ('-no_slaves') 
    no_slaves = .true.
  case default
    call out_io (s_warn$, r_name, "BAD SEARCH SWITCH: " // search_string)
    call re_allocate_eles(eles, 0, exact = .true.)
    return
  end select

  call string_trim (string(ix+1:), string, ix)
enddo

! Find elements

call tao_locate_elements (string, u%ix_uni, eles, err, err_stat_level = s_nooutput$)
if (size(eles) > 0 .and. present(found_one)) found_one = .true.

warn_given = .false.
n_ele = 0
do j = 1, size(eles)
  ele => eles(j)%ele
  t = ele%slave_status
  if ((t == multipass_slave$ .or. t == super_slave$) .and. no_slaves) cycle
  t = ele%lord_status 
  if ((t == girder_lord$ .or. t == overlay_lord$ .or. t == group_lord$ .or. &
       t == super_lord$ .or. t == multipass_lord$) .and. no_lords) cycle
  ! If attribute is not free then don't count it
  if (present(attribute)) then
    ix_attrib = attribute_index(ele, attribute)
    if (ix_attrib < 1) then
      call out_io (s_error$, r_name, &
                      'BAD ATTRIBUTE: ' // attribute, &
                      'FOR VARIABLE SEARCH: ' // search_string, &
                      'FOR ELEMENT: ' // ele%name)
      return
    endif
    if (.not. attribute_free (eles(j)%ele, attribute, .false.)) then
      if (.not. warn_given) then
        call out_io (s_error$, r_name, &
                  'NON-FREE ATTRIBUTE ' // attribute, &
                  'FOR VARIABLE SEARCH: ' // search_string, &
                  'FOR ELEMENT: ' // ele%name)
        warn_given = .true.
      endif
      cycle
    endif
  endif
  ! Passes test so add it to the list
  n_ele = n_ele + 1
  eles(n_ele)%ele => eles(j)%ele
enddo

call re_allocate_eles (eles, n_ele, .true., .true.)

end subroutine tao_init_find_elements
