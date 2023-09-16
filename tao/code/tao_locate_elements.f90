!+
! Subroutine tao_locate_elements (ele_list, ix_universe, eles, err, lat_type, ignore_blank, 
!                                                           print_err, ix_dflt_branch, multiple_eles_is_err) 
!
! Subroutine to find the lattice elements in the lattice
! corresponding to the ele_list argument. 
!
! Note: If ele_list can contain elements from different universes, use
! the routine:
!   tao_locate_all_elements 
!
! Input:
!   ele_list       -- character(*): String with element names using element list format.
!   ix_universe    -- integer: Universe to search. -1 => search s%global%default_universe. -2 (all unis) => error.
!   lat_type       -- integer, optional: model$ (default), design$, or base$.
!   ignore_blank   -- logical, optional: If present and true then do nothing if
!                     ele_list is blank. otherwise treated as an error.
!   print_err      -- integer, optional: Status level for error messages. If not present,
!                       print with level s_error$. Use s_nooutput$ to prevent printing.
!   ix_dflt_branch -- integer, optional: If present and positive then use this as the branch index 
!                       for elements specified using an integer index (EG: "43").
!                       If not present or -1 the default branch is branch 0.
!   multiple_eles_is_err
!                  -- logical, optional: If present and True then matching to more than one element is an error.
!
! Output:
!   eles  -- ele_pointer_struct(:), allocatable: Array of elements in the model lat. 
!   err   -- logical: Set true on error.
!-

subroutine tao_locate_elements (ele_list, ix_universe, eles, err, lat_type, ignore_blank, &
                                           print_err, above_ubound_is_err, ix_dflt_branch, multiple_eles_is_err)

use tao_interface, dummy => tao_locate_elements

implicit none

type (tao_universe_struct), pointer :: u
type (tao_lattice_struct), pointer :: tao_lat
type (ele_pointer_struct), allocatable :: eles(:)

integer, optional :: lat_type, ix_dflt_branch, print_err
integer ios, ix, ix_universe, num, i, i_ix_ele, n_loc, print_lev

character(*) ele_list
character(len(ele_list)) ele_name
character(*), parameter :: r_name = 'tao_locate_elements'

logical err
logical, optional :: ignore_blank, above_ubound_is_err, multiple_eles_is_err

! 

err = .true.
print_lev = integer_option(s_error$, print_err)

call re_allocate_eles (eles, 0, exact = .true.)

call str_upcase (ele_name, ele_list)
call string_trim (ele_name, ele_name, ix)

if (ix == 0 .and. logic_option(.false., ignore_blank)) then
  err = .false.
  return
endif

if (ix == 0) then
  call out_io (print_lev, r_name, 'ELEMENT NAME IS BLANK')
  return
endif

u => tao_pointer_to_universe (ix_universe)
if (.not. associated(u)) return

tao_lat => tao_pointer_to_tao_lat (u, lat_type)

call lat_ele_locator (ele_name, tao_lat%lat, eles, n_loc, err, above_ubound_is_err, ix_dflt_branch)
if (err) return

if (n_loc == 0) then
  call out_io (print_lev, r_name, 'ELEMENT(S) NOT FOUND: ' // ele_list)
  err = .true.
  return
endif

if (logic_option(.false., multiple_eles_is_err) .and. n_loc > 1) then
  call out_io (print_lev, r_name, 'MULTIPLE ELEMENTS FOUND MATCHING: ' // ele_list)
  err = .true.
  return
endif


call re_allocate_eles (eles, n_loc, .true., .true.)

end subroutine tao_locate_elements
