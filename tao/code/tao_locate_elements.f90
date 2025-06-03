!+
! Subroutine tao_locate_elements (ele_list, ix_universe, eles, err, lat_type, ignore_blank, 
!                                                           err_stat_level, ix_branch, multiple_eles_is_err) 
!
! Subroutine to find the lattice elements in the lattice corresponding to the ele_list argument. 
!
! Note: If ele_list can contain elements from different universes, use
! the routine:
!   tao_locate_all_elements
! The disadvantage of tao_locate_all_elements is that the ele_list cannot use Bmad constructs like 
! ranges (EG "q1:q2" to designate all elements between q1 and q2).
!
! Input:
!   ele_list       -- character(*): String with element names using element list format.
!   ix_universe    -- integer: Universe to search. -1 => search s%global%default_universe. -2 (all unis) => error.
!                       ix_universe is ignored if ele_list starts with a universe specifier "N@".
!   lat_type       -- integer, optional: model$ (default), design$, or base$.
!   ignore_blank   -- logical, optional: If present and true then do nothing if
!                     ele_list is blank. otherwise treated as an error.
!   err_stat_level -- integer, optional: Status level for error messages. If not present,
!                       print with level s_error$. Use s_nooutput$ to prevent printing.
!   ix_branch      -- integer, optional: If present and non-negative then use this as the branch index 
!                       for elements specified using an integer index (EG: "43").
!                       If -1 use the default branch, search all branches.
!   multiple_eles_is_err
!                  -- logical, optional: If present and True then matching to more than one element is an error.
!
! Output:
!   eles(:)        -- ele_pointer_struct(:), allocatable: Array of elements in the model lat. 
!   err            -- logical: Set true on error.
!-

subroutine tao_locate_elements (ele_list, ix_universe, eles, err, lat_type, ignore_blank, &
                                           err_stat_level, above_ubound_is_err, ix_branch, multiple_eles_is_err)

use tao_interface, dummy => tao_locate_elements

implicit none

type (tao_universe_struct), pointer :: u
type (tao_lattice_struct), pointer :: tao_lat
type (ele_pointer_struct), allocatable :: eles(:)

integer, optional :: lat_type, ix_branch, err_stat_level
integer ios, ix, ix_universe, num, i, i_ix_ele, n_loc, n_loc_old, print_lev, ixb

character(*) ele_list
character(len(ele_list)) ele_lst
character(*), parameter :: r_name = 'tao_locate_elements'

logical err, err_flag
logical, allocatable :: picked(:)
logical, optional :: ignore_blank, above_ubound_is_err, multiple_eles_is_err

! 

err = .true.
print_lev = integer_option(s_error$, err_stat_level)

call re_allocate_eles (eles, 0, exact = .true.)

call str_upcase (ele_lst, ele_list)
ele_lst = unquote(ele_lst)
call string_trim (ele_lst, ele_lst, ix)

if (ix == 0 .and. logic_option(.false., ignore_blank)) then
  err = .false.
  return
endif

if (ix == 0) then
  call out_io (print_lev, r_name, 'ELEMENT NAME IS BLANK')
  return
endif

!

ix = index(ele_lst, '@')
if (ix == 0) then
  u => tao_pointer_to_universe (ix_universe)
  if (.not. associated(u)) return
  allocate(picked(ubound(s%u, 1)))
  picked = .false.
  picked(u%ix_uni) = .true.
else
  call tao_pick_universe(ele_lst, ele_lst, picked, err_flag)
  if (err_flag) return
endif

n_loc = 0
n_loc_old = 0
do i = 1, ubound(picked,1)
  if (.not. picked(i)) cycle
  u => s%u(i)
  tao_lat => tao_pointer_to_tao_lat (u, lat_type)

  if (present(ix_branch)) then
    call lat_ele_locator (ele_lst, tao_lat%lat, eles, n_loc, err, above_ubound_is_err, tao_branch_index(ix_branch), append_eles = .true.)
  else
    call lat_ele_locator (ele_lst, tao_lat%lat, eles, n_loc, err, above_ubound_is_err, append_eles = .true.)
  endif

  if (err) return
  eles(n_loc_old+1:n_loc)%id = i
  n_loc_old = n_loc
enddo

!

if (logic_option(.false., multiple_eles_is_err) .and. n_loc > 1) then
  call out_io (print_lev, r_name, 'MULTIPLE ELEMENTS FOUND MATCHING: ' // ele_list)
  err = .true.
  return
endif

if (n_loc == 0) then
  call out_io (print_lev, r_name, 'ELEMENT(S) NOT FOUND: ' // ele_list)
  err = .true.
  return
endif

call re_allocate_eles (eles, n_loc, .true., .true.)

end subroutine tao_locate_elements
