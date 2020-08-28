!+
! Subroutine tao_locate_all_elements (ele_list, eles, err, ignore_blank) 
!
! Subroutine to find the lattice elements in the lattice corresponding to the ele_list argument. 
!
! See also tao_locate_elements for a routine that only searches one given universe.
!
! Input:
!   ele_list     -- Character(*): String with element names using element list format.
!   ignore_blank -- Logical, optional: If present and true then do nothing if
!     ele_list is blank. otherwise treated as an error.
!
! Output:
!   eles  -- ele_pointer_struct(:), allocatable: Array of elements in the model lat. 
!     %id  -- Set to universe number.
!   err   -- Logical: Set true on error.
!-

subroutine tao_locate_all_elements (ele_list, eles, err, ignore_blank)

use tao_interface, except_dummy => tao_locate_all_elements

implicit none

type (tao_universe_struct), pointer :: u
type (ele_pointer_struct), allocatable :: eles(:)
type (ele_pointer_struct), allocatable :: this_eles(:)

integer i, ix, ix_word, n_eles, n0

character(*) ele_list
character(100) ele_name, word, string
character(*), parameter :: r_name = 'tao_locate_all_elements'
character(1) delim

logical err, delim_found
logical, optional :: ignore_blank
logical, allocatable :: picked(:)

!

err = .true.

call re_allocate_eles (eles, 0, exact = .true.)

call str_upcase (ele_name, ele_list)
call string_trim (ele_name, ele_name, ix)

if (ix == 0 .and. logic_option(.false., ignore_blank)) return

if (ix == 0) then
  call out_io (s_error$, r_name, 'ELEMENT NAME IS BLANK')
  return
endif

! Loop over all items in the element list

string = ele_list
do 
  call tao_pick_universe (string, string, picked, err);  if (err) return
  call word_read (string, ', ', word, ix_word, delim, delim_found, string)
  if (ix_word == 0) exit
  do i = lbound(s%u, 1), ubound(s%u, 1)
    if (.not. picked(i)) cycle
    call lat_ele_locator (word, s%u(i)%model%lat, this_eles, n_eles, err)
    if (err) return
    if (n_eles == 0) cycle
    n0 = size(eles)
    call re_allocate_eles (eles, n0+n_eles, .true., .true.)
    eles(n0+1:n0+n_eles) = this_eles(1:n_eles)
    eles(n0+1:n0+n_eles)%id = i
  enddo
enddo

if (size(eles) == 0) then
  call out_io (s_error$, r_name, 'ELEMENT(S) NOT FOUND: ' // ele_list)
  err = .true.
  return
endif

end subroutine tao_locate_all_elements

