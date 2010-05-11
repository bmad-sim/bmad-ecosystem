!+
! Subroutine pointers_to_attribute (lat, ele_name, attrib_name, do_allocation,
!                     ptr_array, err_flag, err_print_flag, eles, ix_attrib)
!
! Returns an array of pointers to an attribute with name attrib_name within 
! elements with name ele_name.
! Note: ele_name = 'BEAM_START' corresponds to the lat%beam_start substructure. 
! Note: ele_name can be a list of element indices. For example:
!           ele_name = "3:5"
!  This sets elements 3, 4, and 5 in the lat%ele(:) array.
! Note: ele_name can be in key:name format and include the wild card characters "*" and "%".
!       For example: "quad:q*"
! Note: Use attribute_free to see if the attribute may be varied independently.
! Note: When using wild cards, it is *not* an error if some of the matched elements do 
!       not have the the required attribute as long as at least one does.
! Modules needed:
!   use bmad
!
! Input:
!   lat             -- lat_struct: Lattice.
!   ele_name        -- Character(*): Element name. Must be uppercase
!   attrib_name     -- Character(*): Attribute name. Must be uppercase.
!                       For example: "HKICK".
!   do_allocation   -- Logical: If True then do an allocation if needed.
!                       EG: The multipole An and Bn arrays need to be allocated
!                       before their use.
!   err_print_flag  -- Logical, optional: If present and False then suppress
!                       printing of an error message on error.
!
! Output:
!   ptr_array(:) -- Real_pointer_struct, allocatable: Pointer to the attribute.
!                     Pointer will be deassociated if there is a problem.
!   err_flag     -- Logical: Set True if attribtute not found or attriubte
!                     cannot be changed directly.
!   eles(:)      -- Ele_pointer_struct, optional, allocatable: Array of element pointers.
!   ix_attrib    -- Integer, optional: If applicable then this is the index to the 
!                     attribute in the ele%value(:) array.
!-

#include "CESR_platform.inc"

Subroutine pointers_to_attribute (lat, ele_name, attrib_name, do_allocation, &
                        ptr_array, err_flag, err_print_flag, eles, ix_attrib)

use bmad_struct
use bmad_interface, except_dummy => pointers_to_attribute
use lat_ele_loc_mod

implicit none

type (lat_struct), target :: lat
type (ele_struct), target :: beam_start
type (real_pointer_struct), allocatable :: ptr_array(:)
type (real_pointer_struct), allocatable, save :: ptrs(:)
type (ele_pointer_struct), optional, allocatable :: eles(:)
type (ele_pointer_struct), allocatable, save :: eles2(:)

integer, optional :: ix_attrib
integer n, i, ix, key, ix_a, n_loc

character(*) ele_name
character(100) ele_name_temp
character(*) attrib_name
character(24) :: r_name = 'pointers_to_attribute'

logical err_flag, do_allocation, do_print
logical, optional :: err_print_flag

! init

err_flag = .false.
do_print = logic_option (.true., err_print_flag)

! beam_start

if (ele_name == 'BEAM_START') then

  if (present(eles)) then
    call re_allocate_eles (eles, 0)
  endif

  call re_allocate (ptr_array, 1)

  select case(attrib_name)
  case ('EMITTANCE_A'); ptr_array(1)%r => lat%a%emit 
  case ('EMITTANCE_B'); ptr_array(1)%r => lat%b%emit
  case ('EMITTANCE_Z'); ptr_array(1)%r => lat%z%emit
  case default
    beam_start%key = def_beam_start$
    ix = attribute_index (beam_start, attrib_name)
    if (ix < 1) then
      if (do_print) call out_io (s_error$, r_name, &
             'INVALID ATTRIBUTE: ' // attrib_name, 'FOR ELEMENT: ' // ele_name)
      deallocate (ptr_array)
      err_flag = .true.
      return
    endif
    if (present(ix_attrib)) ix_attrib = ix
    ptr_array(1)%r => lat%beam_start%vec(ix)
  end select

  return

endif

! Locate elements

call lat_ele_locator (ele_name, lat, eles2, n_loc, err_flag)
if (n_loc == 0) then
  if (do_print) call out_io (s_error$, r_name, 'ELEMENT NOT FOUND: ' // ele_name)
  if (allocated(ptr_array)) deallocate (ptr_array)
  err_flag = .true.
  return  
endif

! Locate attributes

call re_allocate (ptrs, n_loc)
n = 0
do i = 1, n_loc
  call pointer_to_attribute (eles2(i)%ele, &
          attrib_name, do_allocation, ptrs(n+1)%r, err_flag, .false., ix_a)
  if (.not. err_flag) then
    n = n + 1
    if (present(ix_attrib)) ix_attrib = ix_a
  endif
enddo

if (n == 0) then
  if (do_print) call out_io (s_error$, r_name, 'ATTRIBUTE: ' // attrib_name, &
                                               'NOT FOUND FOR: ' // ele_name)
  if (allocated(ptr_array)) deallocate (ptr_array)
  err_flag = .true.
  return  
endif

! Transfer pointers to ptr_array

if (present(eles)) then
  call re_allocate_eles (eles, n)
  eles = eles2(1:n)
endif

call re_allocate (ptr_array, n)
do i = 1, n
  ptr_array(i)%r => ptrs(i)%r
enddo

err_flag = .false.

end subroutine

