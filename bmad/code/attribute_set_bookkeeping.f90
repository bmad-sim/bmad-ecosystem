!+
! Subroutine attribute_set_bookkeeping (ele, attrib_name, err_flag, attrib_ptr)
!
! Bookkeeping routine to be used after an element attribute is set.
!
! Input:
!   ele         -- ele_struct: Element containing the attribute
!   attrib_name -- character(*): Name of the attribute. Must be upper case.
!   attrib_ptr  -- all_pointer_struct, optional: Pointer to the attribute.
!                     The presence of this argument saves a small amount of time.
!-

subroutine attribute_set_bookkeeping (ele, attrib_name, err_flag, attrib_ptr)

use bmad_interface, dummy => attribute_set_bookkeeping

implicit none

type (ele_struct) ele
type (all_pointer_struct), optional :: attrib_ptr
type (all_pointer_struct) a_ptr
type (branch_struct), pointer :: branch

character(*) attrib_name
logical err_flag

!

select case (attrib_name)
case ('VOLTAGE')
  if (ele%value(l$) /= 0) then
    ele%value(gradient$) = ele%value(voltage$) / ele%value(l$)
  endif

case ('GRADIENT')
    ele%value(voltage$) = ele%value(gradient$) * ele%value(l$)

case default
  if (.not. field_attribute_free(ele, attrib_name)) then
    ele%field_master = .not. ele%field_master
    branch => pointer_to_branch(ele)
    call attribute_bookkeeper(ele, force_bookkeeping = .true.)
    ele%field_master = .not. ele%field_master
  endif
end select

! Set bookkeeping flags

if (present(attrib_ptr)) then
  call pointer_to_attribute (ele, attrib_name, .true., attrib_ptr, err_flag)
  call set_flags_for_changed_attribute (ele, attrib_ptr)
else
  if (attribute_type(attrib_name) == is_string$) return  ! No bookkeeping needed
  call pointer_to_attribute (ele, attrib_name, .true., a_ptr, err_flag)
  call set_flags_for_changed_attribute (ele, a_ptr)
endif

end subroutine attribute_set_bookkeeping
