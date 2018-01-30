!+
! Function value_of_attribute (ele, attrib_name, err_flag, err_print_flag, err_value) result (value)
!
! Returns the value of an element attribute.
!
! For logical attributes, to convert the value returned by this routine into a logical, use the construction:
!     is_true(value_of_attribute(...))
! Similarly the nint() intrinsic fuction can be used to convert the value returned into an integer if needed.
!
! Note: If you want to set the value of the attribute, use alternatively
!     pointer_to_attribute
!     pointers_to_attribute
!     set_ele_attribute
!
! Input:
!   ele             -- Ele_struct: After this routine finishes Ptr_attrib 
!                        will point to a variable within this element.
!   attrib_name     -- Character(40): Name of attribute. Must be uppercase.
!                       For example: "HKICK".
!   err_print_flag  -- Logical, optional: If present and True then print an error message if there is an  error.
!   err_value       -- real(rp), optional: Value to set value argument if there is an error. Default is 0.
!
! Output:
!   value      -- real(rp): Value of the attribute. Set to err_value if not found.
!   err_flag   -- Logical: Set True if attribtute not found. False otherwise.
!-

function value_of_attribute (ele, attrib_name, err_flag, err_print_flag, err_value) result (value)

use bmad_interface, except_dummy => value_of_attribute

implicit none

type (ele_struct), target :: ele
type (all_pointer_struct) a_ptr

real(rp) value
real(rp), optional :: err_value

character(*) attrib_name
character(24) :: r_name = 'value_of_attribute'

logical err
logical, optional :: err_print_flag, err_flag

!

call pointer_to_attribute (ele, attrib_name, .true., a_ptr, err, logic_option(.false., err_print_flag))
if (present(err_flag)) err_flag = err

if (err) then
  value = real_option(0.0_rp, err_value)
else
  value = value_of_all_ptr(a_ptr)
endif

end function
