module set_ele_attribute_mod

use bmad_parser_struct
use bmad_interface

implicit none

contains

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!+
! Subroutine set_ele_attribute (ele, set_string, err_flag, err_print_flag)
!
! Routine to set an element's attribute.
! This routine is useful for parsing user input.
!
! The set_string argument is in the form:
!   "<attribute_name> = <value>"
! <value> may be any expression that can be evaluated by bmad_parser.
! Examples: 
!   "tracking_method = taylor"
!   "hkick = 0.01 * pi
!
! Note: If using intellegent bookkeeping (bmad_com%auto_bookkeeper = F), the routine
! lattice_bookkeeper must be called after all elements sets are done.
!
! Also see:
!   set_ele_real_attribute
!   
! Input:
!   ele            -- ele_struct: Element with attribute to set.
!   set_string     -- character(*): Attribute and value for set.
!   err_print_flag -- logical, optional: If present and False then suppress printing 
!                       of an error message if attribute is, for example, not free.
!
! Output:
!   ele      -- ele_struct: Element with attribute set.
!   err_flag -- Logical: Set True if there is an error, False otherwise.
!-

subroutine set_ele_attribute (ele, set_string, err_flag, err_print_flag)

type (ele_struct) ele
type (stack_file_struct), target :: current_file
type (bp_common_struct) bp_save

integer ix

logical, optional :: err_print_flag
logical err_flag, delim_found, print_save, file_input_save, is_slaved_field_attribute

character(*) set_string
character(100) string
character(1) delim
character(*), parameter :: r_name = 'set_ele_attribute'
character(40) a_name

! Check if free. Except if we know how to handle the attribute.

err_flag = .true.

call str_upcase (string, set_string)
ix = index(string, '=')
if (ix == 0) then
  call out_io (s_error$, r_name, 'NO "=" SIGN FOUND')
  return
endif

a_name = string(1:ix-1)
if (.not. attribute_free (ele, a_name, err_print_flag, dependent_attribs_free = .true.)) return

! Evaluate and set.
! This essentially is a wrapper for the bmad_parser routine parser_set_attribute.

if (.not. allocated(bp_com2%const)) call init_bmad_parser_common

print_save = bp_com%print_err
file_input_save = bp_com%input_from_file

call init_bmad_parser_common  ! init variable list
bp_com%input_from_file = .false.
bp_com%parser_name = r_name
bp_com%parse_line = string
bp_com%current_file => current_file
bp_com%print_err = logic_option(.true., err_print_flag)
current_file%full_name = ''

call parser_set_attribute (redef$, ele, delim, delim_found, err_flag, set_field_master = .false.)

bp_com%input_from_file = file_input_save
bp_com%print_err       = print_save

if (err_flag) return

call attribute_set_bookkeeping (ele, a_name, err_flag)

end subroutine set_ele_attribute

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!+
! Subroutine set_ele_real_attribute (ele, attrib_name, value, err_flag, err_print_flag)
!
! Routine to set a lattice element's real, or integer named attributes.
! This routine cannot be used with string or logical attributes and should not be used with switch
! attributes (switch attributes are components like ele%mat6_calc_method which are mapped to a
! set of strings ("bmad_standard, etc.).
!
! Note: If using intellegent bookkeeping (bmad_com%auto_bookkeeper = F), the routine
! lattice_bookkeeper must be called after all elements sets are done.
!
! Also see:
!   set_ele_attribute
!   
! Input:
!   ele            -- ele_struct: Element with attribute to set.
!   attrib_name    -- character(*): Attribute name.
!   value          -- real(rp): value to set to.
!   err_print_flag -- logical, optional: If present and False then suppress printing 
!                       of an error message if attribute is, for example, not free.
!
! Output:
!   ele      -- ele_struct: Element with attribute set.
!   err_flag -- Logical: Set True if there is an error, False otherwise.
!-

subroutine set_ele_real_attribute (ele, attrib_name, value, err_flag, err_print_flag)

type (ele_struct) ele
type (all_pointer_struct) a_ptr

real(rp) value

logical, optional :: err_print_flag
logical err_flag

character(*) attrib_name
character(*), parameter :: r_name = 'set_ele_real_attribute'
character(40) a_name

!

err_flag = .true.
call str_upcase (a_name, attrib_name)
if (.not. attribute_free (ele, a_name, err_print_flag, dependent_attribs_free = .true.)) return

call pointer_to_attribute (ele, attrib_name, .true., a_ptr, err_flag)
if (associated(a_ptr%r)) then
  a_ptr%r = value
elseif (associated(a_ptr%i)) then
  a_ptr%i = nint(value)
else
  if (logic_option(.true., err_print_flag)) then
    call out_io (s_error$, r_name, 'BAD ATTRIBUTE: ' // a_name)
  endif
  return
endif

call attribute_set_bookkeeping (ele, a_name, err_flag, a_ptr)

end subroutine set_ele_real_attribute

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
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

end module
