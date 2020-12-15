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

use bmad_parser_mod, dummy => set_ele_attribute
use bmad_interface, dummy2 => set_ele_attribute


implicit none

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

