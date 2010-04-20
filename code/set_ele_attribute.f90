!+
! Subroutine set_ele_attribute (ele, set_string, lat, err_flag)
!
! Routine to set an element's attribute.
! This is a general routine that will set almost any element attribute.
!
! The set_string argument is in the form:
!   <attribute_name> = <value>
! <value> may be any expression that can be evaluated by bmad_parser.
! Examples: 
!   "tracking_method = taylor"
!   "hkick = 0.01 * pi
!
! Modules needed:
!   use bmad
!
! Input:
!   ele        -- ele_struct: Element with attribute to set.
!   set_string -- Character(*): attribute and value for set.
!   lat        -- lat_struct: Lattice containing the element. 
!                   Used if an evaluation is needed
!
! Output:
!   ele      -- ele_struct: Element with attribute set.
!   err_flag -- Logical: Set True if there is an error, False otherwise.
!-

subroutine set_ele_attribute (ele, set_string, lat, err_flag)

use bmad_parser_mod

implicit none

type (ele_struct) ele
type (lat_struct) lat
type (stack_file_struct), target :: current_file

integer ix

character(*) set_string
character(100) string
character(1) delim
character(20) :: r_name = 'set_ele_attribute'

logical err_flag, delim_found

! Check if free

err_flag = .true.

call str_upcase (string, set_string)
ix = index(string, '=')
if (ix == 0) then
  call out_io (s_error$, r_name, 'NO "=" SIGN FOUND')
  return
endif

if (.not. attribute_free (ele, string(1:ix-1), lat)) return

! Evaluate and set.
! This essentially is a wrapper for the bmad_parser routine get_attribute.

if (.not. allocated(bp_com%var_name)) call init_bmad_parser_common
bp_com%input_from_file = .false.
bp_com%parser_name = r_name
bp_com%parse_line = string
bp_com%current_file => current_file
current_file%full_name = ''

call get_attribute (redef$, ele, lat, delim, delim_found, err_flag, .true.)

end subroutine
