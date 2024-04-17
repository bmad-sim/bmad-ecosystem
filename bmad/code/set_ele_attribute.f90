!+
! Subroutine set_ele_attribute (ele, set_string, err_flag, err_print_flag, set_lords)
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
!   ele             -- ele_struct: Element with attribute to set.
!   set_string      -- character(*): Attribute and value for set.
!   err_print_flag  -- logical, optional: If present and False then suppress printing 
!                        of an error message if attribute is, for example, not free.
!   set_lords       -- logical, optional: Default False. If True, set the super_lord(s) or multipass_lord
!                        if the element is a super_slave or multipass_slave.
!
! Output:
!   ele      -- ele_struct: Element with attribute set.
!   err_flag -- Logical: Set True if there is an error, False otherwise.
!-

subroutine set_ele_attribute (ele, set_string, err_flag, err_print_flag, set_lords)

use bmad_parser_mod, dummy => set_ele_attribute
use bmad_interface, dummy2 => set_ele_attribute


implicit none

type (ele_struct), target :: ele
type (ele_struct), pointer :: lord
type (stack_file_struct), target :: current_file
type (bp_common_struct) bp_save

integer i, ix

logical, optional :: err_print_flag, set_lords
logical err_flag, delim_found, is_slaved_field_attribute

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

a_name = adjustl(string(1:ix-1))

select case (a_name)
case ('SLAVE', 'VAR', 'REF_BEGINNING', 'REF_CENTER', 'REF_END', 'ELE_BEGINNING', &
      'ELE_CENTER', 'ELE_END', 'SUPERIMPOSE', 'REFERENCE', 'OFFSET', 'FIELD_OVERLAPS', &
      'CREATE_JUMBO_SLAVE', 'ELE_ORIGIN', 'REF_ORIGIN', 'WRAP_SUPERIMPOSE', 'TO_ELEMENT')
  if (logic_option(.true., err_print_flag)) call out_io (s_error$, r_name, &
                                'ELEMENT PARAMETER NOT SETTABLE AFTER THE LATTICE HAS BEEN READ IN: ' // a_name)
  return
end select

! Evaluate and set.
! This essentially is a wrapper for the bmad_parser routine parser_set_attribute.

if (.not. allocated(bp_com2%const)) call init_bmad_parser_common

bp_save = bp_com

call init_bmad_parser_common  ! init variable list

bp_com%input_from_file = .false.
bp_com%parser_name = r_name
bp_com%parse_line = string
bp_com%current_file => current_file
bp_com%print_err = logic_option(.true., err_print_flag)
bp_com%undefined_vars_evaluate_to_zero = .false.

current_file%full_name = ''

if (ele%slave_status == super_slave$ .or. ele%slave_status == multipass_slave$) then
  do i = 1, ele%n_lord
    lord => pointer_to_lord(ele, i)
    if (lord%slave_status == multipass_slave$) lord => pointer_to_lord(ele, 1)
    if (lord%lord_status /= super_lord$ .and. lord%lord_status /= multipass_lord$) cycle
    if (.not. attribute_free (lord, a_name, err_print_flag, dependent_attribs_free = .true.)) return
    call parser_set_attribute (redef$, lord, delim, delim_found, err_flag, set_field_master = .false.)
    if (err_flag) exit
    call attribute_set_bookkeeping (lord, a_name, err_flag)
  enddo
    
else
  if (.not. attribute_free (ele, a_name, err_print_flag, dependent_attribs_free = .true.)) return
  call parser_set_attribute (redef$, ele, delim, delim_found, err_flag, set_field_master = .false.)
  if (.not. err_flag) call attribute_set_bookkeeping (ele, a_name, err_flag)
endif

bp_com = bp_save

end subroutine set_ele_attribute

