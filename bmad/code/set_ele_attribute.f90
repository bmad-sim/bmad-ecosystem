!+
! Subroutine set_ele_attribute (ele, set_string, err_flag, err_print_flag, set_lords, err_id)
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
! Note: The routine lattice_bookkeeper must be called after all elements sets are done.
!
! Also see:
!   set_ele_real_attribute
!   
! Input:
!   ele             -- ele_struct: Element with attribute to set.
!   set_string      -- character(*): Attribute and value for set.
!   err_print_flag  -- logical, optional: If present and False then suppress printing 
!                        of an error message if attribute is, for example, not free.
!   set_lords       -- logical, optional: Default False. If True, set the super_lord(s)
!                        if the element is a super_slave.
!
! Output:
!   ele      -- ele_struct: Element with attribute set.
!   err_flag -- Logical: Set True if there is an error, False otherwise.
!   err_id   -- integer, optional: Set to an integer which identifies the error type. 0 = no error.
!                 The higher the error the further along the error was encountered.
!-

subroutine set_ele_attribute (ele, set_string, err_flag, err_print_flag, set_lords, err_id)

use bmad_interface, dummy => set_ele_attribute
use parser_set_attribute_mod, dummy2 => set_ele_attribute

implicit none

type (ele_struct), target :: ele
type (ele_struct), pointer :: lord
type (stack_file_struct), target :: current_file
type (bp_common_struct) bp_save

integer, optional :: err_id
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
if (present(err_id)) err_id = 0

ix = index(set_string, '=')
if (ix == 0) then
  if (logic_option(.true., err_print_flag)) call out_io (s_error$, r_name, 'NO "=" SIGN FOUND')
  if (present(err_id)) err_id = 1
  return
endif

a_name = upcase(adjustl(set_string(1:ix-1)))
if (attribute_type(a_name) == is_string$) then
  string = trim(a_name) // set_string(ix:)
else
  call str_upcase (string, set_string)
endif

select case (a_name)
case ('SLAVE', 'VAR', 'REF_BEGINNING', 'REF_CENTER', 'REF_END', 'ELE_BEGINNING', &
      'ELE_CENTER', 'ELE_END', 'SUPERIMPOSE', 'REFERENCE', 'OFFSET', 'FIELD_OVERLAPS', &
      'CREATE_JUMBO_SLAVE', 'ELE_ORIGIN', 'REF_ORIGIN', 'WRAP_SUPERIMPOSE', 'TO_ELEMENT')
  if (logic_option(.true., err_print_flag)) call out_io (s_error$, r_name, &
                                'ELEMENT PARAMETER NOT SETTABLE AFTER THE LATTICE HAS BEEN READ IN: ' // a_name)
  if (present(err_id)) err_id = 2
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

if (logic_option(.false., set_lords) .and. ele%slave_status == super_slave$) then
  do i = 1, ele%n_lord
    lord => pointer_to_lord(ele, i)
    if (lord%slave_status == multipass_slave$) lord => pointer_to_lord(ele, 1)
    if (lord%lord_status /= super_lord$ .and. lord%lord_status /= multipass_lord$) cycle
    if (.not. attribute_free (lord, a_name, err_print_flag, dependent_attribs_free = .true.)) then
      if (present(err_id)) err_id = 3
      goto 8000
      return
    endif
    call parser_set_attribute (redef$, lord, delim, delim_found, err_flag, set_field_master = .false.)
    if (err_flag) then
      if (present(err_id)) err_id = 4
      goto 8000
    endif
    call attribute_set_bookkeeping (lord, a_name, err_flag)
    if (err_flag) then
      if (present(err_id)) err_id = 5
      goto 8000
    endif
  enddo
    
else
  if (.not. attribute_free (ele, a_name, err_print_flag, dependent_attribs_free = .true.)) then
    if (present(err_id)) err_id = 3
    goto 8000
  endif

  call parser_set_attribute (redef$, ele, delim, delim_found, err_flag, set_field_master = .false.)
  if (err_flag) then
    if (present(err_id)) err_id = 4
    goto 8000
  endif

  call attribute_set_bookkeeping (ele, a_name, err_flag)
endif

8000 continue
bp_com = bp_save

end subroutine set_ele_attribute

