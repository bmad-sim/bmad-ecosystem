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

use attribute_mod, dummy => set_ele_real_attribute

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

call pointer_to_attribute (ele, attrib_name, .true., a_ptr, err_flag, do_unlink = .true.)
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
