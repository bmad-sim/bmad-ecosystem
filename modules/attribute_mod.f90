module attribute_mod

use bmad_struct
use basic_bmad_interface
use multipole_mod
use lat_ele_loc_mod

type(ele_struct), private, pointer, save :: ele0 ! For Error message purposes
character(40), private, save :: attrib_name0     ! For Error message purposes

private check_this_attribute_free, print_error

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function attribute_free
!
! Overloaded function for:
!   Function attribute_free1 (ix_ele, attrib_name, lat,
!                                err_print_flag, except_overlay) result (free)
!   Function attribute_free2 (ele, attrib_name, lat, 
!                                err_print_flag, except_overlay) result (free)
!   Function attribute_free3 (ix_ele, ix_branch, attrib_name, lat, 
!                                err_print_flag, except_overlay) result (free)
!
! Routine to check if an attribute is free to vary.
!
! Attributes that cannot be changed directly include super_slave attributes (since
! these attributes are controlled by their super_lords) and attributes that
! are controlled by an overlay_lord.
!
! Also dependent variables such as the angle of a bend cannot be 
!   freely variable.
!
! Modules needed:
!   use bmad
!
! Input:
!   ix_ele          -- Integer: Index of element in element array.
!   ix_branch       -- Integer: Branch index of element. 
!   ele             -- Ele_struct: Element containing the attribute
!   attrib_name     -- Character(*): Name of the attribute. Assumed upper case.
!   lat             -- lat_struct: Lattice structure.
!   err_print_flag  -- Logical, optional: If present and False then suppress
!                       printing of an error message if attribute is not free.
!   except_overlay  -- Logical, optional: If present and True then an attribute that
!                       is an overlay_slave will be treated as free. This is used by,
!                       for example, the create_overlay routine.
!
! Output:
!   free   -- Logical: Set True if attribtute not found or attriubte
!                     cannot be changed directly.
!-

interface attribute_free
  module procedure attribute_free1
  module procedure attribute_free2
  module procedure attribute_free3
end interface

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function ele_attribute_value (ele, attrib_name, do_allocation,
!                            attrib_value, err_flag, err_print_flag, ix_attrib)
!
! Returns the value of an element attribute. 
! This routine has the advantage of being able to handle logical and integer parameters.
! Logical values are returned as 0 => False, and 1 => True  
!
! Note: Use attribute_free to see if the attribute may be varied independently.
! Note: Alternatively, consider the routines:
!     pointer_to_attribute
!     pointers_to_attribute
!
! Modules needed:
!   use bmad
!
! Input:
!   ele             -- Ele_struct: After this routine finishes Ptr_attrib 
!                        will point to a variable within this element.
!   attrib_name     -- Character(40): Name of attribute. Must be uppercase.
!                       For example: "HKICK".
!   do_allocation   -- Logical: If True then do an allocation if needed.
!                       EG: The multipole An and Bn arrays need to be allocated
!                       before their use.
!   err_print_flag  -- Logical, optional: If present and False then suppress
!                       printing of an error message on error.
!
! Output:
!   attrib_value -- Real(rp): Value of the attribute.
!   err_flag     -- Logical: Set True if attribtute not found. False otherwise.
!   ix_attrib    -- Integer, optional: If applicable then this is the index to the 
!                     attribute in the ele%value(:), ele%a_pole(:) or ele%b_pole arrays.
!-

subroutine ele_attribute_value (ele, attrib_name, do_allocation, &
                  attrib_value, err_flag, err_print_flag, ix_attrib)


type (ele_struct), target :: ele

real(rp) attrib_value
real(rp), pointer :: ptr_attrib

integer, optional :: ix_attrib

character(*) attrib_name
character(24) :: r_name = 'ele_attribute_value'
character(40) a_name

logical err_flag, do_allocation
logical, optional :: err_print_flag

! Init

attrib_value = 0
err_flag = .false.
call str_upcase (a_name, attrib_name)

! Special cases


select case (a_name)
case ('NUM_STEPS')
  attrib_value = ele%num_steps
  if (present(ix_attrib)) ix_attrib = num_steps$
  return
case ('IS_ON')
  attrib_value = logic_val(ele%is_on)
  if (present(ix_attrib)) ix_attrib = is_on$
  return
end select

!

call pointer_to_attribute (ele, attrib_name, do_allocation, &
                  ptr_attrib, err_flag, err_print_flag, ix_attrib)
if (associated(ptr_attrib)) attrib_value = ptr_attrib


!---------------------------------------------
contains

function logic_val (l_val) result (real_val)

logical l_val
real(rp) real_val

!

if (l_val) then
  real_val = 1
else
  real_val = 0
endif

end function

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine pointer_to_indexed_attribute (ele, ix_attrib, do_allocation,
!                                      ptr_attrib, err_flag, err_print_flag)
!
! Returns a pointer to an attribute of an element ele with attribute index ix_attrib.
! 
! Use of this routine is restricted to attributes that have an index. That is,
! attributes in the ele%value(:) array and ele%a_pole(:), and ele%b_pole(:) values.
! A more general routine is pointer_to_attribute.
! Alternatively, consider the routine pointers_to_attribute.
! Note: Use attribute_free to see if the attribute may be varied independently.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele             -- Ele_struct: After this routine finishes Ptr_attrib 
!                        will point to a variable within this element.
!   ix_attrib       -- Integer, Attribute index.
!   do_allocation   -- Logical: If True then do an allocation if needed.
!                       EG: The multipole An and Bn arrays need to be allocated
!                       before their use.
!   err_print_flag  -- Logical, optional: If present and False then suppress
!                       printing of an error message on error.
!
! Output:
!   ptr_attrib -- Real(rp), pointer: Pointer to the attribute.
!                     Pointer will be deassociated if there is a problem.
!   err_flag   -- Logical: Set True if attribtute not found. False otherwise.
!-

subroutine pointer_to_indexed_attribute (ele, ix_attrib, do_allocation, &
                                        ptr_attrib, err_flag, err_print_flag)

implicit none

type (ele_struct), target :: ele

real(rp), pointer :: ptr_attrib

integer :: ix_attrib

character(40) :: r_name = 'pointer_to_indexed_attribute'

logical err_flag, do_allocation, do_print
logical, optional :: err_print_flag

! Init

err_flag = .true.
nullify (ptr_attrib)
do_print = logic_option (.true., err_print_flag)

! multipole?

if (ix_attrib >= a0$ .and. ix_attrib <= b20$) then   ! multipole attribute

  if (.not. associated(ele%a_pole)) then
    if (do_allocation) then
      call multipole_init (ele)
    else
      if (do_print) call out_io (s_error$, r_name, &
                      'MULTIPOLE NOT ALLOCATED FOR ELEMENT: ' // ele%name)
      return
    endif
  endif

  if (ix_attrib >= b0$) then
    ptr_attrib => ele%b_pole(ix_attrib-b0$)
  else
    ptr_attrib => ele%a_pole(ix_attrib-a0$)
  endif

elseif (ix_attrib < 1 .or. ix_attrib > n_attrib_maxx) then
  if (do_print) call out_io (s_error$, r_name, &
          'INVALID ATTRIBUTE INDEX: \i0\ ', 'FOR THIS ELEMENT: ' // ele%name, &
          i_array = (/ ix_attrib /))
  return

! otherwise must be in ele%value(:) array

else
  ptr_attrib => ele%value(ix_attrib)
endif


err_flag = .false.
return

end subroutine pointer_to_indexed_attribute 

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function attribute_free1 (ix_ele, attrib_name, 
!                                 lat, err_print_flag, except_overlay) result (free)
!
! This function overloaded by attribute_free. See attribute_free for more details.
!-

function attribute_free1 (ix_ele, attrib_name, lat, &
                                  err_print_flag, except_overlay) result (free)

implicit none

type (lat_struct) :: lat

integer ix_ele

character(*) attrib_name

logical free, do_print, do_except_overlay
logical, optional :: err_print_flag, except_overlay

!

do_print = logic_option (.true., err_print_flag)
do_except_overlay = logic_option(.false., except_overlay)

call check_this_attribute_free (lat%ele(ix_ele), attrib_name, lat, &
                                             do_print, do_except_overlay, free, 0)

end function attribute_free1

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function attribute_free2 (ele, attrib_name, 
!                                 lat, err_print_flag, except_overlay) result (free)
!
! This function overloaded by attribute_free. See attribute_free for more details.
!-

function attribute_free2 (ele, attrib_name, lat, err_print_flag, except_overlay) result (free)

implicit none

type (lat_struct), target :: lat
type (ele_struct) ele

character(*) attrib_name

logical free, do_print, do_except_overlay
logical, optional :: err_print_flag, except_overlay

character(16) :: r_name = 'attribute_free'

! init & check

do_print = logic_option (.true., err_print_flag)
do_except_overlay = logic_option(.false., except_overlay)

call check_this_attribute_free (ele, attrib_name, lat, do_print, do_except_overlay, free, 0)

end function attribute_free2

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function attribute_free3 (ix_ele, ix_branch, 
!                           lat, err_print_flag, except_overlay) result (free)
!
! This function overloaded by attribute_free. See attribute_free for more details.
!-

function attribute_free3 (ix_ele, ix_branch, attrib_name, &
                            lat, err_print_flag, except_overlay) result (free)

implicit none

type (lat_struct) :: lat

integer ix_ele, ix_branch
character(*) attrib_name

logical free, do_print, do_except_overlay
logical, optional :: err_print_flag, except_overlay

!

do_print = logic_option (.true., err_print_flag)
do_except_overlay = logic_option(.false., except_overlay)

call check_this_attribute_free (lat%branch(ix_branch)%ele(ix_ele), attrib_name, &
                                             lat, do_print, do_except_overlay, free, 0)

end function attribute_free3

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

recursive subroutine check_this_attribute_free (ele, attrib_name, lat, &
                            do_print, do_except_overlay, free, ix_recursion, ix_lord)

type (ele_struct), target :: ele
type (lat_struct), target :: lat
type (ele_struct), pointer :: ele_p, lord

integer ix_branch, ix_recursion, i, ir, ix_attrib, ix, ic
integer, optional :: ix_lord

character(*) attrib_name
character(40) a_name

logical free, do_print, do_except_overlay

! If this is first time then set pointers used in error message printing.

if (ix_recursion == 0) then
  ele0 => lat%branch(ele%ix_branch)%ele(ele%ix_ele)
  attrib_name0 = attrib_name
endif

! Check that the name corresponds to an attribute

free = .false.

ix_attrib = attribute_index(ele, attrib_name)
if (ix_attrib < 1) then
  if (do_print) call print_error (ele, ix_attrib, &
          'THIS ATTRIBUTE INDEX DOES NOT CORRESPOND TO A VALID ATTRIBUTE.')
  return
endif

a_name = attribute_name (ele, ix_attrib)

! csr_calc_on, etc. are always free

select case (a_name)
case ('CSR_CALC_ON', 'IS_ON')
  free = .true.
  return
end select

! if the attribute is controled by an overlay lord then it cannot be varied.
! Exception: Multiple overlays can control the same attribute.

if (.not. do_except_overlay) then
  do i = 1, ele%n_lord
    lord => pointer_to_lord (lat, ele, i, ix)
    if (present(ix_lord)) then
      if (ix_lord == lord%ix_ele) cycle
      if (lat%ele(ix_lord)%lord_status == overlay_lord$) cycle
    endif
    if (lord%lord_status == overlay_lord$) then
      if (lat%control(ix)%ix_attrib == ix_attrib) then 
        if (do_print) call print_error (ele, ix_attrib, & 
           'IT IS CONTROLLED BY THE OVERLAY_LORD: ' // lord%name)
        return
      endif
    endif
  enddo
endif

! Super_slaves attributes cannot be varied

if (ele%slave_status == super_slave$) then
  if (do_print) call print_error (ele, ix_attrib, 'THIS ELEMENT IS A SUPER_SLAVE.')
  return
endif

! Check for a multipass_slave.
! Exception: dphi0 can be varied for lcavity and rfcavity slaves, etc.

if (ele%slave_status == multipass_slave$) then
  free = .true.
  select case (ele%key)
  case (lcavity$, rfcavity$) 
    if (ix_attrib == dphi0$) return
  case (patch$)
    lord => pointer_to_lord (lat, ele, 1)
    if (associated (pointer_to_slave(lat, lord, 1), ele) .and. &
                                    ele%ref_orbit == patch_in$) return
  end select

  free = .false.
  if (do_print) call print_error (ele, ix_attrib, 'THIS ELEMENT IS A MULTIPASS_SLAVE.')
  return
endif

! only one particular attribute of an overlay lord is allowed to be adjusted

if (ele%lord_status == overlay_lord$) then
  if (ix_attrib /= ele%ix_value) then
    if (do_print) call print_error (ele, ix_attrib, &
           'FOR THIS OVERLAY ELEMENT THE ATTRIBUTE TO VARY IS: ' // ele%component_name)
    return
  endif
endif

if (ele%lord_status == group_lord$) then
  if (ix_attrib /= command$ .and. ix_attrib /= old_command$) then
    if (do_print) call print_error (ele, ix_attrib, &
          'FOR THIS GROUP ELEMENT THE ATTRIBUTE TO VARY IS: "COMMAND" OR "OLD_COMMAND"')
    return
  endif
endif

! check if it is a dependent variable.

free = .true.

select case (ele%key)
case (sbend$)
  if (any(ix_attrib == (/ angle$, l_chord$, rho$ /))) free = .false.
case (rfcavity$)
  if (ix_attrib == rf_frequency$ .and. ele%value(harmon$) /= 0) free = .false.
case (beambeam$)
  if (ix_attrib == bbi_const$) free = .false.
case (wiggler$)
  if (ix_attrib == k1$ .or. ix_attrib == rho$) free = .false. 
case (lcavity$)
  if (any(ix_attrib == (/ delta_e$, p0c_start$, e_tot_start$ /))) free = .false.
case (elseparator$)
  if (ix_attrib == e_field$ .or. ix_attrib == voltage$) free = .false.
end select

if (has_orientation_attributes(ele%key)) then
  if (ix_attrib == tilt_tot$) free = .false.
  if (ix_attrib == x_pitch_tot$) free = .false.
  if (ix_attrib == y_pitch_tot$) free = .false.
  if (ix_attrib == x_offset_tot$) free = .false.
  if (ix_attrib == y_offset_tot$) free = .false.
  if (ix_attrib == s_offset_tot$) free = .false.
endif

if (ix_attrib == e_tot$) free = .false.
if (ix_attrib == p0c$) free = .false.

if (ele%key == sbend$ .and. ele%lord_status == multipass_lord$ .and. &
    ele%value(n_ref_pass$) == 0 .and. ix_attrib == p0c$) free = .true.

if (.not.free) then
  if (do_print) call print_error (ele, ix_attrib, 'THE ATTRIBUTE IS A DEPENDENT VARIABLE.')
  return
endif

! field_master on means that the b_field and bn_gradient values control
! the strength.

if (ele%field_master) then
  select case (ele%key)
  case (quadrupole$)
    if (ix_attrib == k1$) free = .false.
  case (sextupole$)
    if (ix_attrib == k2$) free = .false.
  case (octupole$)
    if (ix_attrib == k3$) free = .false.
  case (solenoid$)
    if (ix_attrib == ks$) free = .false.
  case (sol_quad$)
    if (ix_attrib == ks$) free = .false.
    if (ix_attrib == k1$) free = .false.
    if (ix_attrib == k2$) free = .false.
  case (sbend$)
    if (ix_attrib == g$) free = .false.
    if (ix_attrib == g_err$) free = .false.
  case (hkicker$, vkicker$)
    if (ix_attrib == kick$) free = .false.
  end select

  if (has_hkick_attributes(ele%key)) then
    if (ix_attrib == hkick$) free = .false.
    if (ix_attrib == vkick$) free = .false.
  endif

else
  select case (ele%key)
  case (quadrupole$)
    if (ix_attrib == b1_gradient$) free = .false.
  case (sextupole$)
    if (ix_attrib == b2_gradient$) free = .false.
  case (octupole$)
    if (ix_attrib == b3_gradient$) free = .false.
  case (solenoid$)
    if (ix_attrib == bs_field$) free = .false.
  case (sol_quad$)
    if (ix_attrib == bs_field$) free = .false.
    if (ix_attrib == b1_gradient$) free = .false.
  case (sbend$)
    if (ix_attrib == b_field$) free = .false.
    if (ix_attrib == b_field_err$) free = .false.
    if (ix_attrib == b1_gradient$) free = .false.
    if (ix_attrib == b2_gradient$) free = .false.
  case (hkicker$, vkicker$)
    if (ix_attrib == bl_kick$) free = .false.
  end select

  if (has_hkick_attributes(ele%key)) then
    if (ix_attrib == bl_hkick$) free = .false.
    if (ix_attrib == bl_vkick$) free = .false.
  endif

endif

if (.not. free) then
  if (do_print) call print_error (ele, ix_attrib, &
       "THE ATTRIBUTE IS A DEPENDENT VARIABLE SINCE", &
       "THE ELEMENT'S FIELD_MASTER IS " // on_off_logic (ele%field_master))
  return
endif

! check slaves

if (ele%lord_status == group_lord$ .or. ele%lord_status == overlay_lord$) then
  do i = 1, ele%n_slave
    ele_p => pointer_to_slave(lat, ele, i, ic)
    call check_this_attribute_free (ele_p, attribute_name(ele_p, lat%control(ic)%ix_attrib), &
                            lat, do_print, do_except_overlay, free, 1, ele%ix_ele)
    if (.not. free) return
  enddo
endif

end subroutine

!-------------------------------------------------------

subroutine print_error (ele, ix_attrib, l1, l2)

type (ele_struct) ele

integer ix_attrib, nl

character(*) l1
character(*), optional :: l2
character(100) li(8)
character (20) :: r_name = 'attribute_free'

!

nl = 0

nl=nl+1; li(nl) =   'THE ATTRIBUTE: ' // attrib_name0
nl=nl+1; li(nl) =   'OF THE ELEMENT: ' // ele0%name

if (ele%ix_branch == 0 .and. ele%ix_ele == ele0%ix_ele) then
  nl=nl+1; li(nl) = 'IS NOT FREE TO VARY SINCE:'
else 
  nl=nl+1; li(nl) = 'IS NOT FREE TO VARY SINCE IT IS TRYING TO CONTROL:'
  nl=nl+1; li(nl) = 'THE ATTRIBUTE: ' // attrib_name0
  nl=nl+1; li(nl) = 'OF THE ELEMENT: ' // ele%name
  nl=nl+1; li(nl) = 'AND THIS IS NOT FREE TO VARY SINCE'
endif

nl=nl+1; li(nl) = l1
if (present(l2)) then
  nl=nl+1; li(nl) = l2
endif

call out_io (s_error$, r_name, li(1:nl))   

end subroutine

end module
