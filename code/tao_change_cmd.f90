!+
! Subroutine tao_change_cmd (who, name, where, num_str)
!
! Routine to change a variable in the model lattice.
!
! Input:
!   who      -- Character(*): 'var', 'ele' or 'begin'
!   name     -- Character(*): Name of variable or element.
!   where    -- Character(*): For variable: Index of variable.
!                             For element: Attribute name of element.
!                             For begin: axis to change
!   num_str  -- Character(*): Change in value. 
!                                A '@' signifies a absolute set.
!                                A 'd' signifies a set relative design.        
!
! Output:
!    %u(s%global%u_view)%model -- model lattice where the variable lives.
!-

subroutine tao_change_cmd (who, name, where, num_str)

use tao_mod
use quick_plot

implicit none

type (tao_universe_struct), pointer :: u
type (tao_var_struct), pointer :: v_ptr

real(rp) change_number, model_value
real(rp), pointer :: attrib_ptr, design_ptr
real(rp) old_value, new_value, design_value

integer i, ix_ele, ixa, ix

character(*) who, name, where, num_str
character(80) num
character(20) :: r_name = 'tao_change_cmd'
character(16) ele_name
character(40) fmt
character(80) line(2)

logical err, absolute_num, rel_to_design

!-------------------------------------------------

call to_number;  if (err) return

! If changing a variable...

select case (who)

case ('var')

  call tao_find_var (err, name, var_number = where, v_ptr = v_ptr)
  if (err) return
  if (.not. v_ptr%exists) then
    call out_io (s_error$, r_name, 'VARIABLE DOES NOT EXIST.')
    return
  endif

  old_value = v_ptr%model_value

  if (absolute_num) then
    call tao_set_var_model_value (v_ptr, change_number)
  elseif (rel_to_design) then
    call tao_set_var_model_value (v_ptr, change_number + v_ptr%design_value)
  else
    call tao_set_var_model_value (v_ptr, v_ptr%model_value + change_number)
  endif

  new_value = v_ptr%model_value
  design_value = v_ptr%design_value

!-------------------------------------------------
! If changing an element attribute...

case ('ele')

  call string_trim (name, name, ix)

  if (ix .gt. 16) then
    call out_io (s_error$, r_name, 'ELEMENT NAME CANNOT BE &
  		GREATER THAN 16 CHARACTERS')
    return
  endif
  
  ele_name = name
  
  u => s%u(s%global%u_view)
  call tao_locate_element (ele_name, u%model, ix_ele)
  if (ix_ele < 0) return

  select case (where)
  case ('x', 'x_p', 'y', 'y_p', 'z', 'z_p')
    call change_orbit (u%design_orb(ix_ele), u%model_orb(ix_ele))
  case default  
    call pointer_to_attribute (u%design%ele_(ix_ele), where, .true., &
                                                         attrib_ptr, ixa, err)
    design_value = attrib_ptr
 
    call pointer_to_attribute (u%model%ele_(ix_ele), where, .true., &
                                                         attrib_ptr, ixa, err)
    if (err) return
 
    old_value = attrib_ptr
 
    if (absolute_num) then
      attrib_ptr = change_number
    elseif (rel_to_design) then
      attrib_ptr = design_value + change_number
    else
      attrib_ptr = attrib_ptr + change_number
    endif
 
    new_value = attrib_ptr
  end select

!

case default

  call out_io (s_error$, r_name, 'UNKNOWN "WHO" (ele or var?): ' // trim(who))

end select

!---------------------------------------------------
! print results

if (max(abs(old_value), abs(new_value), abs(design_value)) > 100) then
  fmt = '(5x, 2(a, f11.0), 6x, a, f11.0)'
else
  fmt = '(5x, 2(a, f11.6), 6x, a, f11.6)'
endif

write (line(1), fmt) 'Change in Value:   ', old_value, ' ->', new_value, &
                                                 'Del =', new_value-old_value
write (line(2), fmt) 'Relative to Design:', old_value-design_value, &
                                  ' ->', new_value-design_value

call out_io (s_blank$, r_name, line)

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
contains

subroutine to_number

  integer ios

!

  num = num_str
  err = .true.

  ix = index (num, '@')
  absolute_num = .false.
  if (ix /= 0) then
    absolute_num = .true.
    num(ix:ix) = ' '
  endif

  ix = index (num, 'd')
  rel_to_design = .false.
  if (ix /= 0) then
    rel_to_design = .true.
    num(ix:ix) = ' '
  endif

  if (absolute_num .and. rel_to_design) then
    call out_io (s_error$, r_name, &
        '"@" AND "d" QUALIFIERS CANNONT BOTH BE USED AT THE SAME TIME.')
    return
  endif

  read (num, *, iostat = ios) change_number
  if (ios /= 0) then
    call out_io (s_error$, r_name, 'BAD NUMBER: ' // num)
    return
  endif

  err = .false.

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! contains

subroutine change_orbit (design_orb, model_orb)

implicit none

type (coord_struct), target :: design_orb, model_orb 

real(rp), pointer :: design_ptr, model_ptr

integer direction

!check if this is a linear lattice

!check if we are changin the beginning element

!point to correct direction
select case (where)
  case ('x')
    direction = 1
  case ('x_p')
    direction = 2
  case ('y')
    direction = 3
  case ('y_p')
    direction = 4
  case ('z')
    direction = 5
  case ('z_p')
    direction = 6
  case default
    err = .true.
    return
end select

design_ptr => design_orb%vec(direction)
model_ptr  => model_orb%vec(direction)

!do the change
design_value = design_ptr
old_value  = model_ptr

if (absolute_num) then
  model_ptr = change_number
elseif (rel_to_design) then
  model_ptr = design_value + change_number
else
  model_ptr = model_ptr + change_number
endif

new_value = model_ptr

err = .false.

end subroutine change_orbit

end subroutine tao_change_cmd




