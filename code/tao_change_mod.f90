module tao_change_mod

use tao_mod
use quick_plot

contains

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine tao_change_var (name, num_str)
!
! Routine to change a variable in the model lattice.
!
! Input:
!   name     -- Character(*): Name of variable or element.
!   num_str  -- Character(*): Change in value. 
!                                A '@' signifies a absolute set.
!                                A 'd' signifies a set relative design.        
!
! Output:
!    %u(s%global%u_view)%model -- model lattice where the variable lives.
!-

subroutine tao_change_var (name, num_str)

implicit none

type (tao_universe_struct), pointer :: u
type (tao_var_array_struct), allocatable :: v_array(:)
type (tao_var_struct), pointer :: var

real(rp) change_number, model_value, old_merit, new_merit, max_val
real(rp) old_value, new_value, design_value, delta

integer nl, i, ixa, ix, err_num, n

character(*) name, num_str
character(20) :: r_name = 'tao_change_var'
character(20) abs_or_rel
character(100) l1, num, fmt
character(200), allocatable, save :: lines(:)

logical err

!-------------------------------------------------

call re_allocate (lines, 200, 200)
call to_number (num_str, change_number, abs_or_rel, err);  if (err) return
old_merit = tao_merit()
nl = 0

call tao_find_var (err, name, v_array = v_array)
if (err) return
if (.not. allocated(v_array)) then
  call out_io (s_error$, r_name, 'BAD VARIABLE NAME: ' // name)
  return
endif

! now change all desired variables

max_val = 0
do i = 1, size(v_array)
  var => v_array(i)%v
  var%old_value = var%model_value
  if (abs_or_rel == 'ABS') then
    call tao_set_var_model_value (var, change_number)
  elseif (abs_or_rel == 'REL') then
    call tao_set_var_model_value (var, var%design_value + change_number)
  else
    call tao_set_var_model_value (var, var%model_value + change_number)
  endif
  max_val = max(max_val, abs(var%old_value))
  max_val = max(max_val, abs(var%design_value)) 
  max_val = max(max_val, abs(var%model_value)) 
enddo

! print results

new_merit = tao_merit()

if (max_val > 100) then
  fmt = '(5x, I5, 2x, f12.0, a, 4f12.0)'
else
  fmt = '(5x, I5, 2x, f12.6, a, 4f12.6)'
endif

l1 = '     Index         Old              New       Delta  Old-Design  New-Design'
nl=nl+1; lines(nl) = l1

n = size(v_array)
call re_allocate (lines, len(lines(1)), n+100)

do i = 1, size(v_array)
  var => v_array(i)%v
  delta = var%model_value - var%old_value
  nl=nl+1; write (lines(nl), fmt) i, var%old_value, '  ->', &
       var%model_value, delta, var%old_value - var%design_value, &
	     var%model_value - var%design_value
enddo
nl=nl+1; lines(nl) = l1

if (max(abs(old_merit), abs(new_merit)) > 100) then
  fmt = '(5x, 2(a, f13.2), f13.2)'
else
  fmt = '(5x, 2(a, f13.6), f13.6)'
endif

nl=nl+1; lines(nl) = ' '
nl=nl+1; write (lines(nl), fmt) 'Merit:      ', &
                        old_merit, '  ->', new_merit, new_merit-old_merit
if (delta /= 0) then
  nl=nl+1; write (lines(nl), '(5x, a, es11.3)') 'dMerit/dVar:', &
                                     (new_merit-old_merit) / delta
endif

call out_io (s_blank$, r_name, lines(1:nl))

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine tao_change_ele (name, attribute, num_str)
!
! Routine to change a variable in the model lattice.
!
! Input:
!   who       -- Character(*): 'var' or 'ele'
!   name      -- Character(*): Name of variable or element.
!   attribute -- Character(*): Attribute name of element.
!   num_str   -- Character(*): Change in value. 
!                                A '@' signifies a absolute set.
!                                A 'd' signifies a set relative design.        
!
! Output:
!    %u(s%global%u_view)%model -- model lattice where the variable lives.
!-

subroutine tao_change_ele (name, attribute, num_str)

use tao_mod
use quick_plot

implicit none

type (tao_universe_struct), pointer :: u
type (tao_v1_var_struct), pointer :: v1_ptr
type (tao_var_struct), pointer :: var

real(rp) change_number, model_value, old_merit, new_merit, max_val
real(rp) old_value, new_value, design_value, delta

real(rp), pointer :: attrib_ptr, design_ptr

integer nl, i, ixa, ix, err_num, n
integer, allocatable, save :: ix_ele(:)

character(*) name, attribute, num_str
character(80) num
character(20) :: r_name = 'tao_change_ele'
character(16) ele_name, abs_or_rel
character(40) fmt
character(100) l1
character(200), allocatable, save :: lines(:)

logical err

!-------------------------------------------------

call re_allocate (lines, 200, 200)
call to_number (num_str, change_number, abs_or_rel, err);  if (err) return
old_merit = tao_merit()
nl = 0

! 

call string_trim (name, name, ix)

if (ix .gt. 16) then
  call out_io (s_error$, r_name, &
                   'ELEMENT NAME CANNOT BE GREATER THAN 16 CHARACTERS')
  return
endif
  
ele_name = name
  
u => s%u(s%global%u_view)
call tao_locate_element (ele_name, s%global%u_view, ix_ele)
if (ix_ele(1) < 0) return

select case (attribute)
case ('x', 'p_x', 'y', 'p_y', 'z', 'p_z')
  call change_orbit (u%design%orb(ix_ele(1)), u%model%orb(ix_ele(1)))
  if (err) return
case default  
  do i = 1, size(ix_ele)
    call pointer_to_attribute (u%design%lat%ele_(ix_ele(i)), attribute, .true., &
                                                           attrib_ptr, ixa, err)
    if (err) then
      if (size(ix_ele) .eq. 1) then
        return
      else
        cycle
      endif
    endif
    if (.not. attribute_free (ix_ele(i), ixa, u%model%lat, .true.)) return
    design_value = attrib_ptr
     
    call pointer_to_attribute (u%model%lat%ele_(ix_ele(i)), attribute, .true., &
                                                           attrib_ptr, ixa, err)
    old_value = attrib_ptr
     
    if (abs_or_rel == 'ABS') then
      attrib_ptr = change_number
    elseif (abs_or_rel == 'REL') then
      attrib_ptr = design_value + change_number
    else
      attrib_ptr = attrib_ptr + change_number
    endif
     
    new_value = attrib_ptr
  enddo

end select

s%global%lattice_recalc = .true.

! don't print results if changing multiple elements

if (size(ix_ele) /= 1) return

!----------------------------------
! print results

new_merit = tao_merit()
delta = new_value - old_value
if (max(abs(old_value), abs(new_value), abs(design_value)) > 100) then
  fmt = '(5x, 2(a, f12.0), f12.0)'
else
  fmt = '(5x, 2(a, f12.6), f12.6)'
endif

nl=nl+1;write (lines(nl), '(27x, a)') 'Old              New      Delta'
nl=nl+1;write (lines(nl), fmt) 'Value:       ', old_value, '  ->', new_value, delta
nl=nl+1;write (lines(nl), fmt) 'Value-Design:', old_value-design_value, &
                                  '  ->', new_value-design_value

if (max(abs(old_merit), abs(new_merit)) > 100) then
  fmt = '(5x, 2(a, f13.2), f13.2)'
else
  fmt = '(5x, 2(a, f13.6), f13.6)'
endif

nl=nl+1;write (lines(nl), fmt) 'Merit:      ', &
                        old_merit, '  ->', new_merit, new_merit-old_merit
nl=nl+1;lines(nl) = ' '
nl=nl+1;if (delta /= 0) write (lines(nl), '(a, es12.3)') &
                         'dMerit/dValue:  ', (new_merit-old_merit) / delta

call out_io (s_blank$, r_name, lines(1:nl))

!----------------------------------------------------------------------------
! this changes the orbit in a save area which is then injected into model_orb(0)
! at lattice calc time

contains

subroutine change_orbit (design_orb, model_orb)

implicit none

type (coord_struct), target :: design_orb, model_orb 

real(rp), pointer :: design_ptr, model_ptr

integer direction

err = .false.

!check if this is a linear lattice

if ((u%model%lat%param%lattice_type .eq. circular_lattice$) .and. &
    attribute .ne.'p_z') then
  call out_io (s_warn$, r_name, "This is a circular lattice!", &
                    "So only changing the p_z orbit will do anything")
  err = .true.
  return
endif
!check if we are changing the beginning element
if (ix_ele(1) .ne. 0 .or. size(ix_ele) .ne. 1) then
  call out_io (s_warn$, r_name, &
      "Changing the orbit is only applicable for the BEGINNING element!")
  err = .true.
  return
endif

!point to correct direction
select case (attribute)
  case ('x')
    direction = 1
  case ('p_x')
    direction = 2
  case ('y')
    direction = 3
  case ('p_y')
    direction = 4
  case ('z')
    direction = 5
  case ('p_z')
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

if (abs_or_rel == 'ABS') then
  model_ptr = change_number
elseif (abs_or_rel == 'REL') then
  model_ptr = design_value + change_number
else
  model_ptr = model_ptr + change_number
endif

new_value = model_ptr

err = .false.

end subroutine change_orbit

end subroutine tao_change_ele

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

subroutine to_number (num_str, change_number, abs_or_rel, err)

real(rp) change_number

integer ix, ios

character(*) num_str
character(*) abs_or_rel
character(80) n_str
character(20) :: r_name = 'to_number'

logical err

!

n_str = num_str
abs_or_rel = ' ' 

ix = index (n_str, '@')
if (ix /= 0) then
  abs_or_rel = 'ABS'
  n_str(ix:ix) = ' '
endif

ix = index (n_str, 'd')
if (ix /= 0) then
  if (abs_or_rel /= ' ') then
    call out_io (s_error$, r_name, &
        '"@" AND "d" QUALIFIERS CANNONT BOTH BE USED AT THE SAME TIME.')
    return
  endif
  abs_or_rel = 'REL'
  n_str(ix:ix) = ' '
endif

call tao_to_real (n_str, change_number, err)

end subroutine

end module
