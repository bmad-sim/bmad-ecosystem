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

logical err, exists

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

! We need at least one variable to exist.

exists = .false.
do i = 1, size(v_array)
  if (v_array(i)%v%exists) exists = .true.
enddo

if (.not. exists) then
  call out_io (s_error$, r_name, 'VARIABLE DOES NOT EXIST')
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

l1 = '     Index           Old             New       Delta  Old-Design  New-Design'
nl=nl+1; lines(nl) = l1

n = size(v_array)
call re_allocate (lines, len(lines(1)), n+100)

do i = 1, size(v_array)
  var => v_array(i)%v
  if (.not. var%exists) cycle
  delta = var%model_value - var%old_value
  nl=nl+1; write (lines(nl), fmt) i, var%old_value, '  ->', &
       var%model_value, delta, var%old_value - var%design_value, &
	     var%model_value - var%design_value
enddo
nl=nl+1; lines(nl) = l1

if (max(abs(old_merit), abs(new_merit)) > 100) then
  fmt = '(5x, 2(a, es15.6), es15.6)'
else
  fmt = '(5x, 2(a, f13.6), f13.6)'
endif

nl=nl+1; lines(nl) = ' '
nl=nl+1; write (lines(nl), fmt) 'Merit:      ', old_merit, '  ->', new_merit
nl=nl+1; write (lines(nl), fmt) 'dMerit:     ', new_merit-old_merit
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
! Subroutine tao_change_ele (ele_name, attrib_name, num_str)
!
! Routine to change a variable in the model lattice.
!
! Input:
!   ele_name    -- Character(*): Name of variable or element.
!   attrib_name -- Character(*): Attribute name of element.
!   num_str     -- Character(*): Change in value. 
!                                A '@' signifies a absolute set.
!                                A 'd' signifies a set relative design.        
!
! Output:
!    %u(s%global%u_view)%model -- model lattice where the variable lives.
!-

subroutine tao_change_ele (ele_name, attrib_name, num_str)

use tao_mod
use quick_plot

implicit none

type (tao_universe_struct), pointer :: u
type (real_array_struct), allocatable, save :: d_ptr(:), m_ptr(:)

real(rp) change_number

character(*) ele_name, attrib_name, num_str
character(40) e_name, a_name, fmt
character(20) :: r_name = 'tao_change_ele'
character(200), allocatable, save :: lines(:)
character(20) abs_or_rel

real(rp) new_merit, old_merit, old_value, new_value, delta

integer i, ix, ix_a, nl
integer, allocatable :: ix_ele(:)

logical err

!-------------------------------------------------

u => s%u(s%global%u_view)

call re_allocate (lines, 200, 200)
call to_number (num_str, change_number, abs_or_rel, err);  if (err) return
old_merit = tao_merit()
nl = 0

! 

call string_trim (ele_name, e_name, ix)
call string_trim (attrib_name, a_name, ix)
call str_upcase (e_name, e_name)
call str_upcase (a_name, a_name)

call pointers_to_attribute (u%design%lat, e_name, a_name, .true., &
                                                d_ptr, err, .true., ix_ele, ix_a)
if (err) return

call pointers_to_attribute (u%model%lat, e_name, a_name, .true., &
                                                m_ptr, err, .true., ix_ele, ix_a)
if (err) return

do i = 1, size(d_ptr)
  if (ix_ele(i) < 0) cycle  ! bunch_start variables are always free
  if (.not. attribute_free (ix_ele(i), ix_a, u%model%lat, .true.)) return
end do

do i = 1, size(d_ptr)

    old_value = m_ptr(i)%r
     
    if (abs_or_rel == 'ABS') then
      m_ptr(i)%r = change_number
    elseif (abs_or_rel == 'REL') then
      m_ptr(i)%r = d_ptr(i)%r + change_number
    else
      m_ptr(i)%r = m_ptr(i)%r + change_number
    endif
     
    new_value = m_ptr(i)%r
  enddo

s%global%lattice_recalc = .true.

! don't print results if changing multiple elements

if (size(d_ptr) /= 1) return

!----------------------------------
! print results

new_merit = tao_merit()
delta = new_value - old_value
if (max(abs(old_value), abs(new_value), abs(d_ptr(1)%r)) > 100) then
  fmt = '(5x, 2(a, f12.0), f12.0)'
else
  fmt = '(5x, 2(a, f12.6), f12.6)'
endif

nl=nl+1;write (lines(nl), '(27x, a)') 'Old              New      Delta'
nl=nl+1;write (lines(nl), fmt) 'Value:       ', old_value, '  ->', new_value, delta
nl=nl+1;write (lines(nl), fmt) 'Value-Design:', old_value-d_ptr(1)%r, &
                                  '  ->', new_value-d_ptr(1)%r

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
