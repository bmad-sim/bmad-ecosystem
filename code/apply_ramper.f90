!+
! Subroutine apply_ramper (slave, ramper, err_flag)
!
! Input:
!   slave       -- ele_struct: Element to apply ramper to.
!   ramper      -- ele_struct: Ramper element.
!
! Output:
!   err_flag    -- logical: Set true if there is an error. False otherwise.
!-

subroutine apply_ramper (slave, ramper, err_flag)

use bmad, except_dummy => apply_ramper

implicit none

type (ele_struct), target :: slave, ramper
type (control_struct), pointer :: c
type (all_pointer_struct) a_ptr

integer iv, key, ix
logical err_flag, ok

character(100) err_str
character(40) name
character(*), parameter :: r_name = 'apply_ramper'

! 

err_flag = .false.

do iv = 1, size(ramper%control%ramp)
  c => ramper%control%ramp(iv)

  ix = index(c%slave_name, '::')
  if (ix == 0) then
    key = 0
    name = c%slave_name
  else
    key = key_name_to_key_index(c%slave_name(1:ix-1), .true.)
    name = c%slave_name(ix+2:)
  endif

  if ((key /= 0 .and. key /= slave%key) .or. .not. match_wild(slave%name, name)) then
    c%value = real_garbage$
    cycle
  endif

  call pointer_to_attribute (slave, c%attribute, .true., a_ptr, err_flag, .false.)
  if (err_flag .or. .not. associated(a_ptr%r)) then
    c%value = real_garbage$
    cycle
  endif

  if (ramper%control%type == expression$) then
    call evaluate_expression_stack(c%stack, c%value, err_flag, err_str, ramper%control%var, .false.)
    if (err_flag) then
      call out_io (s_error$, r_name, err_str, ' OF RAMPER: ' // ramper%name)
      err_flag = .true.
      return
    endif

  else
    call spline_akima_interpolate (ramper%control%x_knot, c%y_knot, ramper%control%var(1)%value, ok, c%value)
    if (.not. ok) then
      call out_io (s_error$, r_name, 'VARIABLE VALUE OUTSIDE OF SPLINE KNOT RANGE.')
      err_flag = .true.
      return
    endif
  endif

  a_ptr%r = c%value
  call set_flags_for_changed_attribute (slave, a_ptr%r, .true.)
enddo

call attribute_bookkeeper(slave, .true.)

end subroutine
