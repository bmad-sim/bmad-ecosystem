!+
! Subroutine apply_ramper (ele, ramper, err_flag)
!
! Input:
!   ele         -- ele_struct: Element to apply ramper to.
!   ramper      -- ele_struct: Ramper element.
!
! Output:
!   start_orb   -- coord_struct: Modified starting position for track1 to use.
!   err_flag    -- logical: Set true if there is an error. False otherwise.
!-

subroutine apply_ramper (ele, ramper, err_flag)

use bmad, except_dummy => apply_ramper

implicit none

type (ele_struct), target :: ele, ramper
type (control_struct), pointer :: c
integer iv
logical err_flag, ok

character(100) err_str
character(*), parameter :: r_name = 'apply_ramper'

!

err_flag = .false.


do iv = 1, size(ramper%control%ramp)
  c => ramper%control%ramp(iv)

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

enddo


end subroutine
