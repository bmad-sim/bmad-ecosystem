!+
! Function ramper_value (ramper, r1, err_flag) result (value)
!
! Low level routine to evaluate one slave function of a ramper lord.
!
! Input:
!   ramper        -- ele_struct: Ramper lord.
!   r1            -- control_ramp1_struct: Slave function.
!
! Output:
!   err_flag      -- logical: Set True if there is an error, False otherwise.
!   value         -- real(rp): Value of the slave function.
!-

function ramper_value (ramper, r1, err_flag) result (value)

use bmad_routine_interface, dummy => ramper_value
use expression_mod, only: expression_stack_value

implicit none

type (ele_struct) ramper
type (control_ramp1_struct) r1

real(rp) value
logical err_flag
character(100) err_str
character(*), parameter :: r_name = 'ramper_value'

!

if (allocated(r1%stack)) then
  value = expression_stack_value(r1%stack, err_flag, err_str, ramper%control%var, .false.)
  if (err_flag) then
    call out_io (s_error$, r_name, err_str, ' OF RAMPER: ' // ramper%name)
    err_flag = .true.
    return
  endif

elseif (allocated(r1%y_knot)) then
  value = knot_interpolate(ramper%control%x_knot, r1%y_knot, &
          ramper%control%var(1)%value, nint(ramper%value(interpolation$)), err_flag)
  if (err_flag) then
    call out_io (s_error$, r_name, 'VARIABLE VALUE (\es12.4\) OF RAMPER ELEMENT: ' // ramper%name, &
                                   'IS OUTSIDE OF SPLINE KNOT RANGE OF SLAVE: ' // r1%slave_name)
    return
  endif
endif

err_flag = .false.

end function ramper_value
