!+
! Subroutine tao_limit_calc (limited)
!
! Subroutine to make sure the target value of a variable does not go outside 
! the variables's limits. If a variable has a correction_value outside a limit
! the variable's model_value is changed so that the correction_value is inside the
! limit. Additionally the variable's good_user attribute is set to False.
!
! Input:
!    %globalvar_limits_on -- Logical: If False then this routine does nothing.
!
! Output:
!   limited -- Logical: Set True if a variable is past a limit.
!-

subroutine tao_limit_calc (limited)

use tao_mod

implicit none

type (tao_var_struct), pointer :: var

real(rp) value
integer i, j

character(20) :: r_name = 'tao_limit_calc'
character(80) line

logical limited

! If the %correction_value is out of bounds then set the model_value so the
! target value is within bounds.

if (.not. s%global%var_limits_on) return
call tao_var_target_calc ()

limited = .false.

do j = 1, size(s%var)

  var => s%var(j)

  if (var%correction_value > var%high_lim) then
    write (line, '(a, 1pe13.4)') &
     'HAS TARGET VALUE GREATER THAN THE HIGH LIMIT OF: ', var%high_lim
    call out_io (s_warn$, r_name, 'VARIABLE: ' // var%name, line, &
      'RESETTING VARIABLE TO BE WITHIN BOUNDS & VETOING FROM OPTIMIZER LIST')
    value = var%model_value - 1.00001 * (var%correction_value - var%high_lim)
    call tao_set_var_model_value (var, value) 
    var%good_user = .false.
    limited = .true.
  endif

  if (var%correction_value < var%low_lim) then
    write (line, '(a, 1pe13.4)') &
       'HAS TARGET VALUE LESS THAN THE LOW LIMIT OF: ', var%low_lim
    call out_io (s_warn$, r_name, 'VARIABLE: ' // var%name, line, &
      'RESETTING VARIABLE TO BE WITHIN BOUNDS & VETOING FROM OPTIMIZER LIST')
    value = var%model_value + 1.00001 * (var%low_lim - var%correction_value)
    call tao_set_var_model_value (var, value) 
    var%good_user = .false.
    limited = .true.
  endif

enddo

end subroutine
