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
character(80) line, v_name

logical limited

limited = .false.

! If the %correction_value is out of bounds then set the model_value so the
! target value is within bounds.

if (.not. s%global%var_limits_on) return
call tao_var_target_calc ()

do j = 1, size(s%var)

  var => s%var(j)

  if (var%model_value > var%high_lim) then
    write (line, '(1pe13.4)') var%high_lim
    v_name = trim(tao_var1_name(var)) // '  ' // trim(var%ele_name) // &
                                '[' // trim(var%attrib_name) // ']'
    call out_io (s_warn$, r_name, &
      'VARIABLE: ' // v_name, &
      'HAS TARGET VALUE GREATER THAN THE HIGH LIMIT OF: ' // line, &
      'RESETTING VARIABLE TO BE WITHIN BOUNDS & VETOING FROM OPTIMIZER LIST')
    value = var%high_lim
    call tao_set_var_model_value (var, value) 
    var%good_user = .false.
    limited = .true.
  endif

  if (var%model_value < var%low_lim) then
    write (line, '(1pe13.4)') var%low_lim
    v_name = trim(tao_var1_name(var)) // '  ' // trim(var%ele_name) // &
                                '[' // trim(var%attrib_name) // ']'
    call out_io (s_warn$, r_name, &
      'VARIABLE: ' // v_name, &
      'HAS TARGET VALUE LESS THAN THE LOW LIMIT OF: ' // line, &
      'RESETTING VARIABLE TO BE WITHIN BOUNDS & VETOING FROM OPTIMIZER LIST')
    value = var%low_lim
    call tao_set_var_model_value (var, value) 
    var%good_user = .false.
    limited = .true.
  endif

enddo

end subroutine
