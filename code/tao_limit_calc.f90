!+
! Subroutine tao_limit_calc (s, limited)
!
! Subroutine to make sure the target value of a variable does not go outside 
! the variables's limits. If a variable has a target_value outside a limit
! the variable's model_value is changed so that the target_value is inside the
! limit. Additionally the variable's good_user attribute is set to False.
!
! Input:
!   s       -- Tao_super_universe_struct:
!    %globalvar_limits_on -- Logical: If False then this routine does nothing.
!
! Output:
!   s       -- Tao_super_universe_struct:
!   limited -- Logical: Set True if a variable is past a limit.
!-

subroutine tao_limit_calc (s, limited)

use tao_mod

implicit none

type (tao_super_universe_struct), target :: s
type (tao_var_struct), pointer :: var

integer i, j

character(20) :: r_name = 'tao_limit_calc'
character(80) line

logical limited

! If the %target_value is out of bounds then set the model_value so the
! target value is within bounds.

if (.not. s%global%var_limits_on) return
call tao_var_target_calc (s)

limited = .false.

do i = 1, size(s%u)
  do j = 1, size(s%u(i)%var)

    var => s%u(i)%var(j)

    if (var%target_value > var%high_lim_value) then
      write (line, '(a, 1pe13.4)') &
       'HAS TARGET VALUE GREATER THAN THE HIGH LIMIT OF: ', var%high_lim_value
      call out_io (s_warn$, r_name, 'VARIABLE: ' // var%name, line, &
        'RESETTING VARIABLE TO BE WITHIN BOUNDS & VETOING FROM OPTIMIZER LIST')
      var%model_value = var%model_value - &
                  1.00001 * (var%target_value - var%high_lim_value)
      var%good_user = .false.
      limited = .true.
    endif

    if (var%target_value < var%low_lim_value) then
      write (line, '(a, 1pe13.4)') &
         'HAS TARGET VALUE LESS THAN THE LOW LIMIT OF: ', var%low_lim_value
      call out_io (s_warn$, r_name, 'VARIABLE: ' // var%name, line, &
        'RESETTING VARIABLE TO BE WITHIN BOUNDS & VETOING FROM OPTIMIZER LIST')
      var%model_value = var%model_value + &
                  1.00001 * (var%low_lim_value - var%target_value)
      var%good_user = .false.
      limited = .true.
    endif

  enddo
enddo

end subroutine
