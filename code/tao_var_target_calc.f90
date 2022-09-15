!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_var_target_calc ()
! 
! Subroutine to calculate the variable target values (the values that they need
! to be set to to do a correction of the orbit, phase, etc.
!
! Input:
!
! Output:
!-

subroutine tao_var_target_calc ()

use tao_struct

implicit none

type (tao_var_struct), pointer :: var

integer i, j

!

do j = 1, s%n_var_used
  var => s%var(j)
  if (.not. var%exists) cycle
  var%correction_value = var%meas_value + (var%design_value - var%model_value)
enddo

end subroutine tao_var_target_calc
