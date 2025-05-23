!+
! Subroutine tao_to_real (expression, value, err_flag)
!
! Mathematically evaluates an expression.
!
! Input:
!   expression    -- character(*): arithmetic expression
!
! Output:
!   value        -- real(rp): Value of arithmetic expression.
!   err_flag     -- Logical: TRUE on error.
!-

subroutine tao_to_real (expression, value, err_flag)

use tao_interface, dummy => tao_to_real

implicit none

character(*) :: expression

real(rp) value
real(rp), allocatable :: vec(:)

logical err_flag

!

call tao_evaluate_expression (expression, 1, .false., vec, err_flag)
if (err_flag) return
value = vec(1)

end subroutine tao_to_real

