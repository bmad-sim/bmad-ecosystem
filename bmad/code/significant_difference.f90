!+
! Function significant_difference (value1, value2, abs_tol, rel_tol) result (is_different)
!
! Routine to determine if two values are significantly different.
! There is a significant difference if:
!   abs(value1 - value2) > abs_tol + rel_tol * (abs(value1) + abs(value(2)))
!
! Input:
!   value1    -- real(rp): First value.
!   value2    -- real(rp): Second value.
!   abs_tol   -- real(rp), optional: Absolute tolerance. Default is 0.
!   rel_tol   -- real(rp), optional: Relative tolerance. Default is 0.
!
! Output:
!   is_different -- logical: Set True if the difference is significant. False otherwise.
!-

function significant_difference (value1, value2, abs_tol, rel_tol) result (is_different)

use utilities_mod

implicit none

real(rp), intent(in) :: value1, value2
real(rp), intent(in), optional :: abs_tol, rel_tol

logical is_different

!

is_different = (abs(value1 - value2) > real_option(0.0_rp, abs_tol) + real_option(0.0_rp, rel_tol) * (abs(value1) + abs(value2)))

end function significant_difference
