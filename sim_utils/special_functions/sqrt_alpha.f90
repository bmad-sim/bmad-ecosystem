!+
! function sqrt_alpha(alpha, x) result (y)
!
! Routine to calculate to machine precision:
!   y = Sqrt[alpha^2 + x] - alpha 
! Also see the sqrt_one() function.
!
! This is usful if x is near zero where the direct evaluation of y is inaccurate.
!
! Input:
!   alpha   -- real(rp): Number
!   x       -- real(rp): Number
!
! Output:
!   y       -- real(rp): Result.
!-

elemental function sqrt_alpha(alpha, x) result (y)

use utilities_mod
implicit none

real(rp), intent(in) :: alpha, x
real(rp) y

!

y = x / (sqrt(alpha*alpha + x) + alpha)

end function
