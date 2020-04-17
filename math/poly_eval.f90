!+
! Function poly_eval(poly, x) result (y)
!
! Routine to evaluate a polynomial at a point.
!   y = poly(0) + x * poly(1) + x^2 * poly(2) + ...
!
! Input:
!   poly(0:)    -- real(rp): Polynomial
!   x           -- real(rp): Point to evaluate at.
!
! Output:
!   y           -- real(rp): Value of polynomial.
!-

function poly_eval(poly, x) result (y)

use precision_def

implicit none

real(rp) poly(0:), x, y
integer i, n

!

n = ubound(poly, 1)
y = poly(n)
do i = n-1, 0, -1
  y = y * x + poly(i)
enddo

end function
