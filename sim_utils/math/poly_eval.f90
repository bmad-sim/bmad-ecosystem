!+
! Function poly_eval(poly, x, diff_coef) result (y)
!
! Routine to evaluate a polynomial at a point.
!
! If diff_coef = False:
!   y = poly(0) + x * poly(1) + x^2 * poly(2) + ...
! If diff_coef = True:
!   y = poly(0)/0! + x * poly(1)/1! + x^2 * poly(2)/2! + ...
!
! Input:
!   poly(0:)    -- real(rp): Polynomial
!   x           -- real(rp): Point to evaluate at.
!   diff_coef   -- logical, optional: poly(:) array are differentials? Default is False.
!
! Output:
!   y           -- real(rp): Value of polynomial.
!-

function poly_eval(poly, x, diff_coef) result (y)

use sim_utils, dummy => poly_eval

implicit none

real(rp) poly(0:), x, y
integer i, n
logical, optional :: diff_coef

!

n = ubound(poly, 1)

if (logic_option(.false., diff_coef)) then
  y = poly(n) / factorial(n)
  do i = n-1, 0, -1
    y = y * x + poly(i) / factorial(i)
  enddo

else
  y = poly(n)
  do i = n-1, 0, -1
    y = y * x + poly(i)
  enddo
endif

end function
