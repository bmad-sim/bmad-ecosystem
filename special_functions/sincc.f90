!+
! function sincc(x, nd) result (y)
!
! Routine to calculate (x - sin(x)) / x^3 and derivatives.
!
! Input:
!   x  -- real(rp): Number.
!   nd -- integer, optional: Derivative order. nd = 0 (default) -> compute (x - sin(x)) / x^3
!           NOTE: Currently only nd = 0 and nd = 1 are implemented.
!
! Output:
!   y  -- real(rp): nd^th derivative of (x - sin(x)) / x^3
!-

elemental function sincc(x, nd) result (y)

use utilities_mod
implicit none

real(rp), intent(in) :: x
real(rp) y, x2
real(rp), parameter :: c0 = 1.0_rp / 6.0_rp, c1 = -1.0_rp / 120.0_rp, c2 = 1.0_rp / 5040.0_rp
real(rp), parameter :: c3 = -1.0_rp / 362880.0_rp

integer, optional, intent(in) :: nd

!

x2 = x*x

select case (integer_option(0, nd))
case (0)
  if (abs(x) < 0.1_rp) then
    y = c0 + x2 * (c1 + x2 * (c2 + x2 * c3))
  else
    y = (x - sin(x)) / (x*x2)
  end if

case (1)
  if (abs(x) < 0.1_rp) then
    y = x * (2.0_rp * c1 + x2 * (4.0_rp * c2 + 6.0_rp * x2 * c3))
  else
    y = (x * (1 - cos(x)) - 3.0_rp * (x - sin(x))) / (x2 * x2)
  endif

end select

end function
