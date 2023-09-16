!+
! function asinc(x, nd) result (y)
!
! Routine to calculate arcsin(x) / x and derivatives.
!
! Input:
!   x  -- real(rp):
!   nd -- integer, optional: Derivative order. nd = 0 (default) -> compute arcsin(x) / x
!           NOTE: Currently only nd = 0 and nd = 1 are implemented.
!
! Output:
!   y  -- real(rp): nd^th derivative. of arcsin(x)/x
!-

elemental function asinc(x, nd) result (y)

use utilities_mod
implicit none

real(rp), intent(in) :: x
real(rp) y, x2
real(rp), parameter :: c1 = 1.0_rp / 3.0_rp, c2 = 0.3_rp
real(rp), parameter :: c3 = 15.0_rp / 56.0_rp, c4 = 35.0_rp / 144.0_rp

integer, optional, intent(in) :: nd

!

x2 = x*x

select case (integer_option(0, nd))
case (0)
  if (abs(x) < 1d-8) then
     y = 1.0_rp
  else
     y = asin(x)/x
  end if

case (1)
  if (abs(x) < 0.03_rp) then
    y = x * (c1 + x2 * (c2 + x2 * (c3 + x2 * c4)))
  else
    y = 1 / (x * sqrt(1 - x2)) - asin(x) / x2
  endif
end select

end function
