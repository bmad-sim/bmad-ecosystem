!+
! function cosc(x, nd) result (y)
!
! Routine to calculate (1 - cos(x)) / x^2 and derivatives.
!
! Input:
!   x  -- real(rp):
!   nd -- integer, optional: Derivative order. nd = 0 (default) -> compute (1 - cos(x)) / x^2
!           NOTE: Currently only nd = 0 and nd = 1 are implemented.
!
! Output:
!   y  -- real(rp): nd^th derivative of (1 - cos(x)) / x^2
!-

elemental function cosc(x, nd) result (y)

use utilities_mod
use sim_utils_interface, dummy => cosc
implicit none

real(rp), intent(in) :: x
real(rp) y

integer, optional, intent(in) :: nd

!

select case (integer_option(0, nd))
case (0)
  if (abs(x) < 1e-8_rp) then
    y = 0.50_rp
  else
    y =  0.5_rp * sinc(x/2)**2
  end if

case (1)
  y =  0.5_rp * sinc(x/2) * sinc(x/2, 1)
end select

end function
