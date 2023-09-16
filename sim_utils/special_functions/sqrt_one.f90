!+
! function sqrt_one(x, nd) result (y)
!
! Routine to calculate Sqrt[1+x] - 1 and derivatives to machine precision.
!
! This is usful if x can be near zero where the direct evaluation of
! sqrt[1+x] - 1 is inaccurate.
!
! Input:
!   x   -- real(rp): Number
!   nd  -- integer, optional: Derivative order. nd = 0 (default) -> compute Sqrt[1+x] - 1.
!           NOTE: Currently only nd = 0 and nd = 1 are implemented.
!
! Output:
!   y   -- real(rp): Result.
!-

elemental function sqrt_one(x, nd) result (y)

use utilities_mod
implicit none

real(rp), intent(in) :: x
real(rp) y, rad, sq

integer, optional, intent(in) :: nd

!

sq = sqrt(1 + x)
rad = sq + 1

select case (integer_option(0, nd))
case (0)
  y = x / rad

case (1)
  y = 1 / rad - x / (2 * rad**2 * sq)
end select

end function
