!+
! faddeeva_function(z, w, dw)
!
! Routine to calculate the Faddeeva function and derivates.
! The Faddeeva function is typically denoted w(z) with
!   w(z) = exp(-z^2) * erfc(-i*z)
! where erfc is the complementary error function.
!
! The faddeeva function is also called the "complex error function" in the literature.
!
! Input:
!   z(2)        -- real(rp): z = (x,y) vector to evaluate the Faddeeva function at.
!
! Output:
!   w(2)        -- real(rp): Faddeeva function at z.
!   dw(2,2)     -- real(rp): Faddeeva derivatives. Eg: dw(2,1) = dw_y/dx.
!-

subroutine faddeeva_function(z, w, dw)

use precision_def
use physical_constants, only: pi
use, intrinsic :: iso_c_binding

implicit none

real(rp) z(2), w(2), dw(2,2)
real(rp) :: sqrt1pi = 1 / sqrt(pi)
real(c_double) :: errmax = 0    ! Maximum error
real(c_double) x, y, wx, wy

interface
  subroutine faddeeva_w(x, y, wx, wy, errmax) bind(c)
    import
    real(c_double), value :: x, y, errmax
    real(c_double) :: wx, wy
  end subroutine
end interface

!

x = z(1); y = z(2)
call faddeeva_w(x, y, wx, wy, errmax)
w = [wx, wy]

dw(1,:) = [2.0_rp * (z(2)*w(2) - z(1)*w(1)), -2.0_rp * (sqrt1pi - z(1)*w(2) - z(2)*w(1))]
dw(2,:) = [-dw(1,2), dw(1,1)]

end subroutine
