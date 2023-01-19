!+
! Subroutine n_spline_create (deriv0, deriv1, x1, n_spline)
!
! Routine to create a non-smoothing polynomial spline interpolating function given derivative 
! vectors deriv0(0:n) and deriv1(0:n) at (knot) points x = 0 and x1.
!
! The interpolating spline will be degree 2*n+1:
!   y_spline(x) = sum_{p in range [0, 2n+1]} n_spline(p) * x^p / factorial(p)
! For n = 1, the n_spline is a cubic spline.
!
! The spline is constructed so that at x = 0 and x1 the derivatives up to order n agree with 
! deriv0 and deriv1 respectively:
!   d^p y_spline / dx^p (0)  = deriv0(p),
!   d^p y_spline / dx^p (x1) = deriv1(p)   For all p in range [0, n]
!
! To evaluate the p^th derivative of the interpolation spline use:
!   p_deriv = poly_eval(n_spline(p:), x_eval, diff_coef = .true.)
!
! Input:
!   deriv0(0:)    -- real(rp): Derivative vector from order 0 to some order n at x = 0.
!   deriv1(0:)    -- real(rp): Derivative vector from order 0 to some order n at x = x1.
!   x1            -- real(rp): Location where deriv1 derivatives have been evaluated.
!
! Output:
!   n_spline(0:)  -- real(rp), Derivative vector from order 0 to order 2*n+1 of the interpolation spline.
!-

subroutine n_spline_create (deriv0, deriv1, x1, n_spline)

use sim_utils, dummy => n_spline_create

implicit none

real(rp) deriv0(0:), deriv1(0:)
real(rp) x1
real(rp) :: n_spline(0:)
real(rp), allocatable :: mat(:,:), v(:), x1p(:)

integer n, j, k

! 

n = ubound(deriv0,1)
allocate(mat(0:n,0:n), v(0:n), x1p(0:2*n+1))

x1p(0) = 1
forall (j = 1:2*n+1) x1p(j) = x1 * x1p(j-1)

do j = 0, n
  v(j) = deriv1(j) - poly_eval(deriv0(j:), x1, diff_coef = .true.)
  do k = 0, n
    mat(j,k) = x1**(n+1+k-j) / factorial(n+1+k-j)
  enddo
enddo

call mat_inverse (mat, mat)
n_spline(0:n) = deriv0
n_spline(n+1:2*n+1) = matmul(mat, v)

end subroutine
