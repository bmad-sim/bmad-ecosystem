!+
! Subroutine linear_fit_2d (x, y, z, coef, covar)
!
! Subroutine to fit to z = coef(1) x + coef(2) y + coef(3)
!
! Input:
!   x(:)        -- real(rp): Array of x-values.
!   y(:)        -- real(rp): Array of y-values.
!   z(:)        -- real(rp): Array of z-values
!
! Output:
!   coef(3)     -- real(rp): Coefficients of the linear fit
!   covar(3,3)  -- real(rp), optional: Covariance matrix.
!-

subroutine linear_fit_2d (x, y, z, coef, covar)

use sim_utils, dummy => linear_fit_2d

implicit none

integer n, n_data

real(rp) x(:), y(:), z(:), coef(3)
real(rp), optional :: covar(3,3)
real(rp) m(3,3), v(3), m_inv(3,3), sig2, det

!

n_data = size(x)

m(1,:) = [sum(x*x), sum(x*y), sum(x)]
m(2,:) = [  m(1,2), sum(y*y), sum(y)]
m(3,:) = [  m(1,3),   m(2,3), real(n_data, rp)] 

v = [sum(x*z), sum(y*z), sum(z)]

!

det = determinant(m)
if (det == 0) then
  coef = real_garbage$
  covar = real_garbage$
  return
endif

call mat_inverse(m, m_inv)
coef = matmul(m_inv, v)

!

if (present(covar)) then
  sig2 = sum((z - x*coef(1) - y*coef(2) - coef(3))*(z - x*coef(1) - y*coef(2) - coef(3))) / n_data
  covar = sig2 * m_inv * transpose(m_inv)
endif

end subroutine

