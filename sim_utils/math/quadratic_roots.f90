!+
! Function quadratic_roots(coefs) result (root)
!
! Routine to return the roots of the quadratic equation to machine precision.
! That is, avoiding the round-off error inherent in the standard formula.
!
! If coefs(3) is zero 0, root(2) is set to real_garbage$.
! If both coefs(2) and coefs(3) are zero, both roots are set to real_garbage$.
!
! When the roots are complex, the roots are ordered so that Imag(root(1)) < Imag(root(2)).
! When the roots are real, the roots are ordered so that root(1) <= root(2).
!
! Input:
!   coefs(:)  -- real(rp): Coefficients of the quadratic equation with 
!                 0 = coefs(1) + coefs(2) * x + coefs(3) * x^2
!
! Output:
!   root(2)   -- complex(rp): Complex roots. 
!-

function quadratic_roots(coefs) result (root)

use sim_utils_struct, only: rp, real_garbage$
implicit none

!

real(rp) coefs(3), rad, den
complex(rp) root(2)

!

if (coefs(3) == 0) then
  if (coefs(2) == 0) then
    root = real_garbage$
  else
    root(1) = -coefs(1) / coefs(2)
    root(2) = real_garbage$
  endif
  return
endif

rad = coefs(2)**2 - 4.0_rp * coefs(1) * coefs(3)
den = 0.5_rp / coefs(3)

if (rad < 0) then
  rad = sqrt(-rad)
  root(1) = den * cmplx(-coefs(2), -rad)
  root(2) = den * cmplx(-coefs(2), rad)

elseif (coefs(2) > 0) then
  rad = sqrt(rad)
  root(1) = den * (-coefs(2) - rad)
  root(2) = 2.0_rp * coefs(1) / (-coefs(2) - rad)

else
  rad = sqrt(rad)
  root(1) = den * (-coefs(2) + rad)
  root(2) = 2.0_rp * coefs(1) / (-coefs(2) + rad)
endif

end function
