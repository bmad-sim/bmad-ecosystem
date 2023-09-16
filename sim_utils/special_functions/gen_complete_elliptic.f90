!+
! Function gen_complete_elliptic(kc, p, c, s, err_tol) result (value)
!
! Routine to calculate the generalized complete elliptic integral.
!
! Algorithm is by Bulirsch. 
! Code is taken from:
!   "Cylindrical Magnets and Ideal Solenoids"
!   Norman Derby & Stanislaw Olbert
!   https://arxiv.org/pdf/0909.3880.pdf
!
! Conversions:
! The complete elliptic integral of the first kind:
!   K(k) = gen_complete_elliptic((kc, 1, 1, 1)
! The complete elliptic integral of the second kind:
!   E(k) = = gen_complete_elliptic(kc, 1, 1, kc^2)
! The complete elliptic integral of the third kind:
!   Pi(n, k) = = gen_complete_elliptic(kc, n+1, 1, 1)
! Where:
!   k^2 = 1 - kc^2
!
! Input:
!   kc, p, c, s -- real(rp): Fuction input values.
!   err_tol     -- real(rp), optional: Relative error tolerance. Default = 1d-12 
!
! Output:
!   value       -- real(rp): Output value.
!-

function gen_complete_elliptic(kc, p, c, s, err_tol) result (value)

use sim_utils_interface, dummy => gen_complete_elliptic
use utilities_mod

implicit none

real(rp) kc, p, c, s, value
real(rp), optional :: err_tol
real(rp) k, pp, cc, ss, em, f, q, g, kk, errtol

!

if (kc == 0) then
  value = real_garbage$
  return
endif

errtol = real_option(1e-12_rp, err_tol)

k = abs(kc)
pp = p
cc = c
ss = s
em = 1

if (p > 0) then
  pp = sqrt(p)
  ss = s/pp
else
  f = kc*kc
  q = 1. - f
  g = 1. - pp
  f = f - pp
  q = q*(ss - c*pp)
  pp = sqrt( f/g )
  cc = (c - ss)/g
  ss = - q/(g*g*pp) + cc*pp
endif

f = cc
cc = cc + ss/pp
g = k/pp
ss = 2*(ss + f*g)
pp = g + pp
g = em
em = k + em
kk = k

do while (abs(g - k) > g*errtol)
 k = 2*sqrt(kk)
 kk = k*em
 f = cc
 cc = cc + ss/pp
 g = kk/pp
 ss = 2*(ss + f*g)
 pp = g + pp
 g = em
 em = k + em
enddo

value = (pi/2.0_rp)*(ss + cc*em)/( em*(em + pp) )

end function
