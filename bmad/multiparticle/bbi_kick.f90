!+
! Subroutine bbi_kick (x, y, sigma, nk, dnk, linear_kick)
!
! This function computes, as a function of normalized position, the normalized kick due
! to the beam-beam interaction using the standard Bassetti-Erskine complex error function formula.
!
! The calculation assumes the kicking beam has opposite charge and has a Gaussian cross-section.
!
! Note: This routine can be used to compute the beam-ion kick as well.
!
! Input:
!   x           -- real(rp): X coordinate.
!   y           -- real(rp): Y coordinate.
!   sigma(2)    -- real(rp): Beam (x,y) sigmas.
!   linear_kick -- logical, optional: Default False. If present and True, kick and dnk are computed
!                       using the extrapolated kick from the linear region.
!
! Output:
!   nk(2)     -- real(rp): Normalized, dimensionless kick component. 
!                In terms of the the actual kick:
!                nk = [kick_x / (xi_x * sigma_x / beta_x), kick_y / (xi_y * sigma_y / beta_y)
!                nk = -4 * pi * [x/sigma_x, y/sigma_y] in the linear region
!   dnk(2,2)  -- real(rp): derivatives of nk. EG: dnk(2,1) = nk(2)/dx
!
! Note:
!   xi_x = beta_x * bbi_const / sig_x     ! Horizontal tune shift parameter
!   xi_y = beta_y * bbi_const / sig_y     ! Vertical   tune shift parameter
! where
!   bbi_const = N_particles_bunch * r_e / (2 * pi * gamma * (sig_x + sig_y))
! And the tune shifts are:
!   dQ_x = xi_x = beta_x * bbi_const / sig_x
!   dQ_y = xi_y = beta_y * bbi_const / sig_y
!
! In the calling routine, the formulas for computing the actual kicks, kick_x and kick_y, should be:
!   kick_x = bbi_const * nk(1) 
!          ~ -4 * pi * bbi_const * x / sig_x                          [linear region]
!          ~ -2 * N_p * r_e * x / (gamma * sig_x * (sig_x + sig_y))   [linear region]
!          ~ -2 * N_p * r_e * x / (gamma * (x^2 + y^2))               [far from beam]
!   kick_y = bbi_const * nk(2)
!          ~ -4 * pi * bbi_const * y / sig_y                          [linear region]
!          ~ -2 * N_p * r_e * y / (gamma * sig_y * (sig_x + sig_y))   [linear region]
!          ~ -2 * N_p * r_e * y / (gamma * (x^2 + y^2))               [far from beam]
!
! For the beam-ion kick, assuming the ion velocity is neglegeble, the formulas are:
!   kick_x = ion_const * nk(1)
!          ~ -4 * pi * ion_const * x / sig_x                               [linear region]
!          ~ -2 * N_p * r_p * c_light * x / (sig_x * (sig_x + sig_y) * A)  [linear region]
!          ~ -2 * N_p * r_p * c_light * x / ((x^2 + y^2) * A)              [far from beam]
!   kick_y = ion_const * nk(2)
!          ~ -4 * pi * ion_const * y / sig_y                               [linear region]
!          ~ -2 * N_p * r_p * c_light * y / (sig_y * (sig_x + sig_y) * A)  [linear region]
!          ~ -2 * N_p * r_p * c_light * y / ((x^2 + y^2) * A)              [far from beam]
!
! where
!   ion_const = N_particles_bunch * r_p * c_light / (2 * pi * (sig_x + sig_y) * A)
!   A = Mass of ion in AMU.
!-

subroutine bbi_kick (x, y, sigma, nk, dnk, linear_kick)

use precision_def
use physical_constants
use sign_of_mod

implicit none

real(rp) x, y, x_norm, y_norm, r, nk(2), dnk(2,2), u, v, amp, denom, f(2), arg, dscale, damp_dx, damp_dy
real(rp) w1(2), w2(2), dw1(2,2), dw2(2,2), expon, scale, sigma(2), sx, sy, sig(2), expa
real(rp) :: emax = 30
logical, optional :: linear_kick
logical flipped

!

x_norm = x / sigma(1)
y_norm = y / sigma(2)

! Linear kick

if (logic_option(.false., linear_kick)) then
  scale = 4 * pi
  nk = -[x_norm, y_norm] * scale
  dnk(1,:) = [-scale/sigma(1), 0.0_rp]
  dnk(2,:) = [0.0_rp, -scale/sigma(2)]
  return
endif

! Round beam case

r = sigma(2)/sigma(1)

if (r > 0.999 .and. r < 1.001) then ! round beam
  amp = (x_norm**2 + y_norm**2) / 2
  if (amp > emax) then
    scale = 4 * pi / amp
    dscale = -4 * pi / amp**2
  else if (amp < 1d-4) then
    scale = twopi * (2.0_rp - amp)
    dscale = -twopi
  else
    expa = exp(-amp)
    scale = 4 * pi * ((1 - expa) / amp)
    dscale = 4 * pi * (amp*expa + expa - 1.0_rp) / amp**2
  endif

  damp_dx = x_norm / sigma(1)
  damp_dy = y_norm / sigma(2)
  nk = -[x_norm, y_norm] * scale
  dnk(1,:) = -[scale/sigma(1) + x_norm*dscale*damp_dx, x_norm*dscale*damp_dy]
  dnk(2,:) = -[y_norm*dscale*damp_dx,                  scale/sigma(2) + y_norm*dscale*damp_dy]
  return
endif

! Equations have been developed assuming r < 1.
! If not true, switch x and y.

flipped = (r > 1)
if (flipped) then
  r = 1/r
  x_norm = y_norm
  y_norm = x / sigma(1)
  sig = [sigma(2), sigma(1)]
else
  sig = sigma
endif

!********************************************************
!
!                              -(1-R^2)(U^2+V^2)
!     F(U,V,R)  =  W(U+iRV) - E                 W(RU+iV)
!
!********************************************************

denom = 1.0 / sqrt(2.0*(1 - r**2))
scale = 4 * sqrt(pi**3) * denom * (1 + r)
u = abs(x_norm) * denom
v = abs(y_norm) * denom
sx = sign_of(x_norm)
sy = sign_of(y_norm)
arg = (1 - r**2) * (u**2 + v**2)

call faddeeva_function([u, r*v], w1, dw1)
call faddeeva_function([r*u, v], w2, dw2)

expon = exp(-arg)
f = w1 - expon*w2

nk = -scale * [sx*f(2), sy*f(1)]

dnk(1,:) =           [dw1(2,1)/sig(1), dw1(2,2)*r/sig(2)] - &
             expon * [dw2(2,1)*r/sig(1), dw2(2,2)/sig(2)] + &
             2.0_rp * (1 - r**2) * w2(2) * expon * [u/sig(1), v/sig(2)]
dnk(2,:) =           [dw1(1,1)/sig(1), dw1(1,2)*r/sig(2)] - &
             expon * [dw2(1,1)*r/sig(1), dw2(1,2)/sig(2)] + &
             2.0_rp * (1 - r**2) * w2(1) * expon * [u/sig(1), v/sig(2)]

dnk(1,:) = -scale * denom * [dnk(1,1), dnk(1,2)*sx*sy]
dnk(2,:) = -scale * denom * [dnk(2,1)*sx*sy, dnk(2,2)]

if (flipped) then
  nk = [nk(2), nk(1)]
  dw1 = dnk
  dnk(1,:) = [dw1(2,2), dw1(2,1)]
  dnk(2,:) = [dw1(1,2), dw1(1,1)]
endif

end subroutine
