!+
! Subroutine bbi_kick (x_norm, y_norm, r, kx, ky)
!
! This function computes, as a function of normalized position, the normalized kick due
! to the beam-beam interaction using the standard Bassetti-Erskine complex error function formula.
!
! The calculation assumes the kicking beam has opposite charge and has a Gaussian cross-section.
!
! Note: This routine can be used to compute the beam-ion kick as well.
!
! Input:
!   x_norm -- Real(rp): Normalized coordinate: x / sigma_x.
!   y_norm -- Real(rp): Normalized coordinate: y / sigma_y.
!   r -- Real(rp): Aspect ratio (sigma_y / sigma_x)
!
! Output:
!   kx -- Real(rp): Normalized, dimensionless horizontal kick component. 
!         In terms of the the actual kick (kick_x), kx is:
!           = kick_x / (xi_x * sigma_x / beta_x)
!           = -4 * pi * x /sigma_x  in the linear region
!
!   ky -- Real(rp): Normalized, dimensionless vertical kick component.
!         In terms of the the actual kick (kick_y), ky is:
!          = kick_y / (xi_y * sigma_y / beta_y)
!          = -4 * pi * y / sigma_y in the linear region
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
! In the calling routine, the formulas for computing the actual kicks, kick_x and kick_y, are:
!   kick_x = bbi_const * kx 
!          ~ -4 * pi * bbi_const * x / sigma_x                        [linear region]
!          ~ -2 * N_p * r_e * x / (gamma * sig_x * (sig_x + sig_y))   [linear region]
!          ~ -2 * N_p * r_e * x / (gamma * (x^2 + y^2))               [far from beam]
!   kick_y = bbi_const * ky
!          ~ -4 * pi * bbi_const * y / sigma_y                        [linear region]
!          ~ -2 * N_p * r_e * y / (gamma * sig_y * (sig_x + sig_y))   [linear region]
!          ~ -2 * N_p * r_e * y / (gamma * (x^2 + y^2))               [far from beam]
!
! For the beam-ion kick, assuming the ion velocity is neglegeble, the formulas are:
!   kick_x = ion_const * kx
!          ~ -4 * pi * ion_const * x / sigma_x                             [linear region]
!          ~ -2 * N_p * r_p * c_light * x / (sig_x * (sig_x + sig_y) * A)  [linear region]
!          ~ -2 * N_p * r_p * c_light * x / ((x^2 + y^2) * A)              [far from beam]
!   kick_y = ion_const * ky
!          ~ -4 * pi * ion_const * y / sigma_y                             [linear region]
!          ~ -2 * N_p * r_p * c_light * y / (sig_y * (sig_x + sig_y) * A)  [linear region]
!          ~ -2 * N_p * r_p * c_light * y / ((x^2 + y^2) * A)              [far from beam]
!
! where
!  ion_const = N_particles_bunch * r_p * c_light / (2 * pi * (sig_x + sig_y) * A)
!  A = Mass of ion in AMU.
!-

subroutine bbi_kick (x_norm, y_norm, r, kx, ky)

use precision_def
use physical_constants

implicit none

real(rp) x_norm, y_norm, r, xx, yy, rr, kx, ky, u, v, amp, denom, fr, fi, arg
real(rp) wr1, wi1, wr2, wi2, expon, scale
real(rp) :: emax = 30

! round beam case

if (r >= 0.998 .and. r <= 1.002) then ! round beam
  amp = (x_norm**2 + y_norm**2) / 2
  scale = 4 * pi
  if (amp > emax) then
    kx = -x_norm * scale / amp
    ky = -y_norm * scale / amp
  else if (amp < 1d-4) then
    kx = -x_norm * scale
    ky = -y_norm * scale
  else
    kx = -(x_norm * scale / amp) * (1 - exp(-amp))
    ky = -(y_norm * scale / amp) * (1 - exp(-amp))
  endif
  return
endif

! The calculation assumes that r < 1.
! If beam has r > 1 then switch x_norm, y_norm, and r

if (r > 1.0) then
  xx = y_norm
  yy = x_norm
  rr = 1 / r
else
  xx = x_norm
  yy = y_norm
  rr = r
endif

!

denom = 1.0/sqrt(2.0*(1 - rr**2))
scale = 4 * sqrt(pi**3) * denom * (1 + rr)
u = abs(xx) * denom
v = abs(yy) * denom

!********************************************************
!
!                              -(1-R^2)(V^2 +U^2)
!     F(U,V,R)  =  W(U+iRV) - E                  W(RU+iV)
!
!********************************************************

call complex_error_function(wr1, wi1, u, rr*v)
call complex_error_function(wr2, wi2, rr*u, v)

arg = (1 - rr**2) * (u**2 + v**2)

! if exponantial exponent is too small use only first term to
! evaluate KR and KI

if (arg > emax) then
  fr = wr1
  fi = wi1
else
  expon = exp(-arg)
  fr = wr1 - expon*wr2
  fi = wi1 - expon*wi2
endif

! if we have switched x_norm, y_norm, and r then switch kx for ky

if (r > 1.0) then  ! switch
  kx = -sign(fr, yy) * scale
  ky = -sign(fi, xx) * scale
else
  kx = -sign(fi, xx) * scale
  ky = -sign(fr, yy) * scale
endif

end subroutine
