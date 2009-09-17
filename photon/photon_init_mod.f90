module photon_init_mod

use bmad_struct
use bmad_interface
use nr, only: erf_s
use random_mod

integer, private, parameter :: sigma_cut = 5, n_pts_per_sigma = 20
integer, private, parameter :: n_g = sigma_cut*n_pts_per_sigma
real(rp), private, save :: s1, s2, c1, c2
real(rp), private, save :: g1(0:n_g), g2(0:n_g)

contains

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!+
! Subroutine photon_init (g_bend, gamma, orbit)
!
! Routine to initalize a photon
!
! Modules needed:
!   use photon_init_mod
!
! Input:
!   g_bend -- Real(rp): Bending 1/rho.
!   gamma  -- Real(rp): Relativistic beam gamma factor.
!
! output:
!   orbit(6) -- Real(rp): Orbit at the emission point.
!-

subroutine photon_init (g_bend, gamma, orbit)

implicit none

real(rp) g_bend, gamma, orbit(6)
real(rp) E_rel, gamma_phi

!

call photon_energy_init (E_rel)
call photon_vert_angle_init (E_rel, gamma_phi)

orbit = 0
orbit(4) = gamma_phi / gamma
orbit(6) = E_rel * 3 * h_bar_planck * c_light * gamma**3 * g_bend / 2

end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!+
! Subroutine photon_vert_angle_init (E_rel, gamma_phi, r_in)
!
! Routine to convert a "random" number in the interval [0,1] to a photon vertical emission 
! angle for a simple bend.
! The verital angle distribution is modeled to be of the form:
!   Prob(x) =  c1 * e^(-x^2/(2*s1^2)) + c2 * x^2 * e^(-x^2/(2*s2^2))
! where x = gamma * phi, and c1, c2, s1, and s2 are functions of E_rel.
! The c1 * e^(-x^2/(2*s1^2)) term models the perpendicular radiation component.
! The c2 * e^(-x^2/(2*s2^2)) term models the parallel radiation component.
!
! Module needed:
!   use photon_init_mod
!
! Input:
!   r_in  -- Real(rp), optional: number in the range [0,1].
!             If not present, a random number will be used.
!   E_rel -- Real(rp): Relative photon energy E/E_crit. 
! 
! Output:
!   gamma_phi -- Real(rp): gamma * phi where gamma is the beam relativistic factor and
!                 phi is the vertical photon angle (in radians).
!-

subroutine photon_vert_angle_init (E_rel, gamma_phi, r_in)

implicit none

real(rp), optional :: r_in
real(rp) e_rel, gamma_phi
real(rp) rr, r, ss, x, log_E, frac

real(rp) :: rel_amp(0:32) = (/ &
                  0.249546, 0.249383, 0.249162, 0.248862, 0.248454, &
                  0.247902, 0.247154, 0.246141, 0.244774, 0.242931, &
                  0.240455, 0.237143, 0.232737, 0.226921, 0.219325, &
                  0.209542, 0.197181, 0.181949, 0.163786, 0.143015, &
                  0.120465, 0.0974421, 0.0755035, 0.0560654, 0.0400306, &
                  0.0276386, 0.018578, 0.0122386, 0.00794693, 0.00510924, &
                  0.00326308, 0.00207494, 0.0013157 /)

integer i, ix
logical, save :: init_needed = .true.

! Init

if (init_needed) then
  do i = 0, n_g-1
    x = real(i, rp) / n_pts_per_sigma
    g1(i) = 0.5 * erf_s (x/sqrt_2)
    g2(i) = 0.5 * erf_s (x/sqrt_2) - x * exp(-x**2/2) / sqrt(twopi)
  enddo
  g1(n_g) = 0.50000001
  g2(n_g) = 0.50000001
  init_needed = .false.
endif

! The following was derived from a fit to the standard synchrotron radiation emissian equation.

! Sigmas of the parallel (s1) and perpendicular (s2) components.
! The sigmas are a function of E_rel.

if (E_rel < 0.1) then
  s1 = 0.667734 * E_rel**(-0.336135)
  s2 = 0.545472 * E_rel**(-0.334665)
else
  s1 = 0.524027 * E_rel**(-0.449124)
  s2 = 0.470202 * E_rel**(-0.413415)
endif

! Parallel (c1) and perpendicular (c2) weights.
! the array, rel_amp, gives the relative amplitude over the range [-4, 2.4] of log(E_rel)
! in steps of 0.2. Use simple interpolation to get the rlative amplitude.

log_E = log10(E_rel)

if (log_E <= -4) then
  r = 0.25
elseif (log_E >= 2.4) then
  r = 0
else
  x = (log_E + 4) / 0.2
  i = int(x)
  frac = x - i
  r = (1 - frac) * rel_amp(i) + frac * rel_amp(i+1)
endif

c1 = 1 / (1 + r)
c2 = r / (1 + r)

! Now use this to generate a Gaussian

if (present(r_in)) then
  rr = r_in
else
  call ran_uniform(rr)
endif

if (rr > 0.5) then
  r = rr - 0.5
  ss = 1
else
  r = 0.5 - rr
  ss = -1
endif

gamma_phi = ss * inverse (g_tot, r, 0.0_rp, sigma_cut * max(s1, s2), 0.0001_rp)

end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------

function g_tot(x) result (g)

real(rp) x, g
real(rp) xx, r

integer i

!

xx = n_pts_per_sigma * x / s1
i = int(xx)
r = xx - i

if (i >= n_g) then
  g = c1 * g1(n_g)
else
  g = c1 * ((1 - r) * g1(i) + r * g1(i+1))
endif

xx = n_pts_per_sigma * x / s2
i = int(xx)
r = xx - i

if (i >= n_g) then
  g = g + c2 * g2(n_g)
else
  g = g + c2 * ((1 - r) * g2(i) + r * g2(i+1))
endif

end function

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!+
! Subroutine photon_energy_init (E_rel, r_in)
!
! Routine to convert a random number in the interval [0,1] to a photon energy.
! The photon probability spectrum is:
!   P(E_rel) = 0.1909859 * Integral_{E_rel}^{Infty} K_{5/3}(x) dx
! Where
!   E_rel = Relative photon energy: E / E_crit, E_crit = Critical energy.
!   K_{5/3} = Modified Bessel function.
! There is a cut-off built into the calculation so that E_rel will be in the 
! range [0, 16]
!
! Module needed:
!   use photon_init_mod
!
! Input:
!   r_in  -- Real(rp), optional: Number in the range [0,1].
!             If not present, a random number will be used.
!
! Output:
!   E_rel -- Real(rp): Relative photon energy E/E_crit. 
!-

subroutine photon_energy_init (E_rel, r_in)

implicit none

! An init_spectrum_struct holds a spline fit of E_rel vs integrated probability
! over a certain range. The spline is of the form:
!   e_rel = c0 + c1 * r + cn * r^n

type init_spectrum_struct
  real(rp) del_p
  real(rp) p_max
  real(rp), allocatable :: c0(:), c1(:), cn(:), n(:) ! Spline fit coefs.
end type

type (init_spectrum_struct), save :: spec(4)

real(rp) E_rel
real(rp), optional :: r_in
real(rp) a, b2, t, rr, r, r_rel, b, c, rr0

integer i, is

logical, save :: init_needed = .true.


! Note: all the arrays are paded with an extra value to prevent
! array out of bounds problems due to rounding.

if (present(r_in)) then
  rr = r_in
  if (rr < 0  .or. rr > 1) then
    print *, 'ERROR: RR IS OUT OF RANGE: ', rr
    stop
  endif
else
  call ran_uniform(rr)
endif

! Init

if (init_needed) then

  spec(:)%del_p = (/ 0.02, 0.01, 0.001, 0.0001 /)
  spec(:)%p_max = (/ 0.8, 0.99, 0.999, 0.9999 /)

  allocate (spec(1)%c0(0:41), spec(1)%c1(0:41), spec(1)%cn(0:41), spec(1)%n(0:41))
  allocate (spec(2)%c0(0:20), spec(2)%c1(0:20), spec(2)%cn(0:20), spec(2)%n(0:20))
  allocate (spec(3)%c0(0:10), spec(3)%c1(0:10), spec(3)%cn(0:10), spec(3)%n(0:10))
  allocate (spec(4)%c0(0:10), spec(4)%c1(0:10), spec(4)%cn(0:10), spec(4)%n(0:10))

  spec(1)%c0(:) = (/ 0.0, &
              4.28341e-6, 0.0000342902, 0.000115858, &
              0.000275057, 0.000538307, 0.0009325, 0.00148513, 0.00222445, &
              0.00317957, 0.00438067, 0.00585915, 0.00764782, 0.00978112, &
              0.0122954, 0.015229, 0.0186228, 0.0225205, 0.026969, 0.0320186, &
              0.037724, 0.0441446, 0.0513452, 0.0593972, 0.0683791, 0.0783781, &
              0.0894916, 0.101829, 0.115512, 0.130682, 0.147499, 0.166146, &
              0.186838, 0.209825, 0.235404, 0.263928, 0.295827, 0.331625, 0.371975, &
              0.417703, 0.469878, 0.0 /)

  spec(1)%c1(:) = (/ 0.0, & 
              0.000642606, 0.00257329, 0.00580068, 0.0103393, 0.0162097, 0.0234387, &
              0.0320599, 0.042114, 0.0536493, 0.0667223, 0.081399, 0.097755, &
              0.115877, 0.135866, 0.157835, 0.181913, 0.20825, 0.237015, 0.268401, &
              0.302633, 0.339965, 0.380692, 0.425158, 0.473758, 0.526956, 0.585298, &
              0.649425, 0.720102, 0.798242, 0.884947, 0.981556, 1.08971, 1.21146, &
              1.34935, 1.50665, 1.68757, 1.89762, 2.1442, 2.43746, 2.79164, 0.0 /)

  spec(2)%c0(:) = (/ &
              0.469878, 0.498807, 0.529911, 0.563449, 0.599724, 0.639103, 0.682028, &
              0.72904, 0.780816, 0.83821, 0.902325, 0.974622, 1.05709, 1.15252, 1.26503, &
              1.40105, 1.57139, 1.79651, 2.12271, 2.69945, 0.0 /)

  spec(2)%c1(:) = (/ & 
              2.79164, 2.9977, 3.22744, 3.48512, 3.77607, 4.10707, 4.48684, 4.92685, &
              5.44242, 6.05453, 6.79263, 7.69945, 8.83931, 10.3137, 12.2927, 15.0837, &
              19.3053, 26.4104, 40.7923, 84.6624, 0.0 /)

  spec(3)%c0(:) = (/ &
              2.69945, 2.78886, 2.88928, 3.0037, 3.13651, 3.29451, 3.48916, 3.742, &
              4.10159, 4.72378, 0.0 /)

  spec(3)%c1(:) = (/ & 
              84.6624, 94.5006, 106.83, 122.726, 143.986, 173.851, 218.823, 294.119, &
              445.575, 903.65, 0.0 /)

  spec(4)%c0(:) = (/ &
              4.72378, 4.81908, 4.92582, 5.04708, 5.18739, 5.35376, 5.55797, 5.82214, &
              6.19606, 6.83908, 0.0 /)

  spec(4)%c1(:) = (/ & 
              903.65, 1005.91, 1133.9, 1298.7, 1518.78, 1827.44, 2291.39, 3066.54, 4621.7, &
              9308.58, 0.0 /)

  ! Fill in rest of the spline fit coefs

  do is = 1, 4
    spec(is)%c1(:) = spec(is)%c1(:) * spec(is)%del_p
    do i = 0, ubound(spec(is)%c0, 1) - 1
      spec(is)%n(i) = (spec(is)%c1(i+1) - spec(is)%c1(i)) / &
                            (spec(is)%c0(i+1) -spec(is)%c0(i) - spec(is)%c1(i))
      spec(is)%cn(i) = spec(is)%c0(i+1) - spec(is)%c0(i) - spec(is)%c1(i)
    enddo
  enddo

  ! In the range above rr = 0.9999: This is such a small part of the spectrum a crude
  ! approximation is use here of the form a + b / (c-rr)
  ! Want to match to the spline fits at rr = 0.9999

  c = 1 + 1e-5
  rr0 = 0.9999
  b = 9308.58 * (c - rr0)**2
  a = 6.83908 - b / (c - rr0)

  init_needed = .false.
endif

! If in range of [0, 0.9999]: Use one of the four spline fits.

do is = 1, 4
  if (rr > spec(is)%p_max) cycle
  if (is == 1) then
    r_rel = rr / spec(is)%del_p
  else
    r_rel = (rr - spec(is-1)%p_max) / spec(is)%del_p
  endif
  i = int(r_rel)
  r = r_rel - i
  E_rel = spec(is)%c0(i) + spec(is)%c1(i) * r + spec(is)%cn(i) * r**spec(is)%n(i)
  return
enddo

! In the range above rr = 0.9999: This is such a small part of the spectrum a crude
! approximation is use here of the form a + b / (c-rr)

E_rel = a + b / (c-rr)

end subroutine

end module
