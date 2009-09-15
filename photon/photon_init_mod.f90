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

real(rp) E_rel
real(rp), optional :: r_in
real(rp) a, b2, t, rr, r, r_rel

real(rp), parameter :: e_rel_arr(0:91) = (/ &
              0.0, 5.35337e-7, 4.28341e-6, 0.0000144605, &
              0.0000342902, 0.0000670066, 0.000115858, 0.000184113, 0.000275057, &
              0.000392008, 0.000538307, 0.000717333, 0.0009325, 0.00118727, &
              0.00148513, 0.00182966, 0.00222445, 0.00267317, 0.00317957, &
              0.00374744, 0.00438067, 0.00508322, 0.00585915, 0.0067126, &
              0.00764782, 0.00866917, 0.00978112, 0.0109883, 0.0122954, 0.0137072, &
              0.015229, 0.0168657, 0.0186228, 0.0205058, 0.0225205, 0.0246729, &
              0.026969, 0.0294153, 0.0320186, 0.0347857, 0.037724, 0.040841, &
              0.0441446, 0.0476431, 0.0513452, 0.0552601, 0.0593972, 0.0637666, &
              0.0683791, 0.0732456, 0.0783781, 0.083789, 0.0894916, 0.0954998, &
              0.101829, 0.108494, 0.115512, 0.122902, 0.130682, 0.138874, 0.147499, &
              0.156581, 0.166146, 0.176222, 0.186838, 0.198027, 0.209825, 0.22227, &
              0.235404, 0.249273, 0.263928, 0.279426, 0.295827, 0.313201, 0.331625, &
              0.351184, 0.371975, 0.394107, 0.417703, 0.442906, 0.469878, 0.498807, &
              0.529911, 0.563449, 0.599724, 0.639103, 0.682028, 0.72904, 0.780816, &
              0.83821, 0.902325, 0.974622 /)


integer i

! e_rel_arr maps rr in the range [0, 0.90] to E_rel in steps of 0.01.
! de_rel_arr is the derivative of e_rel_arr
! Use a cubic spline fit.
! Note: all the arrays are paded with an extra value to prevent
! array out of bounds problems due to rounding.

if (present(r_in)) then
  rr = r_in
else
  call ran_uniform(rr)
endif

if (rr < 0  .or. rr > 1) then
  print *, 'ERROR: RR IS OUT OF RANGE: ', rr
  stop
endif

! If in range of [0, 0.90]:


end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!+
! Subroutine photon_energy_init2 (E_rel, r_in)
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

subroutine photon_energy_init2 (E_rel, r_in)

implicit none

real(rp) E_rel
real(rp), optional :: r_in
real(rp) a, b2, t, rr, r, r_rel

real(rp), parameter :: f2(0:100) = (/ &
              0.0, 5.35337e-7, 4.28341e-6, 0.0000144605, &
              0.0000342902, 0.0000670066, 0.000115858, 0.000184113, 0.000275057, &
              0.000392008, 0.000538307, 0.000717333, 0.0009325, 0.00118727, &
              0.00148513, 0.00182966, 0.00222445, 0.00267317, 0.00317957, &
              0.00374744, 0.00438067, 0.00508322, 0.00585915, 0.0067126, &
              0.00764782, 0.00866917, 0.00978112, 0.0109883, 0.0122954, 0.0137072, &
              0.015229, 0.0168657, 0.0186228, 0.0205058, 0.0225205, 0.0246729, &
              0.026969, 0.0294153, 0.0320186, 0.0347857, 0.037724, 0.040841, &
              0.0441446, 0.0476431, 0.0513452, 0.0552601, 0.0593972, 0.0637666, &
              0.0683791, 0.0732456, 0.0783781, 0.083789, 0.0894916, 0.0954998, &
              0.101829, 0.108494, 0.115512, 0.122902, 0.130682, 0.138874, 0.147499, &
              0.156581, 0.166146, 0.176222, 0.186838, 0.198027, 0.209825, 0.22227, &
              0.235404, 0.249273, 0.263928, 0.279426, 0.295827, 0.313201, 0.331625, &
              0.351184, 0.371975, 0.394107, 0.417703, 0.442906, 0.469878, 0.498807, &
              0.529911, 0.563449, 0.599724, 0.639103, 0.682028, 0.72904, 0.780816, &
              0.83821, 0.902325, 0.974622, 1.05709, 1.15252, 1.26503, 1.40105, &
              1.57139, 1.79651, 2.12271, 2.69945, 0.0 /)

real(rp), parameter :: f3(0:10) = (/ 2.69945, 2.78886, 2.88928, 3.0037, &
                    3.13651, 3.29451, 3.48916, 3.742, 4.10159, 4.72378, 0.0 /)

real(rp), parameter :: f4(0:10) = (/ 4.72378, 4.81908, 4.92582, 5.04708, &
                    5.18739, 5.35376, 5.55797, 5.82214, 6.19606, 6.83908, 0.0 /)

real(rp), parameter :: f5(0:10) = (/ 6.83908, 6.93721, 7.04701, 7.17163, &
                    7.31567, 7.48627, 7.6954, 7.96553, 8.34717, 9.00187, 0.0 /)

real(rp), parameter :: f6(0:10) = (/ 9.00187, 9.10162, 9.2132, 9.33978, &
                    9.48601, 9.65912, 9.87119, 10.1449, 10.5313, 11.1934, 0.0 /)

real(rp), parameter :: f7(0:11) = (/ 11.1934, 11.2942, 11.4069, 11.5347, &
                    11.6823, 11.8571, 18.7481, 12.3472, 12.7366, 13.4036, 16.0, 0.0 /)

integer i

! f2 maps rr in the range [0, 0.99] to E_rel in steps of 0.01.
! f3 maps rr in the range [0.99, 0.999] to E_rel in steps of 0.001.
! f4 maps rr in the range [0.999, 0.9999] to E_rel in steps of 0.0001.
! f5 maps rr in the range [0.9999, 0.99999] to E_rel in steps of 0.00001.
! f6 maps rr in the range [0.99999, 0.999999] to E_rel in steps of 1E-6.
! f7 maps rr in the range [0.999999, 1.0] to E_rel in steps of 1E-7.
! Things are smooth enough to use a simple linear interpolation.
! Note: all the arrays are paded with an extra value to prevent
! array out of bounds problems due to rounding.
! The very last bin is a fudge but the width of this bin is 1e-7.

if (present(r_in)) then
  rr = r_in
else
  call ran_uniform(rr)
endif

if (rr < 0  .or. rr > 1) then
  print *, 'ERROR: RR IS OUT OF RANGE: ', rr
  stop
endif

! If in range of f2:

r_rel = 100 * rr
if (r_rel <= 99) then 
  i = int(r_rel)
  r = r_rel - i
  E_rel = f2(i) * (1 - r) + f2(i+1) * r
  return
endif

! If in f3:

r_rel = 10 * r_rel - 990
if (r_rel <= 9) then
  i = int(r_rel)
  if (i < 0) i = 0  ! to prevent array out of bounds problems
  r = r_rel - i
  E_rel = f3(i) * (1 - r) + f3(i+1) * r
  return
endif
 
! If in f4:

r_rel = 10 * r_rel - 90
if (r_rel <= 9) then
  i = int(r_rel)
  if (i < 0) i = 0  ! to prevent array out of bounds problems
  r = r_rel - i
  E_rel = f4(i) * (1 - r) + f4(i+1) * r
  return
endif
 
! If in f5:

r_rel = 10 * r_rel - 90
if (r_rel <= 9) then
  i = int(r_rel)
  if (i < 0) i = 0  ! to prevent array out of bounds problems
  r = r_rel - i
  E_rel = f5(i) * (1 - r) + f5(i+1) * r
  return
endif
 
! If in f6:

r_rel = 10 * r_rel - 90
if (r_rel <= 9) then
  i = int(r_rel)
  if (i < 0) i = 0  ! to prevent array out of bounds problems
  r = r_rel - i
  E_rel = f6(i) * (1 - r) + f6(i+1) * r
  return
endif
 
! Must be in f7:

r_rel = 10 * r_rel - 90
i = int(r_rel)
if (i < 0) i = 0  ! to prevent array out of bounds problems
r = r_rel - i
E_rel = f7(i) * (1 - r) + f7(i+1) * r
 
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------


end module
