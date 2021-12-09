module photon_init_mod

use coord_mod

! An init_spectrum_struct holds an array of spline fits of E_rel of gamma_phi vs 
! integrated probability r over a certain range. Each spline section is fit to 
! one of two forms:
! 
! The gen_poly_spline$ type is of the form:
! [See the code how which form is chosen.]
!   fit_value = c0 + c1 * t + c2 * t^c3      
! or
!   fit_value = c0 + c1 * t + c2 * t^2 + c3 * t^3
!
! The end_spline$ type is of the form:
!   fit_value = c0 + c1 * t + c2 * t^2 / (1 - t/c3)
!
! where t is in the range [0, 1].
!   t = 0 at the start of the spline section.
!   t = 1 at the end of the spline section.

type photon_init_spline_pt_struct
  real(rp) c0, c1, c2, c3
end type

type photon_init_spline_struct
  real(rp) x_min    ! Lower bound
  real(rp) x_max    ! Upper bound of Region of validity of this spline fit.
                    ! The lower bound is given by the upper bound of the previos struct.
  real(rp) del_x    ! Spacing between spline points
  type (photon_init_spline_pt_struct), allocatable :: pt(:)
  integer spline_type
end type

integer, parameter :: gen_poly_spline$ = 1, end_spline$ = 2

private photon_init_spline_eval, photon_init_spline_coef_calc

contains

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!+
! Subroutine absolute_photon_position (e_orb, photon_orb)
! 
! Routine to calculate the photon phase space coordinates given:
!   1) The phase space coords of the emitting charged particle and
!   2) The photon phase space coords relative to the emitting particle.
!      The photon (x, y, z) position is ignored (it is assumed the photon is emitted at
!      the charged particle position) and only the photon's (vx, vy, vz) velocity matters.
!
! Input:
!   e_orb      -- coord_struct: charged particle position.
!   photon_orb -- coord_struct: Photon position relative to e_orb.
!
! Output:
!   photon_orb -- coord_struct: Absolute photon position.
!-

subroutine absolute_photon_position (e_orb, photon_orb)

implicit none

type (coord_struct) photon_orb, e_orb
real(rp) e_vec(3), w_mat(3,3), theta

! Remember: Phase space description for charged particle is different from photons.

photon_orb%vec(1) = e_orb%vec(1)
photon_orb%vec(3) = e_orb%vec(3)

e_vec(1:2) = e_orb%vec(2:4:2) / (e_orb%p0c * (1 + e_orb%vec(6)))
theta = asin(norm2(e_vec(1:2)))
if (theta == 0) return
call axis_angle_to_w_mat ([-e_vec(2), e_vec(1), 0.0_rp], theta, w_mat)
photon_orb%vec(2:6:2) = matmul(w_mat, photon_orb%vec(2:6:2))

end subroutine absolute_photon_position

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!+
! Subroutine bend_photon_init (g_bend_x, g_bend_y, gamma, orbit, E_min, E_max, E_integ_prob,
!                                         vert_angle_min, vert_angle_max, vert_angle_symmetric, emit_probability)
!
! Routine to initalize a photon for dipole bends and wigglers (but not undulators).
! The photon is initialized using the standard formulas for bending radiation.
!
! The energy of the photon is calculated in one of two ways:
!
!   1) If E_integ_prob is present and non-negative, the photon energy E will be such that the integrated 
!       probability  [E_min, E] relative to the integrated probability in the range [E_min, E_max] is E_integ_prob. 
!       That is, E_integ_prob can be used to to give a set of photon energies equally spaced in terms of the 
!       integrated probability distribution.
!
!   2) If E_integ_prob is not present, or is negative, the photon energy is chosen at random in 
!       the range [E_min, E_max].
!
! An E_integ_prob of zero means that the generated photon will have energy E_min.
! An E_integ_prob of one means that the generated photon will have energy E_max.
!
! The photon's polarization, will have unit amplitude.
!
! This routine assumes that the emitting charged particle is on-axis and moving in 
! the forward direction. To correct for the actual charged particle postion use the routine
!   absolute_photon_position
!
! Input:
!   g_bend_x             -- real(rp): Bending 1/rho component in horizontal plane.
!   g_bend_y             -- real(rp): Bending 1/rho component in vertical plane.
!   gamma                -- real(rp): Relativistic gamma factor of generating charged particle.
!   E_min                -- real(rp), optional: Minimum photon energy. Default is zero. Ignored if negative.
!   E_max                -- real(rp), optional: Maximum photon energy.  Default is Infinity. Ignored if negative.
!                            If non-positive then E_max will be taken to be Infinity.
!   E_integ_prob         -- real(rp):, optional :: integrated energy probability. See above.
!                            If E_integ_prob is non-negative, it must be in the range [0, 1].
!   vert_angle_min       -- real(rp), optional: Minimum vertical angle to emit a photon. 
!                           -pi/2 is used if argument not present or if argument is less than -pi/2.
!   vert_angle_max       -- real(rp), optional: Maximum vertical angle to emit a photon. 
!                           pi/2 is used if argument not present or if argument is greater than pi/2.
!   vert_angle_symmetric -- logical, optional: Default is False. If True, photons will be emitted
!                             in the range [-vert_angle_max, -vert_angle_min] as well as the range
!                             [vert_angle_min, vert_angle_max]. In this case vert_angle_min/max must be positive.
!   emit_probability     -- real(rp), optional: Probability of emitting a photon in the range [E_min, E_max] or 
!                             in the vertical angular range given. The probability is normalized so that the 
!                             probability of emitting if no ranges are given is 1.
!
! Output:
!   orbit            -- coord_struct: Initialized photon.
!-

subroutine bend_photon_init (g_bend_x, g_bend_y, gamma, orbit, E_min, E_max, E_integ_prob, &
                                         vert_angle_min, vert_angle_max, vert_angle_symmetric, emit_probability)


implicit none

type (coord_struct) orbit
real(rp), optional :: E_min, E_max, E_integ_prob, emit_probability, vert_angle_min, vert_angle_max
real(rp) g_bend_x, g_bend_y, g_bend, gamma, phi
real(rp) E_rel, E_photon, r_min, r_max, r, f, phi_min, phi_max, integ_prob
integer sgn
logical, optional :: vert_angle_symmetric

! Photon energy

g_bend = sqrt(g_bend_x**2 + g_bend_y**2)
phi_min = real_option(-pi/2, vert_angle_min)
phi_max = real_option(pi/2, vert_angle_max)

r = real_option(-1.0_rp, E_integ_prob)
if (r < 0) call ran_uniform(r)

if (phi_min > -pi/2 .or. phi_max < pi/2) then
  integ_prob = init_photon_integ_prob (gamma, g_bend, E_min, E_max, vert_angle_min, &
                                         vert_angle_max, vert_angle_symmetric, r, E_photon)
  E_rel = E_photon / E_crit_photon(gamma, g_bend)
  if (present(emit_probability)) emit_probability = integ_prob

else
  r_min = 0
  r_max = 1

  if (real_option(0.0_rp, E_min) > 0) r_min = bend_photon_energy_integ_prob(E_min, g_bend, gamma)
  if (real_option(0.0_rp, E_max) > 0) r_max = bend_photon_energy_integ_prob(E_max, g_bend, gamma)
  if (present(emit_probability)) emit_probability = r_max - r_min

  r = r_min + r * (r_max - r_min)
  E_rel = bend_photon_e_rel_init (r)

  E_photon = E_rel * E_crit_photon(gamma, g_bend)
  if (present(emit_probability)) emit_probability = r_max - r_min
endif

!------------
! Photon vertical angle

r_min = bend_vert_angle_integ_prob(phi_min, E_rel, gamma)
r_max = bend_vert_angle_integ_prob(phi_max, E_rel, gamma)

call ran_uniform(r)
if (logic_option(.false., vert_angle_symmetric)) then
  f = 2
  if (r > 0.5_rp) then
    sgn = 1
    r = 2 * (r - 0.5_rp)
  else
    sgn = -1
    r = 2 * r
  endif
else
  f = 1
  sgn = 1
endif

r = r_min + r * (r_max - r_min)
phi = sgn * bend_photon_vert_angle_init (E_rel, gamma, r)

orbit%vec = 0
orbit%vec(2) =  g_bend_y * sin(phi) / g_bend
orbit%vec(4) = -g_bend_x * sin(phi) / g_bend
orbit%vec(6) = cos(phi)

call init_coord (orbit, orbit%vec, particle = photon$, E_photon = E_photon)

! Polaraization

call bend_photon_polarization_init(g_bend_x, g_bend_y, E_rel, gamma*phi, orbit)

end subroutine bend_photon_init

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!+
! Function bend_photon_energy_integ_prob (E_photon, g_bend, gamma) result (integ_prob)
!
! Routine to find the integrated probability corresponding to emitting a photon
! from a bend in the range [0, E_photon].
!
! Input:
!   E_photon -- real(rp): Photon energy.
!   g_bend   -- real(rp): 1/rho bending strength.
!   gamma    -- Real(rp): Relativistic gamma factor of generating charged particle.
!
! Output:
!   integ_prob -- real(rp): Integrated probability. Will be in the range [0, 1].
!-

function bend_photon_energy_integ_prob (E_photon, g_bend, gamma) result (integ_prob)

use super_recipes_mod, only: super_zbrent

implicit none

real(rp) E_photon, g_bend, gamma, integ_prob, E1, E_rel_target
integer status

! Easy cases. photon_e_rel_init gives a finite energy at integ_prob = 1 (in theory 
! should be infinity) so return 1.0 if E_photon > E (upper bound).

if (E_photon == 0) then
  integ_prob = 0
  return
endif

E_rel_target = E_photon / E_crit_photon(gamma, g_bend)

E1 = bend_photon_e_rel_init (1.0_rp)
if (E_rel_target >= E1) then
  integ_prob = 1
  return
endif

! bend_photon_e_rel_init calculates photon energy given the integrated probability
! so invert using the NR routine zbrent.

integ_prob = super_zbrent(energy_func, 0.0_rp, 1.0_rp, 1e-15_rp, 1d-12, status)

!----------------------------------------------------------------------------------------
contains

function energy_func(integ_prob, status) result (dE)

real(rp), intent(in) :: integ_prob
real(rp) dE, E_rel
integer status

!

E_rel = bend_photon_e_rel_init(integ_prob)
dE = E_rel - E_rel_target

end function energy_func

end function bend_photon_energy_integ_prob

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!+
! Function bend_vert_angle_integ_prob (vert_angle, E_rel, gamma) result (integ_prob)
!
! Routine to find the integrated probability corresponding to emitting a photon
! from a bend and with relative energy E_rel in the vertical angle range [-pi/2, vert_angle/2].
!
! Note: vert_angle is allowed to be out of the range [-pi/2, pi/2]. In this case, integ_prob
! will be set to 0 or 1 as appropriate.
!
! Input:
!   vert_angle  -- real(rp): Vertical angle.
!   E_rel       -- real(rp): Relative photon energy E/E_crit. 
!   gamma       -- real(rp): Relativistic gamma factor of generating charged particle.
!
! Output:
!   integ_prob  -- real(rp): Integrated probability. Will be in the range [0, 1].
!-

function bend_vert_angle_integ_prob (vert_angle, E_rel, gamma) result (integ_prob)

use super_recipes_mod

implicit none

real(rp) vert_angle, E_rel, gamma, integ_prob
integer status
logical invert

! If angle is so large that bend_photon_vert_angle_init is inaccurate, just round off to 0 or 1.

if (abs(bend_photon_vert_angle_init(E_rel, gamma, 1.0e-30_rp)) <= abs(vert_angle)) then
  if (vert_angle < 0) then
    integ_prob = 0
  else
    integ_prob = 1
  endif
  return
endif

! Invert using the NR routine zbrent.
! Use the inverted probability in the calculation for positive angles since it is more accurate.

invert = (vert_angle > 0)
integ_prob = super_zbrent(vert_angle_func, 0.0_rp, 1.0_rp, 1d-12, 1d-20, status)
if (invert) integ_prob = 1.0_rp - integ_prob

!----------------------------------------------------------------------------------------
contains

function vert_angle_func(integ_prob, status) result (d_angle)

real(rp), intent(in) :: integ_prob
real(rp) angle, d_angle
integer status

!

angle = bend_photon_vert_angle_init(E_rel, gamma, integ_prob, invert)
d_angle = angle - vert_angle

end function vert_angle_func

end function bend_vert_angle_integ_prob

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!+
! Subroutine bend_photon_polarization_init (g_bend_x, g_bend_y, E_rel, gamma_phi, orbit)
!
! Routine to set a photon's polarization.
! The photon's polarization will be either in the plane of the bend or out of the plane and
! the magnitude will be 1.
!
! Input:
!   g_bend_x  -- Real(rp): Bending 1/rho component in horizontal plane.
!   g_bend_y  -- Real(rp): Bending 1/rho component in vertical plane.
!   E_rel     -- Real(rp): Relative photon energy E/E_crit. 
!   gamma_phi -- Real(rp): gamma * phi where gamma is the beam relativistic factor and
!                  phi is the vertical photon angle (in radians).
! 
! Output:
!   orbit        -- coord_struct: Photon coords
!     %field(2)     -- (x,y) polaraization. Will have unit magnitude
!     %phase(2)     -- (x,y) phases. Will be [0, pi/2].
!-

subroutine bend_photon_polarization_init (g_bend_x, g_bend_y, E_rel, gamma_phi, orbit)

use fgsl, only: fgsl_sf_bessel_Knu

implicit none

type (coord_struct) :: orbit
real(rp) g_bend_x, g_bend_y, g_bend, gamma_phi
real(rp) gp2, xi, k_23, k_13, pol_x, pol_y
real(rp) E_rel, E_photon

!

g_bend = sqrt(g_bend_x**2 + g_bend_y**2)
gp2 = (gamma_phi)**2
xi = E_rel * sqrt(1+gp2)**3 / 2

if (xi > 100) then
  k_13 = 1 
  k_23 = 1  ! Ratio is 1.
else
  k_13 = fgsl_sf_bessel_Knu(1.0_rp/3, xi)
  k_23 = fgsl_sf_bessel_Knu(2.0_rp/3, xi)
endif

pol_x = k_23
pol_y = k_13 * sqrt(gp2 / (1 + gp2))

orbit%field = [pol_x, sign(pol_y, gamma_phi)] / sqrt(pol_x**2 + pol_y**2)
orbit%phase = [0.0_rp, pi/2]

end subroutine bend_photon_polarization_init

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!+
! Function bend_photon_vert_angle_init (E_rel, gamma, r_in, invert) result (phi)
!
! Routine to convert an integrated probability to a vertical angle for emitting a photon from a bend. 
! The integrated probability is in the range [0,1] with 0 corresponding to a phi = -pi/2 and 
! integrated probability of 1 corresponding to phi = pi/2.
!
! Input:
!   E_rel   -- real(rp): Relative photon energy E/E_crit. 
!		gamma   -- real(rp): beam relativistic factor 
!   r_in    -- real(rp), optional: Integrated probability in the range [0,1].
!               If not present, a random number will be used.
!   invert  -- real(rp), optional: If True then take r_in as the inverse integrated probability 
!               with inverted probability = 1 - probability. This is useful to avoid round-off
!               errors when for looking at the tail of the distribution where the integrated prob is 
!               very close to 1 and small deviations can have large effects. Default is False.
! 
! Output:
!   phi   -- real(rp): The photon vertical emission angle (in radians).
!                  Note: phi is an increasing monotonic function of r_in.
!-

function bend_photon_vert_angle_init (E_rel, gamma, r_in, invert) result (phi)

use super_recipes_mod

implicit none

type photon_vert_angle_pt_struct
  real(rp) gp           ! Gamma-phi
  real(rp) dgp_dprob    ! dGamma-phi/dIntegratedProbability
  real(rp) dgp_dlogE    ! dGamm-phi/dlog_E
  real(rp) dgp_d2       ! d^2Gamma-phi/(dlog_E dIntegratedProbability)
end type

integer, parameter :: n_prob_pt = 63, n_E_pt = 7, n_width = 4

type (photon_init_spline_struct), save :: angle_half_width(n_width)
type (photon_vert_angle_pt_struct), save :: gp(n_E_pt, 0:n_prob_pt)

real(rp), optional :: r_in
real(rp), save :: prob_pt(0:n_prob_pt), inv_prob_pt(0:n_prob_pt), log_E_pt(n_E_pt)
real(rp) rr, inv_rr, phi, E_rel, gamma, gamma_phi, log_E_rel, log_E_prob, dval_dE, dval_dprob
real(rp) f, sig

integer i, sign_phi, ix_E, ix_prob
logical, optional :: invert
logical, save :: init_needed = .true.

! Init

if (init_needed) then
  call init_this()
  init_needed = .false.
endif

!----------------------------------------------------
! Integrated probability number to use.

if (present(r_in)) then
  rr = r_in
else
  call ran_uniform(rr)
endif

! The fit is only for positive phi.
! So make phi negative if rr < 0.5

if (logic_option(.false., invert)) then
  inv_rr = rr
  rr = 1.0_rp - inv_rr
else
  inv_rr = 1.0_rp - rr
endif

if (rr > 0.5_rp) then
  inv_rr = 2.0_rp * inv_rr
  rr = 2.0_rp * rr - 1.0_rp
  sign_phi = 1
else
  inv_rr = 2.0_rp * rr
  rr = 1.0_rp - 2.0_rp * rr
  sign_phi = -1
endif

! Border cases:

if (E_rel == 0 .or. inv_rr == 0) then
  phi = sign_phi * pi / 2
  return
endif

! Find energy range.
! log_E_prob is the value of log10(E_rel) to use with the probability table.

log_E_rel = log10(E_rel)

if (log_E_rel < log_E_pt(1)) then
  log_E_prob = log_E_pt(1)
  ix_E = 1
elseif (log_E_rel > log_E_pt(n_E_pt)) then
  log_E_prob = log_E_pt(n_E_pt)
  ix_E = n_E_pt - 1
else
  log_E_prob = log_E_rel
  do ix_E = 1, n_E_pt - 1
    if (log_E_prob <= log_E_pt(ix_E+1)) exit
  enddo
endif

! Find probability integral range

do ix_prob = 0, n_prob_pt - 1
  if (rr <= prob_pt(ix_prob+1)) exit
enddo

if (ix_prob == n_prob_pt) then  ! Extrapolate
  call super_bicubic_interpolation( &
          [gp(ix_E,ix_prob-1)%gp,        gp(ix_E+1,ix_prob-1)%gp,        gp(ix_E+1,ix_prob)%gp,        gp(ix_E,ix_prob)%gp], &
          [gp(ix_E,ix_prob-1)%dgp_dlogE, gp(ix_E+1,ix_prob-1)%dgp_dlogE, gp(ix_E+1,ix_prob)%dgp_dlogE, gp(ix_E,ix_prob)%dgp_dlogE], &
          [gp(ix_E,ix_prob-1)%dgp_dprob, gp(ix_E+1,ix_prob-1)%dgp_dprob, gp(ix_E+1,ix_prob)%dgp_dprob, gp(ix_E,ix_prob)%dgp_dprob], &
          [gp(ix_E,ix_prob-1)%dgp_d2,    gp(ix_E+1,ix_prob-1)%dgp_d2,    gp(ix_E+1,ix_prob)%dgp_d2,    gp(ix_E,ix_prob)%dgp_d2], &
          log_E_pt(ix_E), log_E_pt(ix_E+1), prob_pt(ix_prob-1), prob_pt(ix_prob), log_E_prob, prob_pt(n_prob_pt), gamma_phi, dval_dE, dval_dprob)
  gamma_phi = gamma_phi * log((1-prob_pt(n_prob_pt))/inv_rr)

else ! Interpolate
  call super_bicubic_interpolation( &
          [gp(ix_E,ix_prob)%gp,        gp(ix_E+1,ix_prob)%gp,        gp(ix_E+1,ix_prob+1)%gp,        gp(ix_E,ix_prob+1)%gp], &
          [gp(ix_E,ix_prob)%dgp_dlogE, gp(ix_E+1,ix_prob)%dgp_dlogE, gp(ix_E+1,ix_prob+1)%dgp_dlogE, gp(ix_E,ix_prob+1)%dgp_dlogE], &
          [gp(ix_E,ix_prob)%dgp_dprob, gp(ix_E+1,ix_prob)%dgp_dprob, gp(ix_E+1,ix_prob+1)%dgp_dprob, gp(ix_E,ix_prob+1)%dgp_dprob], &
          [gp(ix_E,ix_prob)%dgp_d2,    gp(ix_E+1,ix_prob)%dgp_d2,    gp(ix_E+1,ix_prob+1)%dgp_d2,    gp(ix_E,ix_prob+1)%dgp_d2], &
          log_E_pt(ix_E), log_E_pt(ix_E+1), 0.0_rp, prob_pt(ix_prob+1) - prob_pt(ix_prob), log_E_prob, inv_prob_pt(ix_prob) - inv_rr, gamma_phi, dval_dE, dval_dprob)
endif

! Scale result by the sigma of the spectrum at fixed E_rel

if (log_E_rel < angle_half_width(1)%x_min) then
  f = photon_init_spline_eval(angle_half_width(1), log_E_rel, .true.)
elseif (log_E_rel > angle_half_width(n_width)%x_max) then
  f = photon_init_spline_eval(angle_half_width(n_width), log_E_rel, .true.)
else
  do i = 1, n_width - 1
    if (log_E_rel <= angle_half_width(i)%x_max) exit
  enddo
  f = photon_init_spline_eval(angle_half_width(i), log_E_rel)
endif

sig = 10.0_rp**min(100.0_rp, f)
phi = gamma_phi * sig * sign_phi / gamma

if (phi > pi/2) phi = pi/2
if (phi < -pi/2) phi = -pi/2

!------------------------------------------------------------------
contains

subroutine init_this()

real(rp), parameter :: loge10 = log(10.0_rp)
integer i, ai

!

angle_half_width%x_min = [-11.0_rp, -1.0_rp, 0.0_rp, 1.0_rp] ! Log10(E_rel) 
angle_half_width%x_max = [-1.0_rp,  0.0_rp, 1.0_rp, 2.0_rp]  ! Log10(E_rel) 
angle_half_width%del_x = [ 1.0_rp,  0.1_rp, 0.1_rp, 0.1_rp]
angle_half_width%spline_type = gen_poly_spline$

allocate (angle_half_width(1)%pt(0:10))
allocate (angle_half_width(2)%pt(0:10))
allocate (angle_half_width(3)%pt(0:10))
allocate (angle_half_width(4)%pt(0:10))

angle_half_width(1)%pt%c0 = log10([6.2172654D3, 2.8857985D3, 1.3394681D3, 6.2172412D2, 2.8857466D2, 1.3393563D2, 6.2148318D1, 2.8805653D1, 1.3282906D1, 5.9856983D0, 2.4632095D0])
angle_half_width(1)%pt%c1 = [-3.3333336D-1, -3.3333345D-1, -3.3333390D-1, -3.3333594D-1, -3.3334545D-1, -3.3338958D-1, -3.3359428D-1, -3.3454172D-1, -3.3888236D-1, -3.5785640D-1, -4.2524437D-1]

angle_half_width(2)%pt%c0 = log10([2.4632095D0, 2.2307087D0, 2.0150001D0, 1.8152935D0, 1.6308931D0, 1.4611731D0, 1.3055514D0, 1.1634623D0, 1.0343307D0, 9.1755087D-1, 8.1247064D-1])
angle_half_width(2)%pt%c1 = [-4.2524437D-1, -4.3603343D-1, -4.4740856D-1, -4.5921121D-1, -4.7123219D-1, -4.8321372D-1, -4.9485728D-1, -5.0583826D-1, -5.1582748D-1, -5.2451807D-1, -5.3165511D-1]

angle_half_width(3)%pt%c0 = log10([8.1247064D-1, 7.1838437D-1, 6.3453451D-1, 5.6012155D-1, 4.9432034D-1, 4.3630039D-1, 3.8524693D-1, 3.4038027D-1, 3.0097127D-1, 2.6635190D-1, 2.3592088D-1])
angle_half_width(3)%pt%c1 = [-5.3165511D-1, -5.3706333D-1, -5.4066816D-1, -5.4250530D-1, -5.4271650D-1, -5.4153128D-1, -5.3923832D-1, -5.3615168D-1, -5.3257821D-1, -5.2879133D-1, -5.2501420D-1]

angle_half_width(4)%pt%c0 = log10([2.3592088D-1, 2.0914478D-1, 1.8555567D-1, 1.6474640D-1, 1.4636457D-1, 1.3010598D-1, 1.1570824D-1, 1.0294479D-1, 9.1619591D-2, 8.1562582D-2, 7.2625755D-2])
angle_half_width(4)%pt%c1 = [-5.2501420D-1, -5.2141280D-1, -5.1809735D-1, -5.1512946D-1, -5.1253227D-1, -5.1030104D-1, -5.0841283D-1, -5.0683440D-1, -5.0552802D-1, -5.0445544D-1, -5.0358076D-1]

do i = 1, n_width
  call photon_init_spline_coef_calc (angle_half_width(i))
enddo

prob_pt(0) = 0
do i = 1, 9
  ai = real(i, rp)
  prob_pt(i)    =               ai*1.0e-1_rp
  prob_pt(i+ 9) = 0.9_rp      + ai*1.0e-2_rp
  prob_pt(i+18) = 0.99_rp     + ai*1.0e-3_rp
  prob_pt(i+27) = 0.999_rp    + ai*1.0e-4_rp
  prob_pt(i+36) = 0.9999_rp   + ai*1.0e-5_rp
  prob_pt(i+45) = 0.99999_rp  + ai*1.0e-6_rp
  prob_pt(i+54) = 0.999999_rp + ai*1.0e-7_rp
  inv_prob_pt(i)    = 1.0_rp    - ai*1.0e-1_rp
  inv_prob_pt(i+ 9) = 1.0e-1_rp - ai*1.0e-2_rp
  inv_prob_pt(i+18) = 1.0e-2_rp - ai*1.0e-3_rp
  inv_prob_pt(i+27) = 1.0e-3_rp - ai*1.0e-4_rp
  inv_prob_pt(i+36) = 1.0e-4_rp - ai*1.0e-5_rp
  inv_prob_pt(i+45) = 1.0e-5_rp - ai*1.0e-6_rp
  inv_prob_pt(i+54) = 1.0e-6_rp - ai*1.0e-7_rp
enddo

log_E_pt = [-6.0_rp, -3.0_rp, -2.0_rp, -1.0_rp, 0.0_rp, 1.0_rp, 2.0_rp] 

gp(1,  0: 9)%gp         = [0.0D0, 9.9113123D-2, 1.9373264D-1, 2.8223920D-1, 3.6566518D-1, 4.4635979D-1, 5.2742198D-1, 6.1316273D-1, 7.1128116D-1, 8.4271523D-1]
gp(1,  0: 9)%dgp_dprob  = [9.9999907D-1, 9.7394528D-1, 9.1531851D-1, 8.5622304D-1, 8.1590572D-1, 8.0295804D-1, 8.2541869D-1, 9.0197548D-1, 1.0908942D0, 1.6725271D0]
gp(1,  0: 9)%dgp_dlogE  = [9.7833318D-22, 2.3811340D-7, 1.5793575D-6, 4.1753932D-6, 7.6603029D-6, 1.1720158D-5, 1.6240675D-5, 2.1310736D-5, 2.7308501D-5, 3.5480781D-5] * loge10
gp(1,  0: 9)%dgp_d2     = [2.5597753D-10, 6.8724385D-6, 2.0198897D-5, 3.1119649D-5, 3.8112491D-5, 4.2947859D-5, 4.7650016D-5, 5.4391446D-5, 6.7384410D-5, 1.0439830D-4] * loge10

gp(1, 10:18)%gp         = [8.5998548D-1, 8.7863654D-1, 8.9901305D-1, 9.2160980D-1, 9.4717775D-1, 9.7694450D-1, 1.0131451D0, 1.0606159D0, 1.1339932D0]
gp(1, 10:18)%dgp_dprob  = [1.7919194D0, 1.9454862D0, 2.1398573D0, 2.3942485D0, 2.7425730D0, 3.2511291D0, 4.0704525D0, 5.6367994D0, 1.0024797D1]
gp(1, 10:18)%dgp_dlogE  = [3.6558936D-5, 3.7723626D-5, 3.8996194D-5, 4.0407308D-5, 4.2003486D-5, 4.3860722D-5, 4.6117255D-5, 4.9071814D-5, 5.3626793D-5] * loge10
gp(1, 10:18)%dgp_d2     = [1.1188707D-4, 1.2149817D-4, 1.3363859D-4, 1.4949867D-4, 1.7117596D-4, 2.0276891D-4, 2.5357167D-4, 3.5048759D-4, 6.2120197D-4] * loge10

gp(1, 19:27)%gp         = [1.1444241D0, 1.1558851D0, 1.1686344D0, 1.1830451D0, 1.1996870D0, 1.2194978D0, 1.2441987D0, 1.2775507D0, 1.3310960D0]
gp(1, 19:27)%dgp_dprob  = [1.0915040D1, 1.2060789D1, 1.3513344D1, 1.5419762D1, 1.8040994D1, 2.1890400D1, 2.8142387D1, 4.0232913D1, 7.4729365D1]
gp(1, 19:27)%dgp_dlogE  = [5.4273006D-5, 5.4982648D-5, 5.5771577D-5, 5.6662698D-5, 5.7690967D-5, 5.8913876D-5, 6.0436862D-5, 6.2490108D-5, 6.5778960D-5] * loge10
gp(1, 19:27)%dgp_d2     = [6.7600551D-4, 7.4655705D-4, 8.3592107D-4, 9.5318926D-4, 1.1141981D-3, 1.3505970D-3, 1.7339916D-3, 2.4746247D-3, 4.5835695D-3] * loge10

gp(1, 28:36)%gp         = [1.3388916D0, 1.3475072D0, 1.3571519D0, 1.3681289D0, 1.3809027D0, 1.3962410D0, 1.4155598D0, 1.4419731D0, 1.4851198D0]
gp(1, 28:36)%dgp_dprob  = [8.1805567D1, 9.0941620D1, 1.0256550D2, 1.1788229D2, 1.3903775D2, 1.7026987D2, 2.2131945D2, 3.2085503D2, 6.0831404D2]
gp(1, 28:36)%dgp_dlogE  = [6.6257008D-5, 6.6785113D-5, 6.7376020D-5, 6.8048189D-5, 6.8829903D-5, 6.9767873D-5, 7.0948223D-5, 7.2560156D-5, 7.5188689D-5] * loge10
gp(1, 28:36)%dgp_d2     = [5.0147235D-3, 5.5724146D-3, 6.2840525D-3, 7.2159074D-3, 8.5040747D-3, 1.0406571D-2, 1.3516264D-2, 1.9565960D-2, 3.7021484D-2] * loge10

gp(1, 37:45)%gp         = [1.4914736D0, 1.4985163D0, 1.5064254D0, 1.5154587D0, 1.5260125D0, 1.5387431D0, 1.5548653D0, 1.5770607D0, 1.6136777D0]
gp(1, 37:45)%dgp_dprob  = [6.6769791D2, 7.4452318D2, 8.4249062D2, 9.7190962D2, 1.1511776D3, 1.4167230D3, 1.8525304D3, 2.7067397D3, 5.1931886D3]
gp(1, 37:45)%dgp_dlogE  = [7.5575298D-5, 7.6003687D-5, 7.6484615D-5, 7.7033666D-5, 7.7674803D-5, 7.8447797D-5, 7.9426016D-5, 8.0771588D-5, 8.2988429D-5] * loge10
gp(1, 37:45)%dgp_d2     = [4.0635818D-2, 4.5423428D-2, 5.1163156D-2, 5.9022525D-2, 6.9958249D-2, 8.5784746D-2, 1.1269364D-1, 1.6431026D-1, 3.1533729D-1] * loge10

gp(1, 46:54)%gp         = [1.6191060D0, 1.6251334D0, 1.6319152D0, 1.6396776D0, 1.6487685D0, 1.6597654D0, 1.6737395D0, 1.6930618D0, 1.7251428D0]
gp(1, 46:54)%dgp_dprob  = [5.7092210D3, 6.3777055D3, 7.2314372D3, 8.3611515D3, 9.9290128D3, 1.2256681D4, 1.6087300D4, 2.3622468D4, 4.5675303D4]
gp(1, 46:54)%dgp_dlogE  = [8.3316766D-5, 8.3681234D-5, 8.4091282D-5, 8.4560304D-5, 8.5109371D-5, 8.5773444D-5, 8.6616385D-5, 8.7782783D-5, 8.9715588D-5] * loge10
gp(1, 46:54)%dgp_d2     = [3.4957427D-1, 4.0478160D-1, 4.1992958D-1, 4.8577941D-1, 6.0235126D-1, 7.5303943D-1, 1.0458175D0, 1.3405881D0, 2.6509833D0] * loge10

gp(1, 55:63)%gp         = [1.7299196D0, 1.7352295D0, 1.7412118D0, 1.7480686D0, 1.7561122D0, 1.7658609D0, 1.7782774D0, 1.7954984D0, 1.8242175D0]
gp(1, 55:63)%dgp_dprob  = [5.0266819D4, 5.6220367D4, 6.3831782D4, 7.3915207D4, 8.7928899D4, 1.0876706D5, 1.4312717D5, 2.1089407D5, 4.1000335D5]
gp(1, 55:63)%dgp_dlogE  = [9.0002195D-5, 9.0320911D-5, 9.0681626D-5, 9.1095876D-5, 9.1578628D-5, 9.2162177D-5, 9.2914257D-5, 9.3947535D-5, 9.5664389D-5] * loge10
gp(1, 55:63)%dgp_d2     = [3.1215901D0, 4.9569528D0, 3.9379037D0, 1.7469463D0, 5.4302755D0, 9.4968176D0, 4.7899651D0, 6.5800375D0, 2.5087735D1] * loge10

gp(2,  0: 9)%gp         = [0.0D0, 9.9147450D-2, 1.9396090D-1, 2.8284431D-1, 3.6677788D-1, 4.4806532D-1, 5.2978876D-1, 6.1627201D-1, 7.1526951D-1, 8.4790205D-1]
gp(2,  0: 9)%dgp_dprob  = [9.9999911D-1, 9.7493675D-1, 9.1824425D-1, 8.6074847D-1, 8.2146478D-1, 8.0923551D-1, 8.3239334D-1, 9.0994458D-1, 1.1007735D0, 1.6878403D0]
gp(2,  0: 9)%dgp_dlogE  = [-8.0680383D-28, 2.2889804D-5, 1.5259702D-4, 4.0562717D-4, 7.4762850D-4, 1.1480628D-3, 1.5954836D-3, 2.0984831D-3, 2.6944856D-3, 3.5074931D-3] * loge10
gp(2,  0: 9)%dgp_d2     = [2.4555676D-8, 6.6162046D-4, 1.9602469D-3, 3.0440748D-3, 3.7506750D-3, 4.2443047D-3, 4.7224239D-3, 5.4009993D-3, 6.7001305D-3, 1.0390214D-2] * loge10

gp(2, 10:18)%gp         = [8.6536581D-1, 8.8418885D-1, 9.0475327D-1, 9.2755841D-1, 9.5336209D-1, 9.8340311D-1, 1.0199370D0, 1.0678441D0, 1.1418942D0]
gp(2, 10:18)%dgp_dprob  = [1.8084405D0, 1.9634272D0, 2.1595918D0, 2.4163258D0, 2.7678524D0, 3.2810751D0, 4.1079027D0, 5.6885649D0, 1.0116548D1]
gp(2, 10:18)%dgp_dlogE  = [3.6148001D-3, 3.7307285D-3, 3.8574042D-3, 3.9978813D-3, 4.1567928D-3, 4.3417073D-3, 4.5663924D-3, 4.8605999D-3, 5.3142023D-3] * loge10
gp(2, 10:18)%dgp_d2     = [1.1136325D-2, 1.2093879D-2, 1.3303343D-2, 1.4883210D-2, 1.7042444D-2, 2.0189196D-2, 2.5249232D-2, 3.4901873D-2, 6.1863264D-2] * loge10

gp(2, 19:27)%gp         = [1.1524206D0, 1.1639863D0, 1.1768522D0, 1.1913945D0, 1.2081883D0, 1.2281797D0, 1.2531055D0, 1.2867608D0, 1.3407919D0]
gp(2, 19:27)%dgp_dprob  = [1.1014890D1, 1.2171057D1, 1.3636813D1, 1.5560545D1, 1.8205568D1, 2.2089882D1, 2.8398508D1, 4.0598423D1, 7.5406350D1]
gp(2, 19:27)%dgp_dlogE  = [5.3785567D-3, 5.4492281D-3, 5.5277960D-3, 5.6165415D-3, 5.7189457D-3, 5.8407346D-3, 5.9924086D-3, 6.1968915D-3, 6.5244281D-3] * loge10
gp(2, 19:27)%dgp_d2     = [6.7323984D-2, 7.4348571D-2, 8.3249385D-2, 9.4924100D-2, 1.1096522D-1, 1.3450299D-1, 1.7269313D-1, 2.4645032D-1, 4.5646500D-1] * loge10

gp(2, 28:36)%gp         = [1.3486581D0, 1.3573517D0, 1.3670837D0, 1.3781600D0, 1.3910492D0, 1.4065260D0, 1.4260192D0, 1.4526706D0, 1.4962055D0]
gp(2, 28:36)%dgp_dprob  = [8.2546366D1, 9.1764780D1, 1.0349341D2, 1.1894816D2, 1.4029407D2, 1.7180717D2, 2.2331569D2, 3.2374518D2, 6.1378173D2]
gp(2, 28:36)%dgp_dlogE  = [6.5720364D-3, 6.6246303D-3, 6.6834783D-3, 6.7504191D-3, 6.8282688D-3, 6.9216798D-3, 7.0392280D-3, 7.1997558D-3, 7.4615227D-3] * loge10
gp(2, 28:36)%dgp_d2     = [4.9949306D-1, 5.5502631D-1, 6.2564801D-1, 7.1866694D-1, 8.4708061D-1, 1.0365384D0, 1.3459654D0, 1.9486854D0, 3.6865706D0] * loge10

gp(2, 37:45)%gp         = [1.5026164D0, 1.5097224D0, 1.5177025D0, 1.5268169D0, 1.5374654D0, 1.5503102D0, 1.5665769D0, 1.5889710D0, 1.6259154D0]
gp(2, 37:45)%dgp_dprob  = [6.7369758D2, 7.5121087D2, 8.5005540D2, 9.8063265D2, 1.1615044D3, 1.4294241D3, 1.8691260D3, 2.7309625D3, 5.2395853D3]
gp(2, 37:45)%dgp_dlogE  = [7.5000233D-3, 7.5426842D-3, 7.5905759D-3, 7.6452527D-3, 7.7091020D-3, 7.7860781D-3, 7.8834945D-3, 8.0174869D-3, 8.2382410D-3] * loge10
gp(2, 37:45)%dgp_d2     = [4.0455571D0, 4.5091264D0, 5.1004195D0, 5.8815306D0, 6.9624877D0, 8.5637167D0, 1.1189761D1, 1.6332517D1, 3.1281497D1] * loge10

gp(2, 46:54)%gp         = [1.6313922D0, 1.6374734D0, 1.6443158D0, 1.6521474D0, 1.6613194D0, 1.6724144D0, 1.6865131D0, 1.7060075D0, 1.7383740D0]
gp(2, 46:54)%dgp_dprob  = [5.7602163D3, 6.4346562D3, 7.2959924D3, 8.4357659D3, 1.0017583D4, 1.2365964D4, 1.6230650D4, 2.3832794D4, 4.6081447D4]
gp(2, 46:54)%dgp_dlogE  = [8.2709356D-3, 8.3072291D-3, 8.3480536D-3, 8.3947661D-3, 8.4494522D-3, 8.5155754D-3, 8.5995549D-3, 8.7155922D-3, 8.9080436D-3] * loge10
gp(2, 46:54)%dgp_d2     = [3.4398418D1, 3.8388286D1, 4.3527217D1, 5.0296634D1, 5.9721144D1, 7.3705842D1, 9.6644534D1, 1.4180953D2, 2.7371082D2] * loge10

gp(2, 55:63)%gp         = [1.7431932D0, 1.7485504D0, 1.7545858D0, 1.7615036D0, 1.7696186D0, 1.7794538D0, 1.7919806D0, 1.8093544D0, 1.8383281D0]
gp(2, 55:63)%dgp_dprob  = [5.0713716D4, 5.6720085D4, 6.4398999D4, 7.4571939D4, 8.8709892D4, 1.0973270D5, 1.4439753D5, 2.1276473D5, 4.1363718D5]
gp(2, 55:63)%dgp_dlogE  = [8.9366753D-3, 8.9685012D-3, 9.0043436D-3, 9.0454219D-3, 9.0935942D-3, 9.1519542D-3, 9.2262593D-3, 9.3292477D-3, 9.5008889D-3] * loge10
gp(2, 55:63)%dgp_d2     = [3.0283275D2, 3.3614161D2, 3.8291963D2, 4.4276651D2, 5.2706418D2, 6.5031838D2, 8.5645492D2, 1.2606586D3, 2.4482288D3] * loge10

gp(3,  0: 9)%gp         = [0.0D0, 9.9262787D-2, 1.9473650D-1, 2.8492560D-1, 3.7064553D-1, 4.5404365D-1, 5.3814050D-1, 6.2730371D-1, 7.2948636D-1, 8.6647312D-1]
gp(3,  0: 9)%dgp_dprob  = [9.9999923D-1, 9.7827882D-1, 9.2828321D-1, 8.7655599D-1, 8.4115349D-1, 8.3168666D-1, 8.5750542D-1, 9.3876918D-1, 1.1366216D0, 1.7435311D0]
gp(3,  0: 9)%dgp_dlogE  = [4.0963128D-25, 9.2133096D-5, 6.2499896D-4, 1.6932994D-3, 3.1729432D-3, 4.9374202D-3, 6.9344424D-3, 9.1995407D-3, 1.1900030D-2, 1.5599885D-2] * loge10
gp(3,  0: 9)%dgp_d2     = [9.8100235D-8, 2.6763520D-3, 8.1511523D-3, 1.3016103D-2, 1.6390697D-2, 1.8836500D-2, 2.1182471D-2, 2.4404295D-2, 3.0429354D-2, 4.7360165D-2] * loge10

gp(3, 10:18)%gp         = [8.8463536D-1, 9.0408397D-1, 9.2533207D-1, 9.4889548D-1, 9.7555705D-1, 1.0065965D0, 1.0443437D0, 1.0938397D0, 1.1703399D0]
gp(3, 10:18)%dgp_dprob  = [1.8685307D0, 2.0286939D0, 2.2313958D0, 2.4966682D0, 2.8598633D0, 3.3900897D0, 4.2442576D0, 5.8770720D0, 1.0450718D1]
gp(3, 10:18)%dgp_dlogE  = [1.6089075D-2, 1.6617728D-2, 1.7195559D-2, 1.7836529D-2, 1.8561814D-2, 1.9406011D-2, 2.0432051D-2, 2.1775929D-2, 2.3848448D-2] * loge10
gp(3, 10:18)%dgp_d2     = [5.0775857D-2, 5.5158194D-2, 6.0691918D-2, 6.7918689D-2, 7.7793707D-2, 9.2182660D-2, 1.1531727D-1, 1.5944463D-1, 2.8268866D-1] * loge10

gp(3, 19:27)%gp         = [1.1812139D0, 1.1931614D0, 1.2064517D0, 1.2214734D0, 1.2388204D0, 1.2594698D0, 1.2852149D0, 1.3199748D0, 1.3757754D0]
gp(3, 19:27)%dgp_dprob  = [1.1378561D1, 1.2572679D1, 1.4086521D1, 1.6073326D1, 1.8805009D1, 2.2816484D1, 2.9331427D1, 4.1929799D1, 7.7872260D1]
gp(3, 19:27)%dgp_dlogE  = [2.4142524D-2, 2.4465474D-2, 2.4824519D-2, 2.5230084D-2, 2.5698080D-2, 2.6254680D-2, 2.6947875D-2, 2.7882439D-2, 2.9379416D-2] * loge10
gp(3, 19:27)%dgp_d2     = [3.0764942D-1, 3.3975836D-1, 3.8044289D-1, 4.3380655D-1, 5.0712700D-1, 6.1471189D-1, 7.8926660D-1, 1.1263781D0, 2.0862248D0] * loge10

gp(3, 28:36)%gp         = [1.3838987D0, 1.3928765D0, 1.4029264D0, 1.4143643D0, 1.4276741D0, 1.4436555D0, 1.4637837D0, 1.4913023D0, 1.5362512D0]
gp(3, 28:36)%dgp_dprob  = [8.5244713D1, 9.4763116D1, 1.0687327D2, 1.2283051D2, 1.4487010D2, 1.7740664D2, 2.3058668D2, 3.3427191D2, 6.3369600D2]
gp(3, 28:36)%dgp_dlogE  = [2.9597005D-2, 2.9837378D-2, 3.0106335D-2, 3.0412277D-2, 3.0768073D-2, 3.1194986D-2, 3.1732205D-2, 3.2465834D-2, 3.3662093D-2] * loge10
gp(3, 28:36)%dgp_d2     = [2.2828705D0, 2.5366655D0, 2.8594366D0, 3.2845513D0, 3.8714004D0, 4.7372257D0, 6.1513203D0, 8.9055964D0, 1.6847035D1] * loge10

gp(3, 37:45)%gp         = [1.5428701D0, 1.5502065D0, 1.5584453D0, 1.5678551D0, 1.5788484D0, 1.5921090D0, 1.6089019D0, 1.6320199D0, 1.6701566D0]
gp(3, 37:45)%dgp_dprob  = [6.9554926D2, 7.7556828D2, 8.7760709D2, 1.0124026D3, 1.1991148D3, 1.4756816D3, 1.9295661D3, 2.8191786D3, 5.4085484D3]
gp(3, 37:45)%dgp_dlogE  = [3.3838034D-2, 3.4032985D-2, 3.4251837D-2, 3.4501693D-2, 3.4793460D-2, 3.5145206D-2, 3.5590347D-2, 3.6202604D-2, 3.7211259D-2] * loge10
gp(3, 37:45)%dgp_d2     = [1.8485962D1, 2.0605820D1, 2.3307715D1, 2.6876126D1, 3.1816947D1, 3.9131781D1, 5.1128872D1, 7.4624582D1, 1.4292637D2] * loge10

gp(3, 46:54)%gp         = [1.6758100D0, 1.6820872D0, 1.6891501D0, 1.6972341D0, 1.7067014D0, 1.7181536D0, 1.7327058D0, 1.7528269D0, 1.7862328D0]
gp(3, 46:54)%dgp_dprob  = [5.9459244D3, 6.6420519D3, 7.5310774D3, 8.7074809D3, 1.0340117D4, 1.2763914D4, 1.6752648D4, 2.4598640D4, 4.7560198D4]
gp(3, 46:54)%dgp_dlogE  = [3.7360640D-2, 3.7526463D-2, 3.7712986D-2, 3.7926407D-2, 3.8176257D-2, 3.8478354D-2, 3.8862020D-2, 3.9392128D-2, 4.0271280D-2] * loge10
gp(3, 46:54)%dgp_d2     = [1.5710284D2, 1.7545320D2, 1.9887836D2, 2.2982888D2, 2.7283147D2, 3.3661673D2, 4.4155293D2, 6.4782225D2, 1.2506197D3] * loge10

gp(3, 55:63)%gp         = [1.7912066D0, 1.7967356D0, 1.8029646D0, 1.8101042D0, 1.8184794D0, 1.8286297D0, 1.8415577D0, 1.8594875D0, 1.8893877D0]
gp(3, 55:63)%dgp_dprob  = [5.2340786D4, 5.8539439D4, 6.6464178D4, 7.6962669D4, 9.1552893D4, 1.1324808D5, 1.4902101D5, 2.1957209D5, 4.2685519D5]
gp(3, 55:63)%dgp_dlogE  = [4.0402076D-2, 4.0547443D-2, 4.0711173D-2, 4.0898804D-2, 4.1118826D-2, 4.1385392D-2, 4.1724755D-2, 4.2195092D-2, 4.2978781D-2] * loge10
gp(3, 55:63)%dgp_d2     = [1.3758086D3, 1.5389513D3, 1.7478905D3, 2.0220601D3, 2.4041342D3, 2.9717225D3, 3.9097451D3, 5.7576231D3, 1.1182549D4] * loge10

gp(4,  0: 9)%gp         = [0.0d0, 9.9633250D-2, 1.9732491D-1, 2.9218521D-1, 3.8469113D-1, 4.7649976D-1, 5.7039384D-1, 6.7090439D-1, 7.8682229D-1, 9.4283958D-1]
gp(4,  0: 9)%dgp_dprob  = [9.9999962D-1, 9.8912804D-1, 9.6293031D-1, 9.3484921D-1, 9.1788916D-1, 9.2289610D-1, 9.6262321D-1, 1.0620370D0, 1.2923099D0, 1.9881753D0]
gp(4,  0: 9)%dgp_dlogE  = [1.3278291D-22, 2.3226179D-4, 1.6727453D-3, 4.8645037D-3, 9.7377811D-3, 1.6029279D-2, 2.3588264D-2, 3.2544803D-2, 4.3569085D-2, 5.9037015D-2] * loge10
gp(4,  0: 9)%dgp_d2     = [2.4133968D-7, 6.8578237D-3, 2.2986385D-2, 4.0824471D-2, 5.6307503D-2, 6.9380907D-2, 8.2113414D-2, 9.8178061D-2, 1.2578515D-1, 1.9981487D-1] * loge10

gp(4, 10:18)%gp         = [9.6397419D-1, 9.8617378D-1, 1.0104313D0, 1.0373356D0, 1.0677805D0, 1.1032271D0, 1.1463352D0, 1.2028584D0, 1.2902054D0]
gp(4, 10:18)%dgp_dprob  = [2.1326462D0, 2.3158418D0, 2.5476087D0, 2.8508176D0, 3.2658272D0, 3.8715246D0, 4.8469999D0, 6.7111504D0, 1.1930820D1]
gp(4, 10:18)%dgp_dlogE  = [6.1102625D-2, 6.3338847D-2, 6.5787389D-2, 6.8508189D-2, 7.1592199D-2, 7.5188023D-2, 7.9566026D-2, 8.5310534D-2, 9.4186872D-2] * loge10
gp(4, 10:18)%dgp_d2     = [2.1459325D-1, 2.3352604D-1, 2.5740203D-1, 2.8854846D-1, 3.3106941D-1, 3.9298127D-1, 4.9246770D-1, 6.8215810D-1, 1.2118528D0] * loge10

gp(4, 19:27)%gp         = [1.3026192D0, 1.3162579D0, 1.3314284D0, 1.3485744D0, 1.3683728D0, 1.3919380D0, 1.4213150D0, 1.4609718D0, 1.5246172D0]
gp(4, 19:27)%dgp_dprob  = [1.2989500D1, 1.4351922D1, 1.6079019D1, 1.8345537D1, 2.1461539D1, 2.6036932D1, 3.3466818D1, 4.7832172D1, 8.8805569D1]
gp(4, 19:27)%dgp_dlogE  = [9.5447665D-2, 9.6832557D-2, 9.8372579D-2, 1.0011255D-1, 1.0212084D-1, 1.0450996D-1, 1.0748619D-1, 1.1149987D-1, 1.1793092D-1] * loge10
gp(4, 19:27)%dgp_d2     = [1.3191312D0, 1.4571323D0, 1.6319918D0, 1.8613480D0, 2.1764822D0, 2.6388935D0, 3.3891591D0, 4.8381393D0, 8.9637506D0] * loge10

gp(4, 28:36)%gp         = [1.5338809D0, 1.5441184D0, 1.5555779D0, 1.5686192D0, 1.5837937D0, 1.6020125D0, 1.6249562D0, 1.6563196D0, 1.7075378D0]
gp(4, 28:36)%dgp_dprob  = [9.7208696D1, 1.0805730D2, 1.2185915D2, 1.4004449D2, 1.6515980D2, 2.0223407D2, 2.6282508D2, 3.8094444D2, 7.2198396D2]
gp(4, 28:36)%dgp_dlogE  = [1.1886584D-1, 1.1989868D-1, 1.2105438D-1, 1.2236904D-1, 1.2389797D-1, 1.2573255D-1, 1.2804120D-1, 1.3119393D-1, 1.3633468D-1] * loge10
gp(4, 28:36)%dgp_d2     = [9.8089592D0, 1.0899794D1, 1.2287077D1, 1.4114231D1, 1.6636462D1, 2.0357598D1, 2.6434863D1, 3.8271182D1, 7.2395629D1] * loge10

gp(4, 37:45)%gp         = [1.7150787D0, 1.7234367D0, 1.7328224D0, 1.7435416D0, 1.7560639D0, 1.7711677D0, 1.7902931D0, 1.8166191D0, 1.8600407D0]
gp(4, 37:45)%dgp_dprob  = [7.9242542D2, 8.8355168D2, 9.9974954D2, 1.1532424D3, 1.3658422D3, 1.6807351D3, 2.1974778D3, 3.2101875D3, 6.1573770D3]
gp(4, 37:45)%dgp_dlogE  = [1.3709073D-1, 1.3792846D-1, 1.3886889D-1, 1.3994252D-1, 1.4119622D-1, 1.4270760D-1, 1.4462018D-1, 1.4725064D-1, 1.5158367D-1] * loge10
gp(4, 37:45)%dgp_d2     = [7.9437736D1, 8.8545261D1, 1.0015495D2, 1.1548566D2, 1.3671196D2, 1.6813644D2, 2.1967444D2, 3.2059861D2, 6.1395018D2] * loge10

gp(4, 46:54)%gp         = [1.8664767D0, 1.8736227D0, 1.8816628D0, 1.8908649D0, 1.9016411D0, 1.9146759D0, 1.9312379D0, 1.9541359D0, 1.9921466D0]
gp(4, 46:54)%dgp_dprob  = [6.7689521D3, 7.5611752D3, 8.5728907D3, 9.9115945D3, 1.1769394D4, 1.4527324D4, 1.9065636D4, 2.7991907D4, 5.4111428D4]
gp(4, 46:54)%dgp_dlogE  = [1.5222534D-1, 1.5293761D-1, 1.5373878D-1, 1.5465545D-1, 1.5572854D-1, 1.5702598D-1, 1.5867364D-1, 1.6095001D-1, 1.6472474D-1] * loge10
gp(4, 46:54)%dgp_d2     = [6.7480511D2, 7.5357957D2, 8.5414512D2, 9.8717992D2, 1.1717486D3, 1.4456760D3, 1.8961779D3, 2.7816531D3, 5.3700030D3] * loge10

gp(4, 55:63)%gp         = [1.9978055D0, 2.0040959D0, 2.0111823D0, 2.0193047D0, 2.0288321D0, 2.0403785D0, 2.0550837D0, 2.0754764D0, 2.1094800D0]
gp(4, 55:63)%dgp_dprob  = [5.9549048D4, 6.6599513D4, 7.5612720D4, 8.7553530D4, 1.0414699D5, 1.2882006D5, 1.6950095D5, 2.4972440D5, 4.8539673D5]
gp(4, 55:63)%dgp_dlogE  = [1.6528628D-1, 1.6591035D-1, 1.6661326D-1, 1.6741869D-1, 1.6836317D-1, 1.6950738D-1, 1.7096398D-1, 1.7298261D-1, 1.7634572D-1] * loge10
gp(4, 55:63)%dgp_d2     = [5.9085142D3, 6.6057974D3, 7.4989769D3, 8.6837549D3, 1.0319646D4, 1.2762311D4, 1.6783997D4, 2.4715000D4, 4.7963032D4] * loge10

gp(5,  0: 9)%gp         = [0.0D0, 1.0012393D-1, 2.0102183D-1, 3.0362575D-1, 4.0921862D-1, 5.1972735D-1, 6.3827719D-1, 7.7048035D-1, 9.2824622D-1, 1.1467214D0]
gp(5,  0: 9)%dgp_dprob  = [1.0000001D0, 1.0037812D0, 1.0158217D0, 1.0385804D0, 1.0767936D0, 1.1390392D0, 1.2418638D0, 1.4221098D0, 1.7845765D0, 2.8166667D0]
gp(5,  0: 9)%dgp_dlogE  = [-1.7982465D-15, 1.3765591D-4, 1.0882687D-3, 3.6078644D-3, 8.3732115D-3, 1.6021356D-2, 2.7280858D-2, 4.3279444D-2, 6.6389090D-2, 1.0384043D-1] * loge10
gp(5,  0: 9)%dgp_d2     = [1.3822934D-7, 4.1594387D-3, 1.6237174D-2, 3.5468186D-2, 6.1136643D-2, 9.3318250D-2, 1.3405032D-1, 1.9006085D-1, 2.8281439D-1, 5.1500096D-1] * loge10

gp(5, 10:18)%gp         = [1.1769144D0, 1.2084988D0, 1.2430918D0, 1.2815487D0, 1.3251676D0, 1.3760708D0, 1.4381223D0, 1.5196819D0, 1.6460439D0]
gp(5, 10:18)%dgp_dprob  = [3.0306347D0, 3.2986983D0, 3.6372839D0, 4.0796210D0, 4.6843482D0, 5.5660852D0, 6.9850476D0, 9.6952052D0, 1.7280983D1]
gp(5, 10:18)%dgp_dlogE  = [1.0919745D-1, 1.1507682D-1, 1.2160431D-1, 1.2896097D-1, 1.3742205D-1, 1.4743832D-1, 1.5983209D-1, 1.7638659D-1, 2.0251464D-1] * loge10
gp(5, 10:18)%dgp_d2     = [5.6032240D-1, 6.1819805D-1, 6.9099075D-1, 7.8576333D-1, 9.1499470D-1, 1.1031123D0, 1.4056682D0, 1.9839704D0, 3.6072320D0] * loge10

gp(5, 19:27)%gp         = [1.6640266D0, 1.6837894D0, 1.7057781D0, 1.7306372D0, 1.7593505D0, 1.7935370D0, 1.8361676D0, 1.8937334D0, 1.9861469D0]
gp(5, 19:27)%dgp_dprob  = [1.8819410D1, 2.0799182D1, 2.3308819D1, 2.6602214D1, 3.1129844D1, 3.7777788D1, 4.8572743D1, 6.9442718D1, 1.2895948D2]
gp(5, 19:27)%dgp_dlogE  = [2.0627234D-1, 2.1041184D-1, 2.1502913D-1, 2.2026297D-1, 2.2632540D-1, 2.3356560D-1, 2.4262499D-1, 2.5490640D-1, 2.7471951D-1] * loge10
gp(5, 19:27)%dgp_d2     = [3.9371159D0, 4.3619024D0, 4.9007593D0, 5.6084599D0, 6.5822509D0, 8.0135485D0, 1.0340529D1, 1.4846092D1, 2.7721594D1] * loge10

gp(5, 28:36)%gp         = [1.9995993D0, 2.0144661D0, 2.0311076D0, 2.0500462D0, 2.0720825D0, 2.0985394D0, 2.1318563D0, 2.1773962D0, 2.2517530D0]
gp(5, 28:36)%dgp_dprob  = [1.4116420D2, 1.5694264D2, 1.7696428D2, 2.0337304D2, 2.3984308D2, 2.9367431D2, 3.8164182D2, 5.5310336D2, 1.0480207D3]
gp(5, 28:36)%dgp_dlogE  = [2.7761218D-1, 2.8081120D-1, 2.8439477D-1, 2.8847625D-1, 2.9322939D-1, 2.9894144D-1, 3.0614224D-1, 3.1599705D-1, 3.3211315D-1] * loge10
gp(5, 28:36)%dgp_d2     = [3.0364840D1, 3.4914961D1, 3.8121986D1, 4.3847051D1, 5.1756335D1, 6.3435713D1, 8.2530630D1, 1.1977048D2, 2.2733955D2] * loge10

gp(5, 37:45)%gp         = [2.2626990D0, 2.2748307D0, 2.2884533D0, 2.3040103D0, 2.3221830D0, 2.3441002D0, 2.3718498D0, 2.4100408D0, 2.4730152D0]
gp(5, 37:45)%dgp_dprob  = [1.1502275D3, 1.2824396D3, 1.4510166D3, 1.6736841D3, 1.9820698D3, 2.4387886D3, 3.1881697D3, 4.6565335D3, 8.9284860D3]
gp(5, 37:45)%dgp_dlogE  = [3.3448783D-1, 3.3712030D-1, 3.4007699D-1, 3.4345438D-1, 3.4740067D-1, 3.5216145D-1, 3.5819109D-1, 3.6649244D-1, 3.8018646D-1] * loge10
gp(5, 37:45)%dgp_d2     = [2.4956220D2, 2.7831093D2, 3.1497102D2, 3.6339887D2, 4.3047459D2, 5.2982689D2, 6.9286051D2, 1.0123467D3, 1.9419101D3] * loge10

gp(5, 46:54)%gp         = [2.4823476D0, 2.4927088D0, 2.5043656D0, 2.5177061D0, 2.5333275D0, 2.5522207D0, 2.5762235D0, 2.6094028D0, 2.6644650D0]
gp(5, 46:54)%dgp_dprob  = [9.8147976D3, 1.0962842D4, 1.2428869D4, 1.4368570D4, 1.7060173D4, 2.1055474D4, 2.7629065D4, 4.0556118D4, 7.8371409D4]
gp(5, 46:54)%dgp_dlogE  = [3.8221627D-1, 3.8446993D-1, 3.8700552D-1, 3.8990748D-1, 3.9330572D-1, 3.9741587D-1, 4.0263769D-1, 4.0985586D-1, 4.2183387D-1] * loge10
gp(5, 46:54)%dgp_d2     = [2.1348142D3, 2.3845768D3, 2.7036308D3, 3.1256589D3, 3.7113103D3, 4.5805743D3, 6.0107958D3, 8.8229690D3, 1.7047401D4] * loge10

gp(5, 55:63)%gp         = [2.6726608D0, 2.6817706D0, 2.6920328D0, 2.7037941D0, 2.7175889D0, 2.7343052D0, 2.7555920D0, 2.7851064D0, 2.8343073D0]
gp(5, 55:63)%dgp_dprob  = [8.6242483D4, 9.6447634D4, 1.0949310D5, 1.2677428D5, 1.5078747D5, 1.8648914D5, 2.4534684D5, 3.6139369D5, 7.0221984D5]
gp(5, 55:63)%dgp_dlogE  = [4.2361659D-1, 4.2559808D-1, 4.2783012D-1, 4.3038813D-1, 4.3338820D-1, 4.3702335D-1, 4.4165192D-1, 4.4806828D-1, 4.5876140D-1] * loge10
gp(5, 55:63)%dgp_d2     = [1.8759846D4, 2.0984283D4, 2.3814991D4, 2.7566799D4, 3.2791466D4, 4.0553059D4, 5.3346565D4, 7.8569428D4, 1.5258899D5] * loge10

gp(6,  0: 9)%gp         = [0.0D0, 1.0025385D-1, 2.0206489D-1, 3.0717165D-1, 4.1773440D-1, 5.3673379D-1, 6.6875939D-1, 8.2186579D-1, 1.0130813D0, 1.2930532D0]
gp(6,  0: 9)%dgp_dprob  = [1.0000003D0, 1.0077208D0, 1.0316003D0, 1.0745069D0, 1.1422832D0, 1.2461426D0, 1.4087318D0, 1.6819846D0, 2.2162445D0, 3.7204800D0]
gp(6,  0: 9)%dgp_dlogE  = [-1.3922505D-15, 1.1625058D-5, 9.4951144D-5, 3.3214559D-4, 8.3012960D-4, 1.7453305D-3, 3.3350156D-3, 6.0844109D-3, 1.1106553D-2, 2.2073570D-2] * loge10
gp(6,  0: 9)%dgp_d2     = [1.1569362D-8, 3.5389579D-4, 1.4584941D-3, 3.4802666D-3, 6.7658443D-3, 1.2002371D-2, 2.0637323D-2, 3.6169849D-2, 6.9429376D-2, 1.7651448D-1] * loge10

gp(6, 10:18)%gp         = [1.3318805D0, 1.3740111D0, 1.4205356D0, 1.4727144D0, 1.5324661D0, 1.6029386D0, 1.6898871D0, 1.8058275D0, 1.9889084D0]
gp(6, 10:18)%dgp_dprob  = [4.0266212D0, 4.4177768D0, 4.9126991D0, 5.5607160D0, 6.4491255D0, 7.7491358D0, 9.8509771D0, 1.3891022D1, 2.5311273D1]
gp(6, 10:18)%dgp_dlogE  = [2.3944894D-2, 2.6089056D-2, 2.8581614D-2, 3.1533725D-2, 3.5118503D-2, 3.9627296D-2, 4.5608816D-2, 5.4300506D-2, 6.9666526D-2] * loge10
gp(6, 10:18)%dgp_d2     = [1.9986587D-1, 2.3051285D-1, 2.7021136D-1, 3.2355176D-1, 3.9882672D-1, 5.1268395D-1, 7.0412160D-1, 1.0906384D0, 2.2629433D0] * loge10

gp(6, 19:27)%gp         = [2.0152821D0, 2.0443536D0, 2.0768051D0, 2.1136236D0, 2.1563192D0, 2.2073826D0, 2.2713961D0, 2.3584041D0, 2.4993604D0]
gp(6, 19:27)%dgp_dprob  = [2.7641011D1, 3.0644165D1, 3.4458302D1, 3.9474189D1, 4.6386413D1, 5.6564060D1, 7.3146360D1, 1.0534422D2, 1.9775510D2]
gp(6, 19:27)%dgp_dlogE  = [7.2042751D-2, 7.4708757D-2, 7.7742190D-2, 8.1256690D-2, 8.5428262D-2, 9.0551348D-2, 9.7176787D-2, 1.0653727D-1, 1.2254346D-1] * loge10
gp(6, 19:27)%dgp_d2     = [2.5117983D0, 2.8362039D0, 3.2534037D0, 3.8097184D0, 4.5884441D0, 5.7559555D0, 7.6997986D0, 1.1579911D1, 2.3177480D1] * loge10

gp(6, 28:36)%gp         = [2.5200030D0, 2.5428509D0, 2.5684694D0, 2.5976784D0, 2.6317361D0, 2.6727245D0, 2.7244896D0, 2.7955030D0, 2.9120560D0]
gp(6, 28:36)%dgp_dprob  = [2.1677584D2, 2.4135714D2, 2.7266591D2, 3.1397097D2, 3.7109869D2, 4.5556972D2, 5.9390022D2, 8.6426631D2, 1.6478591D3]
gp(6, 28:36)%dgp_dlogE  = [1.2497248D-1, 1.2768577D-1, 1.3075875D-1, 1.3430161D-1, 1.3848474D-1, 1.4359259D-1, 1.5015617D-1, 1.5936073D-1, 1.7495419D-1] * loge10
gp(6, 28:36)%dgp_d2     = [2.5621083D1, 2.8800314D1, 3.2880172D1, 3.8307949D1, 4.5886811D1, 5.7218168D1, 7.6025355D1, 1.1342765D2, 2.2468156D2] * loge10

gp(6, 37:45)%gp         = [2.9292739D0, 2.9483739D0, 2.9698428D0, 2.9943876D0, 3.0230953D0, 3.0577688D0, 3.1017470D0, 3.1624093D0, 3.2627653D0]
gp(6, 37:45)%dgp_dprob  = [1.8100674D3, 2.0200399D3, 2.2879721D3, 2.6421801D3, 3.1332266D3, 3.8613024D3, 5.0575934D3, 7.4059043D3, 1.4256584D4]
gp(6, 37:45)%dgp_dlogE  = [1.7730749D-1, 1.7993272D-1, 1.8290176D-1, 1.8631963D-1, 1.9034855D-1, 1.9525922D-1, 2.0155655D-1, 2.1036655D-1, 2.2524567D-1] * loge10
gp(6, 37:45)%dgp_d2     = [2.4806410D2, 2.7846596D2, 3.1744890D2, 3.6927572D2, 4.4157485D2, 5.4957608D2, 7.2863143D2, 1.0842710D3, 2.1403083D3] * loge10

gp(6, 46:54)%gp         = [3.2776705D0, 3.2942287D0, 3.3128694D0, 3.3342180D0, 3.3592372D0, 3.3895257D0, 3.4280504D0, 3.4813834D0, 3.5700872D0]
gp(6, 46:54)%dgp_dprob  = [1.5680191D4, 1.7525058D4, 1.9882150D4, 2.3002655D4, 2.7335699D4, 3.3772560D4, 4.4373541D4, 6.5246796D4, 1.2642431D5]
gp(6, 46:54)%dgp_dlogE  = [2.2748697D-1, 2.2998611D-1, 2.3281116D-1, 2.3606154D-1, 2.3989078D-1, 2.4455493D-1, 2.5053138D-1, 2.5888392D-1, 2.7296712D-1] * loge10
gp(6, 46:54)%dgp_d2     = [2.3621284D3, 2.6501803D3, 3.0198042D3, 3.5108395D3, 4.1956921D3, 5.2179844D3, 6.9121717D3, 1.0273544D4, 2.0236332D4] * loge10

gp(6, 55:63)%gp         = [3.5833105D0, 3.5980144D0, 3.6145857D0, 3.6335874D0, 3.6558872D0, 3.6829278D0, 3.7173902D0, 3.7652242D0, 3.8450944D0]
gp(6, 55:63)%dgp_dprob  = [1.3917275D5, 1.5570667D5, 1.7685091D5, 2.0487219D5, 2.4382853D5, 3.0178063D5, 3.9738915D5, 5.8607760D5, 1.1411139D6]
gp(6, 55:63)%dgp_dlogE  = [2.7508571D-1, 2.7744699D-1, 2.8011489D-1, 2.8318240D-1, 2.8679317D-1, 2.9118608D-1, 2.9680486D-1, 3.0463214D-1, 3.1772353D-1] * loge10
gp(6, 55:63)%dgp_d2     = [2.2339388D4, 2.5036071D4, 2.8512800D4, 3.3124964D4, 3.9542163D4, 4.9116286D4, 6.4874905D4, 9.6041328D4, 1.8671917D5] * loge10

gp(7,  0: 9)%gp         = [0.0D0, 1.0026255D-1, 2.0213604D-1, 3.0742121D-1, 4.1836056D-1, 5.3805752D-1, 6.7130798D-1, 8.2656494D-1, 1.0217952D0, 1.3108530D0]
gp(7,  0: 9)%dgp_dprob  = [1.0000003D0, 1.0079857D0, 1.0326947D0, 1.0771300D0, 1.1474159D0, 1.2553312D0, 1.4247300D0, 1.7105183D0, 2.2724506D0, 3.8699053D0]
gp(7,  0: 9)%dgp_dlogE  = [-2.0182050D-13, 7.3746879D-7, 6.0421140D-6, 2.1248958D-5, 5.3531101D-5, 1.1381021D-4, 2.2084955D-4, 4.1175366D-4, 7.7640081D-4, 1.6335956D-3] * loge10
gp(7,  0: 9)%dgp_d2     = [9.3338668D-10, 2.2465967D-5, 9.3066137D-5, 2.2406304D-4, 4.4138977D-4, 7.9763817D-4, 1.4069138D-3, 2.5554961D-3, 5.1728179D-3, 1.4424635D-2] * loge10

gp(7, 10:18)%gp         = [1.3510522D0, 1.3949970D0, 1.4436473D0, 1.4983667D0, 1.5612358D0, 1.6356780D0, 1.7279738D0, 1.8518423D0, 2.0493729D0]
gp(7, 10:18)%dgp_dprob  = [4.1951001D0, 4.6135517D0, 5.1441312D0, 5.8405117D0, 6.7979084D0, 8.2035806D0, 1.0485810D1, 1.4897384D1, 2.7479174D1]
gp(7, 10:18)%dgp_dlogE  = [1.7874870D-3, 1.9664171D-3, 2.1778588D-3, 2.4329901D-3, 2.7495491D-3, 3.1580446D-3, 3.7174465D-3, 4.5652036D-3, 6.1632050D-3] * loge10
gp(7, 10:18)%dgp_d2     = [1.6551148D-2, 1.9385701D-2, 2.3122290D-2, 2.8242237D-2, 3.5631125D-2, 4.7104499D-2, 6.7024784D-2, 1.0897483D-1, 2.4470995D-1] * loge10

gp(7, 19:27)%gp         = [2.0780287D0, 2.1096753D0, 2.1450755D0, 2.1853352D0, 2.2321493D0, 2.2883206D0, 2.3590200D0, 2.4556253D0, 2.6133921D0]
gp(7, 19:27)%dgp_dprob  = [3.0059911D1, 3.3391991D1, 3.7631670D1, 4.3218821D1, 5.0936896D1, 6.2333802D1, 8.0969072D1, 1.1732637D2, 2.2246048D2]
gp(7, 19:27)%dgp_dlogE  = [6.4215147D-3, 6.7148728D-3, 7.0532165D-3, 7.4512713D-3, 7.9321693D-3, 8.5352512D-3, 9.3355945D-3, 1.0505552D-2, 1.2612747D-2] * loge10
gp(7, 19:27)%dgp_d2     = [2.7464046D-1, 3.1410115D-1, 3.6550571D-1, 4.3504961D-1, 5.3403503D-1, 6.8540157D-1, 9.4365290D-1, 1.4761470D0, 3.1503199D0] * loge10

gp(7, 28:36)%gp         = [2.6366289D0, 2.6623879D0, 2.6913201D0, 2.7243716D0, 2.7629962D0, 2.8096049D0, 2.8686621D0, 2.9500330D0, 3.0844764D0]
gp(7, 28:36)%dgp_dprob  = [2.4419929D2, 2.7233140D2, 3.0821821D2, 3.5564606D2, 4.2137524D2, 5.1879984D2, 6.7882333D2, 9.9284837D2, 1.9087529D3]
gp(7, 28:36)%dgp_dlogE  = [1.2944311D-2, 1.3318358D-2, 1.3746689D-2, 1.4246716D-2, 1.4845670D-2, 1.5589632D-2, 1.6566020D-2, 1.7974080D-2, 2.0463525D-2] * loge10
gp(7, 28:36)%dgp_d2     = [3.5139131D0, 3.9912722D0, 4.6101494D0, 5.4430528D0, 6.6218448D0, 8.4124143D0, 1.1444541D1, 1.7634528D1, 3.6831199D1] * loge10

gp(7, 37:45)%gp         = [3.1044311D0, 3.1265958D0, 3.1515449D0, 3.1801148D0, 3.2135930D0, 3.2541185D0, 3.3056608D0, 3.3770174D0, 3.4957282D0]
gp(7, 37:45)%dgp_dprob  = [2.0990837D3, 2.3457434D3, 2.6609017D3, 3.0781622D3, 3.6576156D3, 4.5185367D3, 5.9367232D3, 8.7301630D3, 1.6923554D4]
gp(7, 37:45)%dgp_dlogE  = [2.0850656D-2, 2.1286069D-2, 2.1783036D-2, 2.2361089D-2, 2.3050699D-2, 2.3903287D-2, 2.5016081D-2, 2.6609763D-2, 2.9399804D-2] * loge10
gp(7, 37:45)%dgp_d2     = [4.0950218D1, 4.6393472D1, 5.3399508D1, 6.2805572D1, 7.6068555D1, 9.6132572D1, 1.3001688D2, 1.9892242D2, 4.1044836D2] * loge10

gp(7, 46:54)%gp         = [3.5134300D0, 3.5331161D0, 3.5553051D0, 3.5807524D0, 3.6106225D0, 3.6468520D0, 3.6930414D0, 3.7571851D0, 3.8643813D0]
gp(7, 46:54)%dgp_dprob  = [1.8631809D4, 2.0847736D4, 2.3682106D4, 2.7439209D4, 3.2663964D4, 4.0439172D4, 5.3272447D4, 7.8615735D4, 1.5324133D5]
gp(7, 46:54)%dgp_dlogE  = [2.9830895D-2, 3.0314922D-2, 3.0866327D-2, 3.1506410D-2, 3.2268225D-2, 3.3207566D-2, 3.4429541D-2, 3.6172316D-2, 3.9204794D-2] * loge10
gp(7, 46:54)%dgp_d2     = [4.5499725D2, 5.1494300D2, 5.9304845D2, 6.9536844D2, 8.3973247D2, 1.0544056D3, 1.4222571D3, 2.1671186D3, 4.4432781D3] * loge10

gp(7, 55:63)%gp         = [3.8804158D0, 3.8982625D0, 3.9183965D0, 3.9415107D0, 3.9686737D0, 4.0016651D0, 4.0437962D0, 4.1024309D0, 4.2007361D0]
gp(7, 55:63)%dgp_dprob  = [1.6883662D5, 1.8907822D5, 2.1499248D5, 2.4937023D5, 2.9722818D5, 3.6852126D5, 4.8637160D5, 7.1952501D5, 1.4080893D6]
gp(7, 55:63)%dgp_dlogE  = [3.9671360D-2, 4.0194842D-2, 4.0790174D-2, 4.1480649D-2, 4.2300878D-2, 4.3310518D-2, 4.4620932D-2, 4.6483830D-2, 4.9710565D-2] * loge10
gp(7, 55:63)%dgp_d2     = [5.0182470D3, 5.4072400D3, 6.5108714D3, 7.5116126D3, 9.0615626D3, 1.1139837D4, 1.5435667D4, 2.3362818D4, 4.7381880D4] * loge10

end subroutine init_this

end function bend_photon_vert_angle_init

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!+
! Function bend_photon_e_rel_init (r_in) result (E_rel)
!
! Routine to convert a random number in the interval [0,1] to a photon energy.
! The photon probability spectrum is:
!   P(E_rel) = (3 / (5 * Pi)) * Integral_{E_rel}^{Infty} K_{5/3}(x) dx
! Where
!   P(E_rel)) = Probability of finding a photon at relative energy E_rel.
!   E_rel     = Relative photon energy: E / E_crit, E_crit = Critical energy.
!   K_{5/3}   = Modified Bessel function.
!
! Notice that the P(E) is not the same as the distribution radiation energy since
! the photons must be energy weighted.
!
! There is a cut-off built into the calculation so that E_rel will be in the 
! range [0, 31.4]. The error in neglecting photons with E_rel > 31.4 translates
! to neglecting one photon for every 10^15 generated.
! If r_in is present: 
!   r_in = 0 => E_rel = 0 
!   r_in = 1 => E_rel = 31.4
!
! Input:
!   r_in  -- Real(rp), optional: Integrated probability in the range [0,1].
!             If not present, a random number will be used.
!
! Output:
!   E_rel -- Real(rp): Relative photon energy E/E_crit. 
!-

function bend_photon_e_rel_init (r_in) result (E_rel)

implicit none

! Four spline fit arrays are used. 
! Each fit array has a different del_x and range of validity.

integer, parameter :: n_spline = 11
type (photon_init_spline_struct), save :: spline(n_spline) 

real(rp) E_rel
real(rp), optional :: r_in
real(rp) rr, rr1, dr, alpha, c0

integer i, n, is

logical, save :: init_needed = .true.
character(*), parameter :: r_name = 'bend_photon_e_rel_init'

! Check for r_in

if (present(r_in)) then
  rr = r_in
  if (rr < 0  .or. rr > 1) then
    call out_io (s_fatal$, r_name, 'R_IN IS OUT OF RANGE: \es12.4\ ', rr)
    stop
  endif
else
  call ran_uniform(rr)
endif

! Init

if (init_needed) then
  call init_this()
  init_needed = .false.
endif


! Spline fit

do is = 1, n_spline
  if (rr > spline(is)%x_max) cycle
  E_rel = photon_init_spline_eval (spline(is), rr, .true.)

  ! Spline(1) is spline of E_rel(rr) / rr^3
  if (is == 1) then
    E_rel = E_rel * rr**3
  ! Spline(2:) are splines of E_rel(rr) * (1 - rr + 24(1-rr)^2)
  else
    rr1 = 1.0_rp - rr
    E_rel = E_rel / (rr1 + 24.0_rp * rr1**2)
  endif

  return
enddo

! In the range above the last spline use the approximation that 1 - P_integrated(E_rel) ~ e^-E_rel / Sqrt[E_rel]

rr1 = 1 - rr
if (rr1 < 1d-15) then
  E_rel = 31.4_rp
  return
endif

n = ubound(spline(n_spline)%pt, 1)
dr = (1.0_rp - spline(n_spline)%x_max)
c0 = spline(n_spline)%pt(n)%c0 / (dr + 24.0_rp * dr**2)
alpha = dr * sqrt(c0) * exp(c0) 
E_rel = inverse(p_func, rr1, c0, 50.0_rp, 1e-10_rp) 

!---------------------------------------------------------------------------
contains

function p_func(E_in) result(rr1)
real(rp) :: E_in, rr1
rr1 = alpha * exp(-E_in) / sqrt(E_in)
end function p_func

!-----------------------------------------------
! contains

subroutine init_this()

integer is

! Init. 
! The values for c0 and c1 were obtained from a Mathematica calculation See:
!   bmad/calculations/bend_radiation_distribution_spline.nb

spline(:)%spline_type = gen_poly_spline$
spline(1:2)%del_x = [0.02_rp, 0.01_rp]
spline(1:2)%x_min = [0.00_rp, 0.80_rp]
spline(1:2)%x_max = [0.80_rp, 0.90_rp]

do i = 3, n_spline
  spline(i)%del_x = 10.0_rp**(1-i)
  spline(i)%x_min = 1.0_rp - 10.0_rp * spline(i)%del_x
  spline(i)%x_max = 1.0_rp - spline(i)%del_x
  allocate(spline(i)%pt(0:9))
enddo

allocate (spline(1)%pt(0:40), spline(2)%pt(0:10))

spline(1)%pt(:)%c0 = [0.53530664D0, 5.3542580D-1, 5.3578369D-1, 5.3638158D-1, 5.3722163D-1, 5.3830686D-1, 5.3964123D-1, 5.4122962D-1, 5.4307796D-1, 5.4519321D-1, &
                      5.4758347D-1, 5.5025808D-1, 5.5322768D-1, 5.5650439D-1, 5.6010187D-1, 5.6403557D-1, 5.6832287D-1, 5.7298331D-1, 5.7803888D-1, 5.8351430D-1, &
                      5.8943743D-1, 5.9583969D-1, 6.0275655D-1, 6.1022821D-1, 6.1830033D-1, 6.2702493D-1, 6.3646151D-1, 6.4667842D-1, 6.5775452D-1, 6.6978129D-1, &
                      6.8286551D-1, 6.9713253D-1, 7.1273070D-1, 7.2983690D-1, 7.4866399D-1, 7.6947077D-1, 7.9257547D-1, 8.1837460D-1, 8.4736963D-1, 8.8020587D-1, 9.1773074D-1]
spline(1)%pt(:)%c1 = [0.0000000D0,  1.1918866D-2, 2.3880242D-2, 3.5927047D-2, 4.8103023D-2, 6.0453156D-2, 7.3024142D-2, 8.5864864D-2, 9.9026927D-2, 1.1256524D-1, &
                      1.2653865D-1, 1.4101073D-1, 1.5605055D-1, 1.7173369D-1, 1.8814331D-1, 2.0537152D-1, 2.2352083D-1, 2.4270602D-1, 2.6305629D-1, 2.8471780D-1, &
                      3.0785685D-1, 3.3266361D-1, 3.5935678D-1, 3.8818930D-1, 4.1945539D-1, 4.5349939D-1, 4.9072695D-1, 5.3161914D-1, 5.7675074D-1, 6.2681380D-1, &
                      6.8264874D-1, 7.4528554D-1, 8.1599929D-1, 8.9638612D-1, 9.8846871D-1, 1.0948456D0, 1.2189070D0, 1.3651535D0, 1.5396801D0, 1.7509330D0, 2.0109362D0]

spline(2)%pt(:)%c0 = [5.4505864D-1, 5.2693936D-1, 5.0744281D-1, 4.8659418D-1, 4.6442647D-1, 4.4098127D-1, 4.1630981D-1, 3.9047405D-1, 3.6354797D-1, 3.3561917D-1, 3.0679058D-1]
spline(2)%pt(:)%c1 = [-1.7424031D0, -1.8811506D0, -2.0177413D0, -2.1514358D0, -2.2814157D0, -2.4067717D0, -2.5264894D0, -2.6394331D0, -2.7443254D0, -2.8397233D0, -2.9239908D0]

spline(3)%pt(:)%c0 = [3.0679058D-1, 2.7718263D-1, 2.4693555D-1, 2.1621211D-1, 1.8520027D-1, 1.5411505D-1, 1.2319719D-1, 9.2700067D-2, 6.2832118D-2, 3.3473224D-2]
spline(3)%pt(:)%c1 = [-2.9239908D0, -2.9952677D0, -3.0514396D0, -3.0901183D0, -3.1086655D0, -3.1043497D0, -3.0749282D0, -3.0207144D0, -2.9530530D0, -2.9453773D0]

spline(4)%pt(:)%c0 = [3.3473224D-2, 3.0521288D-2, 2.7552183D-2, 2.4558263D-2, 2.1528977D-2, 1.8449250D-2, 1.5296480D-2, 1.2034274D-2, 8.5969339D-3, 4.8371545D-3]
spline(4)%pt(:)%c1 = [-2.9453773D0, -2.9594339D0, -2.9800353D0, -3.0095372D0, -3.0515015D0, -3.1116260D0, -3.1997606D0, -3.3349614D0, -3.5614189D0, -4.0251877D0]

spline(5)%pt(:)%c0 = [4.8371545D-3, 4.4308560D-3, 4.0163160D-3, 3.5923103D-3, 3.1572513D-3, 2.7090009D-3, 2.2445295D-3, 1.7592180D-3, 1.2451594D-3, 6.8554971D-4]
spline(5)%pt(:)%c1 = [-4.0251877D0, -4.1023912D0, -4.1904312D0, -4.2922987D0, -4.4123925D0, -4.5575601D0, -4.7393250D0, -4.9793935D0, -5.3267615D0, -5.9388177D0]

spline(6)%pt(:)%c0 = [6.8554971D-4, 6.2569719D-4, 5.6484320D-4, 5.0285761D-4, 4.3957247D-4, 3.7476284D-4, 3.0811164D-4, 2.3913803D-4, 1.6702372D-4, 9.0040425D-5]
spline(6)%pt(:)%c1 = [-5.9388177D0, -6.0334049D0, -6.1395428D0, -6.2603451D0, -6.4003771D0, -6.5667262D0, -6.7712821D0, -7.0363711D0, -7.4122167D0, -8.0595166D0]

spline(7)%pt(:)%c0 = [9.0040425D-5, 8.1932372D-5, 7.3719822D-5, 6.5389516D-5, 5.6924365D-5, 4.8301489D-5, 3.9488672D-5, 3.0437089D-5, 2.1063787D-5, 1.1193742D-5]
spline(7)%pt(:)%c1 = [-8.0595166D0, -8.1583455D0, -8.2689461D0, -8.3944800D0, -8.5395802D0, -8.7114339D0, -8.9220937D0, -9.1941712D0, -9.5784946D0, -1.0237525D1]

spline(8)%pt(:)%c0 = [1.1193746D-5, 1.0165064D-5, 9.1257720D-6, 8.0745305D-6, 7.0096135D-6, 5.9287108D-6, 4.8285707D-6, 3.7042684D-6, 2.5474352D-6, 1.3404298D-6]
spline(8)%pt(:)%c1 = [-1.0237520D1, -1.0337895D1, -1.0450162D1, -1.0577514D1, -1.0724623D1, -1.0898740D1, -1.1112019D1, -1.1387257D1, -1.1775677D1, -1.2440922D1]

spline(9)%pt(:)%c0 = [1.3404298D-6, 1.2155233D-6, 1.0895466D-6, 9.6236497D-7, 8.3380457D-7, 7.0363288D-7, 5.7152261D-7, 4.3697821D-7, 2.9915778D-7, 1.5628710D-7]
spline(9)%pt(:)%c1 = [-1.2440922D1, -1.2542169D1, -1.2655392D1, -1.2783801D1, -1.2932099D1, -1.3107583D1, -1.3322477D1, -1.3599712D1, -1.3990797D1, -1.4660240D1]

spline(10)%pt(:)%c0 = [1.5628710D-7, 1.4157683D-7, 1.2675892D-7, 1.1181981D-7, 9.6742042D-8, 8.1502243D-8, 6.6067533D-8, 5.0388124D-8, 3.4379451D-8, 1.7863331D-8]
spline(10)%pt(:)%c1 = [-1.4660240D1, -1.4762089D1, -1.4875975D1, -1.5005123D1, -1.5154258D1, -1.5330709D1, -1.5546756D1, -1.5825429D1, -1.6218455D1, -1.6891000D1]

spline(11)%pt(:)%c0 = [1.7863331D-8, 1.6169206D-8, 1.4464270D-8, 1.2747161D-8, 1.1016127D-8, 9.2688230D-9, 7.5019478D-9, 5.7105054D-9, 3.8860096D-9, 2.0105829D-9]
spline(11)%pt(:)%c1 = [-1.6891000D1, -1.6993298D1, -1.7107680D1, -1.7237382D1, -1.7387147D1, -1.7564328D1, -1.7781248D1, -1.8061016D1, -1.8455527D1, -1.9130467D1]

! Fill in rest of the spline fit coefs.

do is = 1, ubound(spline, 1)
  call photon_init_spline_coef_calc (spline(is))
enddo

end subroutine init_this

end function bend_photon_e_rel_init

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!+
! Function bend_photon_energy_normalized_probability (E_rel) result (prob)
!
! Routine to return the normalized probability that a photon will be emitted in a bend with energy
! E_rel relative to the critical energy. The probability is normalized such that
!   Integral[0,Infinity] dE_rel P(E_rel) = 1
!
! Input:
!   E_rel     -- real(rp): Photon energy relative to the critical energy.
!
! Output:
!   prob      -- real(rp): Normalized probability.
!-

function bend_photon_energy_normalized_probability (E_rel) result (prob)

implicit none

type (photon_init_spline_struct), save :: p(6)

real(rp) E_rel, prob
real(rp) E0, log_E

integer i, n

logical, save :: init_needed = .true.

!

if (init_needed) then
  call init_this()
  init_needed = .false.
endif

if (E_rel < 1.0e-160_rp) then
  prob = 1e100_rp ! In theory prob -> Infinity at E_rel = 0
  return
endif

! p(1) is a spline of E_rel^(2/3) * prob vs log(E_rel)

if (E_rel < p(2)%x_min) then
  log_E = log10(E_rel)
  if (log_E > p(1)%x_min) then
    prob = photon_init_spline_eval(p(1), log_E) * E_rel**(-2.0_rp/3.0_rp)
  else  ! Extrapolate: prob ~ E_rel^(-2/3)
    prob = p(1)%pt(0)%c0 * (E_rel / 10**p(1)%x_min)**(-2.0_rp/3.0_rp)
  endif
  return
endif

! p(2) and p(3) is a spline of E_rel^(2/3) * prob vs E_rel

do i = 2, 3
  if (E_rel > p(i)%x_max) cycle
  prob = photon_init_spline_eval(p(i), E_rel) * E_rel**(-2.0_rp/3.0_rp)
  return
enddo

! p(4:6) is a spline of Exp[E_rel] * Sqrt[E_rel] * prob vs E_rel

do i = 4, 6
  if (E_rel > p(i)%x_max) cycle
  prob = photon_init_spline_eval(p(i), E_rel) / (exp(E_rel) * sqrt(E_rel))
  return
enddo

! Extrapolate past p(6)

n = size(p(6)%pt)
E0 = p(6)%x_max
prob = p(6)%pt(n)%c0 * exp(-E_rel) / sqrt(E_rel)

!--------------------------------
contains

subroutine init_this()

integer i

!

p%spline_type = gen_poly_spline$
p%x_min = [-7.0_rp, 0.01_rp, 0.1_rp,  1.0_rp,  3.0_rp,  10.0_rp]
p%x_max = [-2.0_rp, 0.1_rp,  1.0_rp,  3.0_rp,  10.0_rp, 40.0_rp]
p%del_x = [ 1.0_rp, 0.01_rp, 0.02_rp, 0.05_rp, 0.2_rp,  1.0_rp]

allocate(p(1)%pt(0:5))
p(1)%pt%c0 = [4.1052219D-1, 4.1049501D-1, 4.1036886D-1, 4.0978334D-1, 4.0706563D-1, 3.9445841D-1]
p(1)%pt%c1 = [-1.1456421D-5, -5.3175995D-5, -2.4682107D-4, -1.1456386D-3, -5.3172381D-3, -2.4646674D-2]

allocate(p(2)%pt(0:9))
p(2)%pt%c0 = [3.9445841D-1, 3.8503667D-1, 3.7715318D-1, 3.7013587D-1, 3.6370582D-1, 3.5771314D-1, 3.5206552D-1, 3.4670102D-1, 3.4157550D-1, 3.3665612D-1]
p(2)%pt%c1 = [-1.0703915D0, -8.4772113D-1, -7.3863972D-1, -6.6916168D-1, -6.1924820D-1, -5.8078618D-1, -5.4975032D-1, -5.2387509D-1, -5.0176995D-1, -4.8252543D-1]

allocate(p(3)%pt(0:45))
p(3)%pt%c0 = [3.3665612D-1, 3.2733982D-1, 3.1860491D-1, 3.1035234D-1, 3.0251133D-1, 2.9502902D-1, 2.8786450D-1, 2.8098528D-1, 2.7436495D-1, 2.6798161D-1, &
              2.6181686D-1, 2.5585499D-1, 2.5008243D-1, 2.4448736D-1, 2.3905935D-1, 2.3378917D-1, 2.2866855D-1, 2.2369004D-1, 2.1884691D-1, 2.1413302D-1, &
              2.0954277D-1, 2.0507099D-1, 2.0071291D-1, 1.9646414D-1, 1.9232056D-1, 1.8827836D-1, 1.8433395D-1, 1.8048397D-1, 1.7672528D-1, 1.7305489D-1, &
              1.6946999D-1, 1.6596792D-1, 1.6254615D-1, 1.5920228D-1, 1.5593402D-1, 1.5273919D-1, 1.4961572D-1, 1.4656161D-1, 1.4357497D-1, 1.4065395D-1, &
              1.3779681D-1, 1.3500187D-1, 1.3226751D-1, 1.2959217D-1, 1.2697435D-1, 1.2441259D-1]
p(3)%pt%c1 = [-4.8252543D-1, -4.5029908D-1, -4.2400202D-1, -4.0183669D-1, -3.8269981D-1, -3.6586975D-1, -3.5085095D-1, -3.3728967D-1, -3.2492535D-1, -3.1356101D-1, &
              -3.0304439D-1, -2.9325558D-1, -2.8409857D-1, -2.7549540D-1, -2.6738195D-1, -2.5970488D-1, -2.5241938D-1, -2.4548744D-1, -2.3887657D-1, -2.3255877D-1, &
              -2.2650975D-1, -2.2070831D-1, -2.1513581D-1, -2.0977581D-1, -2.0461371D-1, -1.9963648D-1, -1.9483244D-1, -1.9019109D-1, -1.8570296D-1, -1.8135944D-1, &
              -1.7715270D-1, -1.7307560D-1, -1.6912160D-1, -1.6528469D-1, -1.6155934D-1, -1.5794043D-1, -1.5442323D-1, -1.5100336D-1, -1.4767671D-1, -1.4443949D-1, &
              -1.4128812D-1, -1.3821929D-1, -1.3522985D-1, -1.3231687D-1, -1.2947759D-1, -1.2670939D-1]

allocate(p(4)%pt(0:40))
p(4)%pt%c0 = [3.3818849D-1, 3.3516487D-1, 3.3233973D-1, 3.2969287D-1, 3.2720686D-1, 3.2486656D-1, 3.2265876D-1, 3.2057185D-1, 3.1859560D-1, 3.1672095D-1, &
              3.1493984D-1, 3.1324506D-1, 3.1163017D-1, 3.1008936D-1, 3.0861740D-1, 3.0720955D-1, 3.0586152D-1, 3.0456940D-1, 3.0332963D-1, 3.0213895D-1, &
              3.0099437D-1, 2.9989316D-1, 2.9883280D-1, 2.9781096D-1, 2.9682551D-1, 2.9587445D-1, 2.9495594D-1, 2.9406829D-1, 2.9320991D-1, 2.9237932D-1, &
              2.9157514D-1, 2.9079610D-1, 2.9004099D-1, 2.8930869D-1, 2.8859815D-1, 2.8790838D-1, 2.8723846D-1, 2.8658751D-1, 2.8595472D-1, 2.8533932D-1, 2.8474057D-1]
p(4)%pt%c1 = [-6.2608102D-2, -5.8415126D-2, -5.4657613D-2, -5.1274686D-2, -4.8216055D-2, -4.5439837D-2, -4.2910890D-2, -4.0599528D-2, -3.8480510D-2, -3.6532254D-2, &
              -3.4736203D-2, -3.3076320D-2, -3.1538682D-2, -3.0111146D-2, -2.8783081D-2, -2.7545137D-2, -2.6389069D-2, -2.5307573D-2, -2.4294163D-2, -2.3343061D-2, &
              -2.2449100D-2, -2.1607654D-2, -2.0814563D-2, -2.0066084D-2, -1.9358833D-2, -1.8689751D-2, -1.8056063D-2, -1.7455246D-2, -1.6885003D-2, -1.6343240D-2, &
              -1.5828040D-2, -1.5337649D-2, -1.4870460D-2, -1.4424995D-2, -1.3999893D-2, -1.3593903D-2, -1.3205867D-2, -1.2834719D-2, -1.2479468D-2, -1.2139199D-2, -1.1813062D-2]

allocate(p(5)%pt(0:35))
p(5)%pt%c0 = [2.8474057D-1, 2.8249896D-1, 2.8047519D-1, 2.7863827D-1, 2.7696291D-1, 2.7542827D-1, 2.7401700D-1, 2.7271453D-1, 2.7150852D-1, 2.7038847D-1, &
              2.6934536D-1, 2.6837140D-1, 2.6745984D-1, 2.6660478D-1, 2.6580106D-1, 2.6504414D-1, 2.6432999D-1, 2.6365504D-1, 2.6301612D-1, 2.6241038D-1, &
              2.6183528D-1, 2.6128852D-1, 2.6076805D-1, 2.6027198D-1, 2.5979864D-1, 2.5934648D-1, 2.5891409D-1, 2.5850020D-1, 2.5810363D-1, 2.5772332D-1, &
              2.5735827D-1, 2.5700758D-1, 2.5667040D-1, 2.5634598D-1, 2.5603358D-1, 2.5573256D-1]
p(5)%pt%c1 = [-1.1813062D-2, -1.0634850D-2, -9.6285326D-3, -8.7616650D-3, -8.0092079D-3, -7.3515685D-3, -6.7732249D-3, -6.2617422D-3, -5.8070577D-3, -5.4009530D-3, &
              -5.0366595D-3, -4.7085591D-3, -4.4119559D-3, -4.1428985D-3, -3.8980419D-3, -3.6745382D-3, -3.4699500D-3, -3.2821809D-3, -3.1094194D-3, -2.9500932D-3, &
              -2.8028323D-3, -2.6664382D-3, -2.5398586D-3, -2.4221666D-3, -2.3125427D-3, -2.2102605D-3, -2.1146740D-3, -2.0252072D-3, -1.9413450D-3, -1.8626257D-3, &
              -1.7886346D-3, -1.7189981D-3, -1.6533789D-3, -1.5914720D-3, -1.5330010D-3, -1.4777149D-3]

allocate(p(6)%pt(0:30))
p(6)%pt%c0 = [2.5573256D-1, 2.5437801D-1, 2.5323167D-1, 2.5224879D-1, 2.5139659D-1, 2.5065056D-1, 2.4999195D-1, 2.4940622D-1, 2.4888186D-1, 2.4840969D-1, &
              2.4798226D-1, 2.4759351D-1, 2.4723838D-1, 2.4691270D-1, 2.4661294D-1, 2.4633612D-1, 2.4607970D-1, 2.4584150D-1, 2.4561965D-1, 2.4541251D-1, &
              2.4521867D-1, 2.4503689D-1, 2.4486607D-1, 2.4470524D-1, 2.4455357D-1, 2.4441028D-1, 2.4427470D-1, 2.4414622D-1, 2.4402430D-1, 2.4390845D-1, 2.4379822D-1]
p(6)%pt%c1 = [-1.4777149D-3, -1.2417257D-3, -1.0583022D-3, -9.1286524D-4, -7.9557438D-4, -6.9958881D-4, -6.2003118D-4, -5.5334694D-4, -4.9689529D-4, -4.4868069D-4, &
              -4.0717218D-4, -3.7117922D-4, -3.3976460D-4, -3.1218231D-4, -2.8783255D-4, -2.6622856D-4, -2.4697197D-4, -2.2973422D-4, -2.1424234D-4, -2.0026803D-4, &
              -1.8761915D-4, -1.7613308D-4, -1.6567140D-4, -1.5611572D-4, -1.4736428D-4, -1.3932922D-4, -1.3193437D-4, -1.2511345D-4, -1.1880850D-4, -1.1296875D-4, -1.0754950D-4]

do i = 1, size(p)
  call photon_init_spline_coef_calc(p(i))
enddo

end subroutine init_this

end function bend_photon_energy_normalized_probability

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!+
! Subroutine photon_init_spline_coef_calc (spline)
!
! Routine to calculate spline coefs c2 and c3 given coefs c0 and c1.
! Note: coef c1 is scaled.
! 
! Input:
!   spline -- photon_init_spline_struct: Spline with c0 and c1 coefs set.
!
! Output:
!   spline -- photon_init_spline_struct: Spline with c2 and c3 coefs calculated.
!-

subroutine photon_init_spline_coef_calc (spline)

implicit none

type (photon_init_spline_struct) spline
integer i, ns
real(rp) v, vp

!

select case (spline%spline_type)

case (gen_poly_spline$)
  spline%pt(:)%c1 = spline%pt(:)%c1 * spline%del_x
  ns = ubound(spline%pt, 1) 

  do i = 0, ns-1
    v  = spline%pt(i+1)%c0 - spline%pt(i)%c0 - spline%pt(i)%c1
    vp = spline%pt(i+1)%c1 - spline%pt(i)%c1
    if (v * vp > 0 .and. abs(vp) > abs(3 * v)) then   ! c0 + c1 x + c2 * x^c3 spline
      spline%pt(i)%c2 = v
      spline%pt(i)%c3 = vp / v
    else    ! Cubic spline
      spline%pt(i)%c2 = 3 * v - vp 
      spline%pt(i)%c3 = vp - 2 * v
    endif
  enddo

  spline%pt(ns)%c2 = 0  ! Need to set this due to roundoff errors
  spline%pt(ns)%c3 = 0  ! Need to set this due to roundoff errors

end select

end subroutine photon_init_spline_coef_calc

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!+
! Function photon_init_spline_eval (spline, rr, extrapolate) result (fit_val)
!
! Function to evaluate a spline fit at rr.
! 
! Input:
!   spline      -- photon_init_spline_struct: spline section. spline%pt(:) must start at 0.
!   rr          -- real(rp): Value to evaluate the fit at.
!   extrapolate -- logical, optional: If False (the default) and rr is not in the x-range of the spline,
!                   the value of rr used in the calculation will be the x-value at the nearest edge of 
!                   the spline range. That is, if extrapolate = F and rr is not in range, fit_val will be 
!                   set to the value at the edge of the spline range.
!
! Output:
!   fit_val -- real(rp): Spline fit evaluated at rr.
!-

function photon_init_spline_eval (spline, rr, extrapolate) result (fit_val)

implicit none

type (photon_init_spline_struct) spline
real(rp) rr, r_rel, x, fit_val, v, vp
integer i, ix, np
logical, optional :: extrapolate

!  Programming error...

if (lbound(spline%pt, 1) /= 0) call err_exit

! Find in which spline section rr is in

r_rel = (rr - spline%x_min) / spline%del_x
i = int(r_rel)   ! Index of which spline section to use.
x = r_rel - i    ! Notice that x will be in the range [0, 1].

! If out of range then adjust point to be at the edge of the range

if (i < 0) then
  i = 0
  if (logic_option(.false., extrapolate)) then
    x = r_rel
  else
    x = 0
  endif
endif

np = ubound(spline%pt, 1)
if (i > np - 1) then 
  i = np - 1
  if (logic_option(.false., extrapolate)) then
    x = r_rel - i
  else
    x = 1
  endif
endif

select case (spline%spline_type)

case (gen_poly_spline$)

  if (i == ubound(spline%pt, 1)) then
    fit_val = spline%pt(i)%c0
    return
  endif

  v  = spline%pt(i+1)%c0 - spline%pt(i)%c0 - spline%pt(i)%c1
  vp = spline%pt(i+1)%c1 - spline%pt(i)%c1

  if (v * vp > 0 .and. abs(vp) > abs(3 * v)) then   ! c0 + c1 x + c2 * x^c3 spline
    fit_val = spline%pt(i)%c0 + spline%pt(i)%c1 * x + spline%pt(i)%c2 * x**spline%pt(i)%c3
  else    ! Cubic spline
    fit_val = spline%pt(i)%c0 + spline%pt(i)%c1 * x + spline%pt(i)%c2 * x**2 + spline%pt(i)%c3 * x**3
  endif

case (end_spline$)

  fit_val = spline%pt(i)%c0 + spline%pt(i)%c1 * x + spline%pt(i)%c2 * x**2 / (1 - x/spline%pt(i)%c3)

! Coding error if here.

case default
  if (global_com%exit_on_error) call err_exit

end select

end function photon_init_spline_eval

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!+
! Function E_crit_photon (gamma, g_bend) result (E_crit)
!
! Routine to calculate the photon critical energy in a bend.
!
! Input:
!   gamma   -- real(rp): Gamma factor of charged particle emitting photon.
!   g_bend  -- real(rp): 1/radius bending strength.
!
! Output:
!   E_crit  -- real(rp): Critical photon energy.
!-

function E_crit_photon (gamma, g_bend) result (E_crit)

real(rp) gamma, g_bend, E_crit
real(rp), parameter :: e_factor = 3.0_rp * h_bar_planck * c_light / 2.0_rp

!

E_crit = e_factor * gamma**3 * g_bend

end function E_crit_photon

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Function init_photon_integ_prob(gamma, g, E_min, E_max, vert_angle_min, 
!              vert_angle_max, vert_angle_symmetric, energy_integ_prob, E_photon) result (integ_prob)
!
! Routine to calcuate the integrated probability of emitting a photon in a given vertical angle range
! and in a given energy range
!
! Input:
!   gamma                 -- real(rp): Gamma factor of charged particle emitting photon.
!   g                     -- real(rp): 1/rho bending strength.
!   E_min                 -- real(rp): Minimum photon energy.
!   E_max                 -- real(rp): Maximum photon energy.
!   vert_angle_min        -- real(rp), optional: Lower bound of vertical angle range.
!   vert_angle_max        -- real(rp), optional: Upper bound of vertical angle range.
!   vert_angle_symmetric  -- real(rp), optional: Use two symmetric ranges [-vert_angle_max, -vert_angle_min] and 
!                             [vert_angle_min, vert_angle_max] instead of just [vert_angle_min, vert_angle_max]?
!   energy_integ_prob     -- real(rp), optional: If present, E_photon will be set to the photon energy such that
!                             the integrated probability of generating a photon in the given angle and energy
!                             range in the interval [E_min, E_photon] is energy_integ_prob. That is,
!                             energy_integ_prob = 0 => E_photon = E_min and energy_integ_prob = 1 => E_photon = E_max.
!
! Output:
!   E_photon              -- real(rp): See energy_integ_prob. E_photon must be present if energy_integ_prob is.
!   integ_prob            -- real(rp): Integrated probablility of emitting a photon in given angle and energy range.
!-

function init_photon_integ_prob (gamma, g, E_min, E_max, vert_angle_min, &
                   vert_angle_max, vert_angle_symmetric, energy_integ_prob, E_photon) result (integ_prob)

type adaptive_integration_interval_struct
  real(rp) x(0:4)         ! [x(0), x(4)] = Interval range
  real(rp) y(0:4)         ! Function value at points x(0:4)
  real(rp) Integral       ! Integral over interval
  real(rp) int_err        ! Integral estimated error
  integer :: ix_previous = -1 ! Next interval
  integer :: ix_next = -1     ! previous interval
end type

type (adaptive_integration_interval_struct), target :: interval(2000)
type (adaptive_integration_interval_struct), pointer :: inter, inter2


real(rp) gamma, g, E_min, E_max, integ_prob
real(rp), optional :: vert_angle_min, vert_angle_max, energy_integ_prob, E_photon
real(rp) r_energy_min, r_energy_max, v_angle_min, v_angle_max, int_tot, r, ai, dx
real(rp) i1, i2, in_prob, r_target, r_left, eps, a, b, c
real(rp) :: rel_tol = 1e-6_rp, abs_tol = 1e-12_rp

integer i, n_int, ix_stack(2000), n_stack, ixs

logical, optional :: vert_angle_symmetric
logical err_flag

character(*), parameter :: r_name = 'init_photon_integ_prob'

!

v_angle_min = real_option(-pi/2, vert_angle_min)
v_angle_max = real_option(pi/2, vert_angle_max)

r_energy_min = 0
r_energy_max = 1

if (E_min > 0) r_energy_min = bend_photon_energy_integ_prob(E_min, g, gamma)
if (E_max > 0) r_energy_max = bend_photon_energy_integ_prob(E_max, g, gamma)

! Simple case where all angles are within range.

if (v_angle_min <= -pi/2 .and. v_angle_max >= pi/2) then
  integ_prob = r_energy_max - r_energy_min
  if (present(energy_integ_prob)) then
    r = r_energy_min + energy_integ_prob * (r_energy_max - r_energy_min)
    E_photon = bend_photon_e_rel_init(r) * E_crit_photon(gamma, g)
  endif
  return
endif

! Use adaptive quadrature to compute the integrated probability.
! The stack is a list of intervals where the computed integral is not accurate enough.
! Initally, the stack has one interval that spans the range to be integrated over.

n_int = 1
n_stack = 1
ix_stack(1) = 1
inter => interval(1)
do i = 0, 4, 2
  ai = real(i, rp)
  inter%x(i) = ((4.0_rp - ai) * r_energy_min + ai * r_energy_max) / 4.0_rp
  inter%y(i) = vert_prob_func(inter%x(i))
enddo

integ_prob = 0

do
  if (n_stack == 0) exit
  ixs = ix_stack(n_stack)
  inter => interval(ixs)

  do i = 1, 3, 2
    ai = real(i, rp)
    inter%x(i) = ((4.0_rp - ai) * inter%x(0) + ai * inter%x(4)) / 4.0_rp
    inter%y(i) = vert_prob_func(inter%x(i))
  enddo

  dx  = inter%x(4) - inter%x(0)
  I1 = dx * (inter%y(0) + 4.0_rp * inter%y(2) + inter%y(4)) / 6.0_rp
  I2 = dx * (inter%y(0) + 4.0_rp * inter%y(1) + 2.0_rp * inter%y(2) + 4.0_rp * inter%y(3) + inter%y(4)) / 12.0_rp
  inter%integral = i2
  inter%int_err = abs(i2 - i1) / 15.0_rp

  in_prob = integ_prob + interval(1)%integral
  ! If integral error is too large then split interval and put both halves on the stack. 
  if (inter%int_err > dx * (in_prob * rel_tol + abs_tol)) then
    n_int = n_int + 1
    if (n_int > size(interval)) then
      call out_io (s_error$, r_name, 'PROBABILITY INTEGRATION FAILED! PLEASE GET HELP!')
      exit
    endif
    inter2 => interval(n_int)
    inter2%x(0:4:2) = inter%x(2:4)
    inter2%y(0:4:2) = inter%y(2:4)
    inter2%ix_previous = ixs
    inter2%ix_next = inter%ix_next

    inter%x(0:4:2) = inter%x(1:3)
    inter%y(0:4:2) = inter%y(1:3)
    inter%ix_next = n_int

    n_stack = n_stack + 1
    ix_stack(n_stack) = n_int

  ! else integral error is acceptable and so remove this interval from the stack
  else
    integ_prob = integ_prob + inter%integral
    n_stack = n_stack - 1
  endif

enddo  

! If we need to compute E_photon then keep dividing until the desired accuracy is reached.

if (present(energy_integ_prob)) then
  r_target = energy_integ_prob * integ_prob
  r_left = 0
  ixs = 1
  do
    inter => interval(ixs)
    r_left = r_left + inter%integral
    if (r_left >= r_target) exit
    ixs = inter%ix_next
  enddo

  r_left = r_left - inter%integral ! Integral up to left edge of inverval
  dx = inter%x(4) - inter%x(0)
  eps = (inter%integral / dx - (inter%y(0) + inter%y(4)) / 2.0_rp)
  a = (inter%y(4) - inter%y(0)) / (2.0_rp * dx)
  b = inter%y(0) + eps  ! eps is a fudge factor since the true function is not exactly quadratic.
  c = -(r_target - r_left)
  r = inter%x(0) - 2.0_rp * c / (b + sqrt(b*b - 4.0_rp * a * c))  ! Quadriatic formula
  E_photon = bend_photon_e_rel_init(r) * E_crit_photon(gamma, g)
endif

if (logic_option(.false., vert_angle_symmetric)) integ_prob = 2 * integ_prob

!-------------------
contains

function vert_prob_func(r_energy) result (prob)

real(rp) :: r_energy, prob
real(rp) E_rel, r_min, r_max

!

E_rel = bend_photon_e_rel_init(r_energy)
r_min = bend_vert_angle_integ_prob(v_angle_min, E_rel, gamma)
r_max = bend_vert_angle_integ_prob(v_angle_max, E_rel, gamma)

prob = r_max - r_min

end function vert_prob_func

end function init_photon_integ_prob

end module
