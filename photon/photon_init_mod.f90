module photon_init_mod

use bmad_struct
use bmad_interface
use random_mod

! An init_spectrum_struct holds an array of spline fits of E_rel of gamma_phi vs 
! integrated probability r over a certain range. Each spline section is fit to 
! one of four forms. 
! The poly_n_spline$ is of the form:
!   fit_value = c0 + c1 * t + c2 * t^c3
! or
!   fit_value = c0 + c1 * t + c2 * t^2 + c3 * t^3
! where t = a + b * r is in the range [0, 1].
!   t = 0 at the start of the spline section.
!   t = 1 at the end of the spline section.
! The end_spline$ is of the form:
!   fit_value = c0 + c1 * t + c2 * t^2 / (1 - t/c3)

type photon_init_spline_pt_struct
  real(rp) c0, c1, c2, c3
end type

type photon_init_spline_struct
  real(rp) del_x    ! Spacing between spline points
  real(rp) x_min    ! Lower bound
  real(rp) x_max    ! Upper bound of Region of validity of this spline fit.
                    ! The lower bound is given by the upper bound of the previos struct.
  type (photon_init_spline_pt_struct), allocatable :: pt(:)
  integer spline_type
end type

type photon_init_spline2_struct
  real(rp) del_y    ! Spacing between spline arrayss
  real(rp) y_min    ! Lower bound
  real(rp) y_max    ! Upper bound
  type (photon_init_spline_struct), allocatable :: int_prob(:)
end type

integer, parameter :: gen_poly_spline$ = 1, end_spline$ = 2

private photon_init_spline_fit, photon_init_spline_coef_calc

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
!   orbit(4) -- Real(rp): Vertical angle. Limited to be less than 0.1
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
if (abs(orbit(4)) > 0.1) orbit(4) = sign(0.1_rp, orbit(4))

orbit(6) = E_rel * 3 * h_bar_planck * c_light * gamma**3 * g_bend / 2

end subroutine photon_init

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!+
! Subroutine photon_vert_angle_init (E_rel, gamma_phi, r_in)
!
! Routine to convert a "random" number in the interval [0,1] to a photon vertical emission 
! angle for a simple bend.
!
! Module needed:
!   use photon_init_mod
!
! Input:
!   r_in  -- Real(rp), optional: number in the range [0,1].
!             If not present, a random number will be used.
!             Note: gamma_phi is an increasing monotonic function of r_in except at one point
!             where it is discontinuous.
!   E_rel -- Real(rp): Relative photon energy E/E_crit. 
! 
! Output:
!   gamma_phi -- Real(rp): gamma * phi where gamma is the beam relativistic factor and
!                 phi is the vertical photon angle (in radians).
!-

subroutine photon_vert_angle_init (E_rel, gamma_phi, r_in)

implicit none

type (photon_init_spline2_struct), target, save :: p(2), dp(2)
type (photon_init_spline_struct), save :: spline

real(rp), optional :: r_in
real(rp) e_rel, gamma_phi, p_perp, sig, log_E_rel, x0, xp0
real(rp) rr, r, ss, x, log_E, frac, rro, drr, log_E_min, log_E_max, del_log_E
real(rp) x1, xp1, v, vp

real(rp) :: rel_amp(0:32) = (/ &
                  0.249546, 0.249383, 0.249162, 0.248862, 0.248454, &
                  0.247902, 0.247154, 0.246141, 0.244774, 0.242931, &
                  0.240455, 0.237143, 0.232737, 0.226921, 0.219325, &
                  0.209542, 0.197181, 0.181949, 0.163786, 0.143015, &
                  0.120465, 0.0974421, 0.0755035, 0.0560654, 0.0400306, &
                  0.0276386, 0.018578, 0.0122386, 0.00794693, 0.00510924, &
                  0.00326308, 0.00207494, 0.0013157 /)

integer i, j, n, ip, ix, sign_phi
logical, save :: init_needed = .true.

! Init

if (init_needed) then

  p%y_min = (/ 0.0, 0.90 /)
  p%y_max = (/ 0.9, 0.99 /)
  p%del_y = (/ 0.1, 0.01 /)

  dp%y_min = p%y_min
  dp%y_max = p%y_max
  dp%del_y = p%del_y

  allocate (spline%pt(0:1))

  do i = 1, size(p)

    allocate (p(i)%int_prob(0:9))
    allocate (dp(i)%int_prob(0:9))

    do j = 0, ubound(p(i)%int_prob, 1)

      allocate (p(i)%int_prob(j)%pt(0:5), dp(i)%int_prob(j)%pt(0:5))

      p(i)%int_prob(j)%spline_type  = gen_poly_spline$
      dp(i)%int_prob(j)%spline_type = gen_poly_spline$

      p(i)%int_prob(j)%x_min = -3   ! min log_E_rel of spline fit
      p(i)%int_prob(j)%x_max =  2   ! max log_E_rel of spline fit
      p(i)%int_prob(j)%del_x =  1   ! delta log_E_rel of spline fit

      dp(i)%int_prob(j)%x_min = -3
      dp(i)%int_prob(j)%x_max =  2
      dp(i)%int_prob(j)%del_x =  1

    enddo

  enddo

  p(1)%int_prob(0)%pt%c0 = 0
  p(1)%int_prob(1)%pt%c0 = (/ 0.213077, 0.209797, 0.197551, 0.172532, 0.159164, 0.157629 /)
  p(1)%int_prob(2)%pt%c0 = (/ 0.41684, 0.411585, 0.391252, 0.346399, 0.3208, 0.317791 /) 
  p(1)%int_prob(3)%pt%c0 = (/ 0.607859, 0.602204, 0.579339, 0.523204, 0.487668, 0.483317 /) 
  p(1)%int_prob(4)%pt%c0 = (/ 0.788241, 0.783377, 0.762758, 0.705161, 0.663198, 0.657732 /) 
  p(1)%int_prob(5)%pt%c0 = (/ 0.962935, 0.959643, 0.944795, 0.895588, 0.852122, 0.845915 /) 
  p(1)%int_prob(6)%pt%c0 = (/ 1.13857, 1.13739, 1.13097, 1.09987, 1.06173, 1.05541 /) 
  p(1)%int_prob(7)%pt%c0 = (/ 1.32443, 1.32584, 1.33026, 1.32768, 1.3048, 1.2995 /) 
  p(1)%int_prob(8)%pt%c0 = (/ 1.53718, 1.5418, 1.5601, 1.59954, 1.60838, 1.60643 /) 
  p(1)%int_prob(9)%pt%c0 = (/ 1.82222, 1.83133, 1.86944, 1.97602, 2.05286, 2.06088 /)

  p(1)%int_prob(0)%pt%c1 = 0
  p(1)%int_prob(1)%pt%c1 = (/ -0.00144017, -0.00616179, -0.0199551, -0.0239961, -0.00429677, -0.000326423 /)
  p(1)%int_prob(2)%pt%c1 = (/ -0.00228913, -0.00994776, -0.0340395, -0.0449884, -0.00840142, -0.000641762 /)
  p(1)%int_prob(3)%pt%c1 = (/ -0.002441, -0.010802, -0.039611, -0.0602353, -0.0120918, -0.000932808 /)
  p(1)%int_prob(4)%pt%c1 = (/ -0.00207945, -0.00937684, -0.0370852, -0.0673929, -0.0150743, -0.00118121 /)
  p(1)%int_prob(5)%pt%c1 = (/ -0.00139017, -0.00642462, -0.0280016, -0.0644311, -0.0169119, -0.00135801 /)
  p(1)%int_prob(6)%pt%c1 = (/ -0.000477501, -0.00239835, -0.013628, -0.0492536, -0.016872, -0.00141169 /)
  p(1)%int_prob(7)%pt%c1 = (/ 0.000632962, 0.00258251,   0.00566825, -0.0188059, -0.013547, -0.00123701 /)
  p(1)%int_prob(8)%pt%c1 = (/ 0.00200317,   0.00879492, 0.0310633, 0.0333115, -0.00364551, -0.000570791 /)
  p(1)%int_prob(9)%pt%c1 = (/ 0.00391001, 0.0175067, 0.0680914, 0.12689, 0.023931, 0.00155287 /)


  p(2)%int_prob(0)%pt%c0 = (/ 1.82222, 1.83133, 1.86944, 1.97602, 2.05286, 2.06088 /)
  p(2)%int_prob(1)%pt%c0 = (/ 1.85968, 1.86938, 1.91015, 2.02619, 2.11413, 2.12405 /)
  p(2)%int_prob(2)%pt%c0 = (/ 1.90013, 1.91048, 1.95413, 2.08052, 2.18098, 2.19314 /)
  p(2)%int_prob(3)%pt%c0 = (/ 1.94432, 1.95538, 2.00218, 2.14002, 2.25481, 2.26962 /)
  p(2)%int_prob(4)%pt%c0 = (/ 1.99333, 2.00516, 2.05547, 2.20616, 2.3376, 2.35564 /)
  p(2)%int_prob(5)%pt%c0 = (/ 2.04878, 2.0615, 2.11577, 2.28118, 2.43241, 2.45448 /)
  p(2)%int_prob(6)%pt%c0 = (/ 2.11333, 2.12709, 2.18599, 2.36873, 2.54422, 2.57151 /)
  p(2)%int_prob(7)%pt%c0 = (/ 2.19184, 2.20684, 2.27137, 2.47544, 2.68217, 2.7166 /)
  p(2)%int_prob(8)%pt%c0 = (/ 2.2948, 2.31143, 2.38333, 2.6157, 2.8661, 2.91133 /)
  p(2)%int_prob(9)%pt%c0 = (/ 2.45393, 2.47307, 2.55635, 2.83301, 3.15652, 3.22186 /)

  p(2)%int_prob(0)%pt%c1 = (/ 0.00391001, 0.0175067, 0.0680914, 0.12689, 0.023931, 0.00155287 /)
  p(2)%int_prob(1)%pt%c1 = (/ 0.00416291, 0.0186658, 0.0731002, 0.140798, 0.0290337, 0.00197257 /)
  p(2)%int_prob(2)%pt%c1 = (/ 0.00443623, 0.0199192, 0.0785332, 0.156168, 0.0349716, 0.00246967 /)
  p(2)%int_prob(3)%pt%c1 = (/ 0.00473494, 0.0212899, 0.0844925, 0.173348, 0.0419805, 0.00306785 /)
  p(2)%int_prob(4)%pt%c1 = (/ 0.00506619, 0.0228108, 0.0911251, 0.192841, 0.050408, 0.00380271 /)
  p(2)%int_prob(5)%pt%c1 = (/ 0.00544082, 0.0245318, 0.0986538, 0.21541, 0.0607979, 0.00473099 /)
  p(2)%int_prob(6)%pt%c1 = (/ 0.00587652, 0.0265346, 0.107443, 0.242307, 0.0740684, 0.00595081 /)
  p(2)%int_prob(7)%pt%c1 = (/ 0.00640543, 0.0289673, 0.118154, 0.275817, 0.0919565, 0.00765291 /)
  p(2)%int_prob(8)%pt%c1 = (/ 0.00709691, 0.0321499, 0.132216, 0.320901, 0.1184, 0.0102849 /)
  p(2)%int_prob(9)%pt%c1 = (/ 0.00816008, 0.0370469, 0.153944, 0.392638, 0.166096, 0.0153642 /)

  dp(1)%int_prob(0)%pt%c0 = (/ 2.1491, 2.11355, 1.98278, 1.72319, 1.58761, 1.57217 /)
  dp(1)%int_prob(1)%pt%c0 = (/ 2.09573, 2.06807, 1.96143, 1.72964, 1.59974, 1.58459 /)
  dp(1)%int_prob(2)%pt%c0 = (/ 1.97406, 1.96255, 1.90957, 1.75031, 1.63752, 1.6233 /)
  dp(1)%int_prob(3)%pt%c0 = (/ 1.85037, 1.85313, 1.85384, 1.78942, 1.70547, 1.69299 /)
  dp(1)%int_prob(4)%pt%c0 = (/ 1.7657, 1.77807, 1.82005, 1.8551, 1.81284, 1.80325 /)
  dp(1)%int_prob(5)%pt%c0 = (/ 1.73909, 1.75775, 1.82971, 1.9621, 1.97737, 1.97254 /)
  dp(1)%int_prob(6)%pt%c0 = (/ 1.7884, 1.81187, 1.90804, 2.13883, 2.23489, 2.23823 /)
  dp(1)%int_prob(7)%pt%c0 = (/ 1.95429, 1.98284, 2.10435, 2.44847, 2.66746, 2.68623 /)
  dp(1)%int_prob(8)%pt%c0 = (/ 2.36249, 2.39907, 2.55887, 3.07053, 3.51229, 3.56615 /)
  dp(1)%int_prob(9)%pt%c0 = (/ 3.61504, 3.67255, 3.92872, 4.83658, 5.88372, 6.05988 /)

  dp(1)%int_prob(0)%pt%c1 = (/ -0.0156567, -0.0665479, -0.210823, -0.245065, -0.0432786, -0.00328204 /)
  dp(1)%int_prob(1)%pt%c1 = (/ -0.0120549, -0.0523352, -0.177843, -0.229804, -0.0423414, -0.00322828 /)
  dp(1)%int_prob(2)%pt%c1 = (/ -0.00480825, -0.0226506, -0.0995454, -0.185447, -0.0393877, -0.00305697 /)
  dp(1)%int_prob(3)%pt%c1 = (/ 0.00141279, 0.00427155, -0.0129626, -0.115597, -0.0339413, -0.00273419 /)
  dp(1)%int_prob(4)%pt%c1 = (/ 0.00549944, 0.0229272, 0.0606838, -0.024231, -0.0250039, -0.00218733 /)
  dp(1)%int_prob(5)%pt%c1 = (/ 0.00811466, 0.0353651, 0.118791, 0.0868289, -0.0105808, -0.00126711 /)
  dp(1)%int_prob(6)%pt%c1 = (/ 0.010096, 0.0449573, 0.167992, 0.22149, 0.0135693, 0.000355742 /)
  dp(1)%int_prob(7)%pt%c1 = (/ 0.0122038, 0.0550616, 0.219601, 0.397048, 0.0578117, 0.0035227 /)
  dp(1)%int_prob(8)%pt%c1 = (/ 0.0155698, 0.0708746, 0.295621, 0.671762, 0.154618, 0.0110212 /)
  dp(1)%int_prob(9)%pt%c1 = (/ 0.0243979, 0.11179, 0.482324, 1.32675, 0.474496, 0.0386947 /)

  dp(2)%int_prob(0)%pt%c0 = (/ 3.61504, 3.67255, 3.92872, 4.83658, 5.88372, 6.05988 /)
  dp(2)%int_prob(1)%pt%c0 = (/ 3.88478, 3.94666, 4.22274, 5.21077, 6.3868, 6.59219 /)
  dp(2)%int_prob(2)%pt%c0 = (/ 4.21748, 4.2847, 4.58521, 5.67121, 7.00664, 7.24929 /)
  dp(2)%int_prob(3)%pt%c0 = (/ 4.6385, 4.71247, 5.04371, 6.25269, 7.79078, 8.08234 /)
  dp(2)%int_prob(4)%pt%c0 = (/ 5.18941, 5.27217, 5.64343, 7.01222, 8.81726, 9.17548 /)
  dp(2)%int_prob(5)%pt%c0 = (/ 5.94354, 6.03829, 6.46406, 8.05031, 10.2241, 10.6779 /)
  dp(2)%int_prob(6)%pt%c0 = (/ 7.04412, 7.1563, 7.6613, 9.56334, 12.2819, 12.883 /)
  dp(2)%int_prob(7)%pt%c0 = (/ 8.81611, 8.95623, 9.58825, 11.9967, 15.6067, 16.4608 /)
  dp(2)%int_prob(8)%pt%c0 = (/ 12.1997, 12.3929, 13.2664, 16.6389, 21.9895, 23.3684 /)
  dp(2)%int_prob(9)%pt%c0 = (/ 21.6485, 21.9891, 23.5329, 29.5915, 39.9738, 43.0058 /)

  dp(2)%int_prob(0)%pt%c1 = (/ 0.0243979, 0.11179, 0.482324, 1.32675, 0.474496, 0.0386947 /)
  dp(2)%int_prob(1)%pt%c1 = (/ 0.0262411, 0.120306, 0.520659, 1.45904, 0.548779, 0.0455162 /)
  dp(2)%int_prob(2)%pt%c1 = (/ 0.0285039, 0.130756, 0.567609, 1.62067, 0.642666, 0.0542909 /)
  dp(2)%int_prob(3)%pt%c1 = (/ 0.0313552, 0.143918, 0.626653, 1.82358, 0.764771, 0.0659247 /)
  dp(2)%int_prob(4)%pt%c1 = (/ 0.0350713, 0.161068, 0.703478, 2.08733, 0.929534, 0.0819626 /)
  dp(2)%int_prob(5)%pt%c1 = (/ 0.0401387, 0.184448, 0.808106, 2.44651, 1.16311, 0.105256 /)
  dp(2)%int_prob(6)%pt%c1 = (/ 0.0475064, 0.218436, 0.960088, 2.96879, 1.51815, 0.141676 /)
  dp(2)%int_prob(7)%pt%c1 = (/ 0.059322, 0.272936, 1.2037, 3.80805, 2.11831, 0.205381 /)
  dp(2)%int_prob(8)%pt%c1 = (/ 0.0817827, 0.376533, 1.66676, 5.4108, 3.33731, 0.340647 /)
  dp(2)%int_prob(9)%pt%c1 = (/ 0.144122, 0.664085, 2.95279, 9.90065, 7.05788, 0.782273 /)

  do i = 1, 2
    do j = 0, ubound(p(i)%int_prob, 1)
      call photon_init_spline_coef_calc (p(i)%int_prob(j))
      call photon_init_spline_coef_calc (dp(i)%int_prob(j))
    enddo
  enddo

  init_needed = .false.
endif

!----------------------------------------------------
! Integrated probability number to use.

if (present(r_in)) then
  rr = r_in
else
  call ran_uniform(rr)
endif

! The spline fit is only for positive phi.
! So make phi negative if rr < 0.5

if (rr > 0.5) then
  rr = 2 * rr - 1
  sign_phi = 1
else
  rr = 1 - 2 * rr
  sign_phi = -1
endif

! gamma_phi is a function of E_rel and rr.
! gamma_phi is calculated by first finding the spline coefficients for the given E_rel

log_E_rel = log10(E_rel)

! In the range above rr = 0.99 we use an extrapolation that matches gamma_phi and
! it's first two derivatives at rr = 0.99.

if (rr >= p(2)%y_max) then

  spline%x_min = p(2)%y_max
  spline%x_max = 1.00000001
  spline%del_x = 1 - p(2)%y_max
  spline%spline_type = end_spline$

  n = ubound(p(2)%int_prob, 1)
  x1  = photon_init_spline_fit (p(2)%int_prob(n), log_E_rel)
  xp1 = photon_init_spline_fit (dp(2)%int_prob(n), log_E_rel) * spline%del_x

  x0  = photon_init_spline_fit (p(2)%int_prob(n-1), log_E_rel)
  xp0 = photon_init_spline_fit (dp(2)%int_prob(n-1), log_E_rel) * spline%del_x

  v = x1 - x0 - xp0
  vp = xp1 - xp0

  spline%pt(0)%c0 = x1
  spline%pt(0)%c1 = xp1 
  spline%pt(0)%c2 = max(vp/2, 2*vp - 3*v) 
  spline%pt(0)%c3 = 1.1

  gamma_phi = photon_init_spline_fit (spline, rr)

else
  
  do ip = 1, 2
    if (rr > p(ip)%y_max) cycle
    r = (rr - p(ip)%y_min) / p(ip)%del_y
    i = int (r)

    spline%pt(0)%c0 = photon_init_spline_fit (p(ip)%int_prob(i),  log_E_rel)
    spline%pt(0)%c1 = photon_init_spline_fit (dp(ip)%int_prob(i), log_E_rel)

    if (i == ubound(p(ip)%int_prob, 1)) then  ! Can happen due to roundoff errors.
      gamma_phi = spline%pt(0)%c0
      exit            ! This prevents an array out-of-bounds problem
    endif

    spline%pt(1)%c0 = photon_init_spline_fit (p(ip)%int_prob(i+1),  log_E_rel)
    spline%pt(1)%c1 = photon_init_spline_fit (dp(ip)%int_prob(i+1), log_E_rel)

    spline%x_min = i * p(ip)%del_y + p(ip)%y_min
    spline%x_max = spline%x_min + p(ip)%del_y 
    spline%del_x = p(ip)%del_y
    spline%spline_type = gen_poly_spline$

    call photon_init_spline_coef_calc (spline)
    gamma_phi = photon_init_spline_fit (spline, rr)
    exit

  enddo

endif

! Scale result by the sigma of the spectrum at fixed E_rel

if (E_rel < 0.1) then
  sig = 0.597803 * E_rel**(-0.336351)
else
  sig = 0.451268 * E_rel**(-0.469377)
endif

gamma_phi = gamma_phi * sig * sign_phi



end subroutine photon_vert_angle_init

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!+
! Subroutine photon_energy_init (E_rel, r_in)
!
! Routine to convert a random number in the interval [0,1] to a photon energy.
! The photon probability splinetrum is:
!   P(E_rel) = 0.1909859 * Integral_{E_rel}^{Infty} K_{5/3}(x) dx
! Where
!   P(E_rel)) = Probability of finding a photon at relative energy E_rel.
!   E_rel     = Relative photon energy: E / E_crit, E_crit = Critical energy.
!   K_{5/3}   = Modified Bessel function.
! There is a cut-off built into the calculation so that E_rel will be in the 
! range [0, ~17]. The error in neglecting photons with E_rel > ~17 is very small. 
! If r_in is given: 
!   r_in = 0 => E_rel = 0 
!   r_in = 1 => E_rel = ~17
!
! Module needed:
!   use photon_init_mod
!
! Input:
!   r_in  -- Real(rp), optional: Integrated probability in the range [0,1].
!             If not present, a random number will be used.
!
! Output:
!   E_rel -- Real(rp): Relative photon energy E/E_crit. 
!-

subroutine photon_energy_init (E_rel, r_in)

implicit none

! Four spline fit arrays are used. 
! Each fit array has a different del_x and range of validity.

type (photon_init_spline_struct), save :: spline(5)  

real(rp) E_rel
real(rp), optional :: r_in
real(rp) rr, x, r_rel, rr0, x0, xp0, x1, xp1, v, vp
real(rp), save :: a, b, c, n

integer i, is, ns

logical, save :: init_needed = .true.

! Check for r_in

if (present(r_in)) then
  rr = r_in
  if (rr < 0  .or. rr > 1) then
    print *, 'ERROR: RR IS OUT OF RANGE: ', rr
    stop
  endif
else
  call ran_uniform(rr)
endif

! Init. 
! The values for c0 and c1 were obtained from a Mathematica calculation.

if (init_needed) then

  spline(:)%del_x = (/ 0.02, 0.01, 0.001, 0.0001, 0.0001 /)
  spline(:)%x_min = (/ 0.0,  0.80, 0.990, 0.9990, 0.9999 /)
  spline(:)%x_max = (/ 0.8,  0.99, 0.999, 0.9999, 1.00000001 /)
  spline(:)%spline_type = gen_poly_spline$

  allocate (spline(1)%pt(0:40), spline(2)%pt(0:19), spline(3)%pt(0:9), spline(4)%pt(0:9))
  allocate (spline(5)%pt(0:1))

  spline(1)%pt(:)%c0 = (/ 0.0, &
              4.28341e-6, 0.0000342902, 0.000115858, &
              0.000275057, 0.000538307, 0.0009325, 0.00148513, 0.00222445, &
              0.00317957, 0.00438067, 0.00585915, 0.00764782, 0.00978112, &
              0.0122954, 0.015229, 0.0186228, 0.0225205, 0.026969, 0.0320186, &
              0.037724, 0.0441446, 0.0513452, 0.0593972, 0.0683791, 0.0783781, &
              0.0894916, 0.101829, 0.115512, 0.130682, 0.147499, 0.166146, &
              0.186838, 0.209825, 0.235404, 0.263928, 0.295827, 0.331625, 0.371975, &
              0.417703, 0.469878 /)

  spline(1)%pt(:)%c1 = (/ 0.0, & 
              0.000642606, 0.00257329, 0.00580068, 0.0103393, 0.0162097, 0.0234387, &
              0.0320599, 0.042114, 0.0536493, 0.0667223, 0.081399, 0.097755, &
              0.115877, 0.135866, 0.157835, 0.181913, 0.20825, 0.237015, 0.268401, &
              0.302633, 0.339965, 0.380692, 0.425158, 0.473758, 0.526956, 0.585298, &
              0.649425, 0.720102, 0.798242, 0.884947, 0.981556, 1.08971, 1.21146, &
              1.34935, 1.50665, 1.68757, 1.89762, 2.1442, 2.43746, 2.79164 /)

  spline(2)%pt(:)%c0 = (/ &
              0.469878, 0.498807, 0.529911, 0.563449, 0.599724, 0.639103, 0.682028, &
              0.72904, 0.780816, 0.83821, 0.902325, 0.974622, 1.05709, 1.15252, 1.26503, &
              1.40105, 1.57139, 1.79651, 2.12271, 2.69945 /)

  spline(2)%pt(:)%c1 = (/ & 
              2.79164, 2.9977, 3.22744, 3.48512, 3.77607, 4.10707, 4.48684, 4.92685, &
              5.44242, 6.05453, 6.79263, 7.69945, 8.83931, 10.3137, 12.2927, 15.0837, &
              19.3053, 26.4104, 40.7923, 84.6624 /)

  spline(3)%pt(:)%c0 = (/ &
              2.69945, 2.78886, 2.88928, 3.0037, 3.13651, 3.29451, 3.48916, 3.742, &
              4.10159, 4.72378 /)

  spline(3)%pt(:)%c1 = (/ & 
              84.6624, 94.5006, 106.83, 122.726, 143.986, 173.851, 218.823, 294.119, &
              445.575, 903.65 /)

  spline(4)%pt(:)%c0 = (/ &
              4.72378, 4.81908, 4.92582, 5.04708, 5.18739, 5.35376, 5.55797, 5.82214, &
              6.19606, 6.83908 /)

  spline(4)%pt(:)%c1 = (/ & 
              903.65, 1005.91, 1133.9, 1298.7, 1518.78, 1827.44, 2291.39, 3066.54, 4621.7, &
              9308.58 /)

  ! In the range above rr = 0.9999 we use a simple extrapolation.

  n = ubound(spline(4)%pt, 1)
  x1 = spline(4)%pt(n)%c0
  xp1 = spline(4)%pt(n)%c1 * spline(5)%del_x

  x0 = spline(4)%pt(n-1)%c0
  xp0 = spline(4)%pt(n-1)%c1 * spline(5)%del_x

  v = x1 - x0 - xp0
  vp = xp1 - xp0

  spline(5)%pt(0)%c0 = x1
  spline(5)%pt(0)%c1 = xp1
  spline(5)%pt(0)%c2 = max(vp/2, 2*vp - 3*v) 
  spline(5)%pt(0)%c3 = 1.1


  ! Fill in rest of the spline fit coefs.

  do is = 1, ubound(spline, 1) - 1
    call photon_init_spline_coef_calc (spline(is))
  enddo

  init_needed = .false.
endif

! Spline fit

do is = 1, ubound(spline, 1)
  if (rr > spline(is)%x_max) cycle
  E_rel = photon_init_spline_fit (spline(is), rr)
  return
enddo

end subroutine photon_energy_init

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!+
! Subroutine photon_init_spline_coef_calc (spline)
!
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
! Function photon_init_spline_fit (spline, rr) result (fit_val)
!
! Function to evaluate a spline fit for some value of the integrated probability rr.
! 
! Module needed:
!   use photon_init_mod
!
! Input:
!   spline    -- photon_init_spline_struct: spline section.
!   r_in      -- Real(rp), Integrated probability in the range [0,1].
!                 If not present, a random number will be used.
!
! Output:
!   fit_val -- Real(rp): Spline fit evaluated at r_in.
!               Note: if r_in is out of range, fit_val will be 
!               set to the value at the edge of the spline range.
!-

function photon_init_spline_fit (spline, rr) result (fit_val)

implicit none

type (photon_init_spline_struct) spline
real(rp) rr, r_rel, x, fit_val, v, vp
integer i, ix, np

! Find in which spline section rr is in

r_rel = (rr - spline%x_min) / spline%del_x
i = int(r_rel)   ! Index of which spline section to use.
x = r_rel - i    ! Notice that x will be in the range [0, 1].

! If out of range then adjust point to be at the edge of the range

if (i < 0) then
  i = 0
  x = 0
endif

np = ubound(spline%pt, 1)
if (i > np - 1) then 
  i = np - 1
  x = 1
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
    fit_val = spline%pt(i)%c0 + spline%pt(i)%c1 * x + &
                                    spline%pt(i)%c2 * x**2 + spline%pt(i)%c3 * x**3
  endif

case (end_spline$)

  fit_val = spline%pt(i)%c0 + spline%pt(i)%c1 * x + &
                                    spline%pt(i)%c2 * x**2 / (1 - x/spline%pt(i)%c3)

! Coding error if here.

case default
  call err_exit

end select

end function photon_init_spline_fit

end module
