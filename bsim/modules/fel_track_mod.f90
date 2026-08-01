!+
! Module fel_track_mod
!
! The steady-state FEL step: period-averaged particle push, source deposition and FFT
! field solve inside an undulator segment, operating on the Bmad-coordinate packed beam of
! fel_beam_mod. The physics is transcribed from Genesis 1.3 Version 4, which GPL permits;
! sources, by routine:
!
!   fel_transverse_track   <- TrackBeam::track, applyDrift, applyFQuad, applyDQuad
!                             (src/Core/TrackBeam.cpp:8-116)
!   fel_advance            <- BeamSolver::advance (src/Core/BeamSolver.cpp:13-87)
!   fel_runge_kutta        <- BeamSolver::RungeKutta (BeamSolver.cpp:89-141), verbatim
!   fel_ode                <- BeamSolver::ODE (BeamSolver.cpp:144-163)
!   fel_field_step         <- FieldSolverFFT::advance, FFT, init, unfiltered path
!   fel_und_coupling       <- Undulator::fc (src/Core/Undulator.cpp:141-167)
!   faw, faw2              <- Undulator::faw, faw2 (Undulator.cpp:170-182)
!
! and the step composition (transverse half step, longitudinal advance, transverse half
! step, field step) is Beam::track plus Gencore's step 4.
!
! Coordinates. The stored state is Bmad's (x, px/p0, y, py/p0, z, pz) plus weight; see
! fel_beam_mod. The longitudinal RK4, however, runs in Genesis's chart (theta, gamma) as
! per-particle working variables, derived at step entry and converted back at step exit:
!
!   entry:  gamma = sqrt((p0_mc*(1+pz))^2 + 1),  theta = phi0 + ks*z/beta
!   exit:   pz = (sqrt(gamma^2-1) - p0_mc)/p0_mc,  z = -beta*(phi0_new - theta)/ks
!
! Two facts make this clean rather than a compromise. First, theta <-> (phi0, z) is an
! affine change of variables with a state-independent shift, and RK4 is exactly invariant
! under affine maps, so integrating in theta is not a different discretization from
! integrating in z -- it lets the audited, verbatim RK4 transcription survive unchanged.
! Second, the per-step gamma <-> pz round trip costs ~1 ulp of gamma per step (the
! subtraction in pz is exact once P is formed), which is the same order as the step's own
! rounding. The common reference phase phi0 advances once per step per beam
! (fel_phi0_rate in fel_beam_mod); the split keeps z small and bunch-scale.
!
! Weights: the source deposition scales per particle as c*w_j/slicelength where Genesis
! has current/N -- identical algebra for uniform weights, and correct for nonuniform ones
! (each macroparticle radiates its own charge). All in Genesis internal field units; see
! fel_field_unit_scale.
!
! Deliberately absent, per the deliverable: time dependence, slippage, space charge,
! wakes, incoherent synchrotron radiation, harmonics beyond the coupling formula, the
! source filter, and orbit or field errors.
!-

module fel_track_mod

use fel_beam_mod
use wavefront_mod

implicit none

!+
! Struct fel_und_struct
!
! One undulator segment, constant parameters along it: one contiguous run of aw > 0 steps
! in Genesis's unrolled lattice, divided into nstep equal steps of dz = l/nstep with
! nstep = round(l/delz_target). kx, ky carry Genesis's unroll scaling by ku^2
! (Lattice.cpp:412-413).
!-

type fel_und_struct
  real(rp) :: aw = 0          ! rms undulator parameter.
  real(rp) :: ku = 0          ! Undulator wavenumber twopi/lambdau [1/m].
  real(rp) :: kx = 0          ! Natural focusing, deck kx * ku^2 [1/m^2].
  real(rp) :: ky = 0          ! Natural focusing, deck ky * ku^2 [1/m^2].
  real(rp) :: ax = 0, ay = 0  ! Transverse offset of the undulator field [m].
  logical :: helical = .false.
  integer :: nstep = 0        ! Number of integration steps over the segment.
  real(rp) :: dz = 0          ! Step length [m].
end type

! Cached kernel state for fel_field_step, mirroring FieldSolverFFT::init. Module state,
! single threaded, same caveat as the wavefront plan cache.

complex(rp), allocatable, private, save :: fel_k2(:,:)        ! -i (kx^2+ky^2)/(2 ks), FFT order.
complex(rp), allocatable, private, save :: fel_crsource(:,:)  ! Source accumulator.
integer, private, save :: fel_cache_ngrid = 0
real(rp), private, save :: fel_cache_dgrid = 0, fel_cache_ks = 0

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function fel_field_unit_scale (wf) result (scale)
!
! Routine to return the factor converting a field in V/m (the wavefront_struct
! convention) to Genesis internal units. Multiply on entry to FEL tracking, divide on
! exit. Composition of the dump relation dfl = u*dgrid*eev/(ks*sqrt(vacimp))
! (writeFieldHDF5.cpp:70) with the reader's E = dfl*sqrt(2 Z0)/dgrid; the dgrid factors
! cancel exactly. Z0 is Bmad's mu_0_vac*c_light because that is what the wavefront reader
! used; vacimp and eev are Genesis's because they define the internal units.
!-

function fel_field_unit_scale (wf) result (scale)

type (wavefront_struct) wf
real(rp) scale, ks

!

ks = twopi / wf%wavelength
scale = ks * sqrt(fel_vacimp) / (sqrt(2 * (mu_0_vac * c_light)) * fel_eev)

end function fel_field_unit_scale

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function fel_und_coupling (und, h) result (fc)
!
! The coupling factor fc(h), transcribed from Undulator::fc: helical couples the
! fundamental with strength aw and nothing else; planar couples odd harmonics through
! JJ = J_{(h-1)/2}(xi) - J_{(h+1)/2}(xi), xi = h/2 * aw^2/(1+aw^2), sign (-1)^((h-1)/2);
! even harmonics not at all (a polarisation limitation, brief section 5).
!-

function fel_und_coupling (und, h) result (fc)

type (fel_und_struct) und
integer h, h0, h1
real(rp) fc, xi, coup

!

coup = und%aw

if (und%helical) then
  if (h == 1) then
    fc = coup
  else
    fc = 0
  endif
  return
endif

if (mod(h, 2) == 1) then
  xi = und%aw**2
  xi = 0.5_rp * xi / (1 + xi) * h
  h0 = (h - 1) / 2
  h1 = h0 + 1
  fc = coup * (bessel_jn(h0, xi) - bessel_jn(h1, xi)) * (-1)**h0
else
  fc = 0
endif

end function fel_und_coupling

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function faw (und, x, y) result (value)
!
! Transverse dependence of the undulator field, first order form used in the particle
! gather. Undulator::faw.
!-

function faw (und, x, y) result (value)

type (fel_und_struct) und
real(rp) x, y, value, dx, dy

!

dx = x - und%ax
dy = y - und%ay
value = 1 + 0.5_rp * (und%kx * dx*dx + und%ky * dy*dy)

end function faw

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function faw2 (und, x, y) result (value)
!
! Square of the transverse dependence, used in the source deposition. Undulator::faw2.
!-

function faw2 (und, x, y) result (value)

type (fel_und_struct) und
real(rp) x, y, value, dx, dy

!

dx = x - und%ax
dy = y - und%ay
value = 1 + und%kx * dx*dx + und%ky * dy*dy

end function faw2

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_track_und_step (und, beam, wf, err_flag)
!
! Routine to advance the whole beam and field one integration step of length und%dz
! inside an undulator segment, in Genesis's order: transverse half step, longitudinal RK4
! (gathering the field of this moment and holding it through the stages), transverse half
! step, then source deposition and field solve. The common phase phi0 advances once per
! step, with the undulator's ku as the reference wavenumber.
!
! Input:
!   und         -- fel_und_struct: Segment parameters.
!   beam        -- fel_beam_struct: The beam. Steady state: one slice, periodic.
!   wf          -- wavefront_struct: The field, in Genesis internal units.
!
! Output:
!   beam, wf    -- Advanced one step.
!   err_flag    -- logical: Set True on error.
!-

subroutine fel_track_und_step (und, beam, wf, err_flag)

type (fel_und_struct) und
type (fel_beam_struct), target :: beam
type (wavefront_struct) wf
logical err_flag

real(rp) ks, phi0_new
integer is

!

err_flag = .true.
ks = twopi / wf%wavelength
phi0_new = beam%phi0 + und%dz * fel_phi0_rate(ks, und%ku, beam%p0_mc)

do is = 1, size(beam%slice)
  call fel_transverse_track (und, beam, beam%slice(is), und%dz/2)
  call fel_advance (und, beam, beam%slice(is), wf, und%dz, phi0_new)
  call fel_transverse_track (und, beam, beam%slice(is), und%dz/2)
enddo

beam%phi0 = phi0_new

do is = 1, size(beam%slice)
  call fel_field_step (und, beam, beam%slice(is), wf, und%dz, err_flag)
  if (err_flag) return
enddo

end subroutine fel_track_und_step

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_track_interlude_genesis (qf, length, beam, wf, err_flag)
!
! Routine to advance one field-free interlude element -- a drift or a quadrupole -- the
! way Genesis does it, as one integration step: transverse half step, the collapsed
! theta advance with the path-length term sampled at mid element and Genesis's drift
! reference xku = ks*0.5/gamma0/gamma0 (BeamSolver.cpp:35-38, division order kept),
! transverse half step, field diffraction with zero source.
!
! This is the Genesis interlude model, transcribed, so that running the full lattice with
! it isolates what the Bmad seam changes: the seam integrates the path-length term
! exactly through the quad map where this samples it once. See the benchmark README; the
! production configuration is the seam.
!-

subroutine fel_track_interlude_genesis (qf, length, beam, wf, err_flag)

type (fel_beam_struct), target :: beam
type (fel_slice_struct), pointer :: sl
type (wavefront_struct), target :: wf
real(rp) qf, length
logical err_flag

type (fel_und_struct) und0
real(rp) xks, xku, qquad, phi0_new
real(rp) px_g, py_g, gam, beta, theta, btpar, btpar0, slope, gz_hat, q_hat
integer is, ip
logical err

!

err_flag = .true.

xks = twopi / wf%wavelength
xku = xks * 0.5_rp / beam%gamma0 / beam%gamma0    ! Genesis's division order, kept.
qquad = qf * beam%gamma0                          ! TrackBeam.cpp:26.
q_hat = qquad / beam%p0_mc
phi0_new = beam%phi0 + length * fel_phi0_rate(xks, xku, beam%p0_mc)

do is = 1, size(beam%slice)
  sl => beam%slice(is)

  call interlude_transverse_half (sl, q_hat, length/2)

  ! theta advance, the collapsed RK4: with no field and no space charge the slope is
  ! theta independent and constant through the stages, so RK4 reduces to one exact step.
  ! btpar = 1 + px_g^2 + py_g^2 (aw = 0), px_g = gamma*beta_x = px * p0_mc.

  do ip = 1, sl%n
    gam = fel_gamma_of(beam%p0_mc, sl%pz(ip))
    beta = fel_beta_of(beam%p0_mc, sl%pz(ip))
    theta = beam%phi0 + xks * sl%z(ip) / beta

    px_g = sl%px(ip) * beam%p0_mc
    py_g = sl%py(ip) * beam%p0_mc
    btpar = 1 + px_g*px_g + py_g*py_g
    btpar0 = sqrt(1 - btpar / (gam * gam))
    slope = xks * (1 - 1/btpar0) + xku
    theta = theta + length * slope

    ! Back to z: tau = (phi0_new - theta)/ks; pz unchanged (no energy change), so beta
    ! is unchanged too.
    sl%z(ip) = -beta * (phi0_new - theta) / xks
  enddo

  call interlude_transverse_half (sl, q_hat, length/2)
enddo

beam%phi0 = phi0_new

! Field diffraction: Genesis's kernel, zero source (aw = 0 makes the coupling exactly
! zero, matching Genesis skipping the deposition when not in an undulator).

und0%aw = 0
do is = 1, size(beam%slice)
  call fel_field_step (und0, beam, beam%slice(is), wf, length, err);  if (err) return
enddo

err_flag = .false.

!------------------------------------------------------------------------------
contains

subroutine interlude_transverse_half (sl, q_hat, dz)

type (fel_slice_struct) sl
real(rp) q_hat, dz, gz_hat
integer ip

do ip = 1, sl%n
  gz_hat = sqrt((1 + sl%pz(ip))**2 - sl%px(ip)**2 - sl%py(ip)**2)
  call fel_apply_focus (dz,  q_hat, sl%x(ip), sl%px(ip), gz_hat)
  call fel_apply_focus (dz, -q_hat, sl%y(ip), sl%py(ip), gz_hat)
enddo

end subroutine interlude_transverse_half

end subroutine fel_track_interlude_genesis

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_transverse_track (und, beam, sl, delz)
!
! Routine to advance the transverse coordinates over delz inside an undulator: drift plus
! natural focusing, transcribed from TrackBeam::track (quad, corrector and chicane
! branches dropped; they never overlap undulator fields here).
!
! In the stored normalization the per particle longitudinal factor is
!
!   gz_hat = gammaz/p0 = sqrt((1+pz)^2 - px^2 - py^2 - (aw/p0_mc)^2)
!
! (Genesis's gammaz = sqrt(gamma^2 - 1 - aw^2 - px_g^2 - py_g^2), divided by p0), and the
! maps read identically with px in Bmad units:
!
!   drift:  x += px*delz/gz_hat
!   quad:   foc^2 = q_hat/gz_hat with q_hat = q/p0_mc;  x' = a1 x + a2 px/gz_hat;
!           px' = a3 x gz_hat + a1 px
!
! The effective strengths (TrackBeam.cpp:23-31): qnat_{x,y} = k_{x,y}*aw^2/(gamma0*
! betpar0), betpar0 = sqrt(1 - (1+aw^2)/gamma0^2), with Genesis's reference gamma.
!-

subroutine fel_transverse_track (und, beam, sl, delz)

type (fel_und_struct) und
type (fel_beam_struct) beam
type (fel_slice_struct) sl
real(rp) delz
real(rp) betpar0, qx_hat, qy_hat, aw2, aw_p0_sq, gz_hat
integer ip

!

aw2 = und%aw**2
betpar0 = sqrt(1 - (1 + aw2)/beam%gamma0**2)

qx_hat = und%kx * aw2 / beam%gamma0 / betpar0 / beam%p0_mc
qy_hat = und%ky * aw2 / beam%gamma0 / betpar0 / beam%p0_mc
aw_p0_sq = (und%aw / beam%p0_mc)**2

do ip = 1, sl%n
  gz_hat = sqrt((1 + sl%pz(ip))**2 - sl%px(ip)**2 - sl%py(ip)**2 - aw_p0_sq)
  call fel_apply_focus (delz, qx_hat, sl%x(ip), sl%px(ip), gz_hat)
  call fel_apply_focus (delz, qy_hat, sl%y(ip), sl%py(ip), gz_hat)
enddo

end subroutine fel_transverse_track

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_apply_focus (delz, q_hat, x, px, gz_hat)
!
! One transverse plane's map over delz: TrackBeam::applyDrift for q_hat = 0, applyFQuad
! for q_hat > 0, applyDQuad for q_hat < 0, in the normalized variables of
! fel_transverse_track (px = P_x/p0, gz_hat = gammaz/p0, q_hat = q/p0).
!-

subroutine fel_apply_focus (delz, q_hat, x, px, gz_hat)

real(rp) delz, q_hat, x, px, gz_hat
real(rp) foc, omg, a1, a2, a3, xtmp

!

if (q_hat == 0) then
  x = x + px * delz / gz_hat
elseif (q_hat > 0) then
  foc = sqrt(q_hat/gz_hat)
  omg = foc * delz
  a1 = cos(omg)
  a2 = sin(omg)/foc
  a3 = -a2 * foc * foc
  xtmp = x
  x  = a1 * xtmp + a2 * px / gz_hat
  px = a3 * xtmp * gz_hat + a1 * px
else
  foc = sqrt(-q_hat/gz_hat)
  omg = foc * delz
  a1 = cosh(omg)
  a2 = sinh(omg)/foc
  a3 = a2 * foc * foc
  xtmp = x
  x  = a1 * xtmp + a2 * px / gz_hat
  px = a3 * xtmp * gz_hat + a1 * px
endif

end subroutine fel_apply_focus

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_advance (und, beam, sl, wf, delz, phi0_new)
!
! Routine to advance the longitudinal plane of every particle over delz. Transcribed from
! BeamSolver::advance: gather the field at (x, y) by bilinear interpolation, form
! rpart = (fc(1)/ks)*faw*conj(E), and integrate (theta, gamma) by the verbatim RK4 with
! rpart, px, py and faw held fixed through the stages.
!
! (theta, gamma) are derived at entry from the stored (z, pz) and written back at exit
! using phi0_new, the common phase at the end of this step; see the module header for why
! this chart change is exact for RK4 and ~1 ulp for the energy.
!-

subroutine fel_advance (und, beam, sl, wf, delz, phi0_new)

type (fel_und_struct) und
type (fel_beam_struct) beam
type (fel_slice_struct) sl
type (wavefront_struct) wf
real(rp) delz, phi0_new

real(rp) xks, xku, aw, rtmp, awloc, btpar, gamma, theta, beta, wx, wy, px_g, py_g, p_mc
complex(rp) cpart, rpart
integer ip, ix, iy
logical on_grid

!

xks = twopi / wf%wavelength
xku = und%ku
aw = und%aw
rtmp = fel_und_coupling(und, 1) / xks     ! fc(harm)/field->xks, harm = 1.

do ip = 1, sl%n
  gamma = fel_gamma_of(beam%p0_mc, sl%pz(ip))
  beta = fel_beta_of(beam%p0_mc, sl%pz(ip))
  theta = beam%phi0 + xks * sl%z(ip) / beta

  awloc = faw(und, sl%x(ip), sl%y(ip))
  px_g = sl%px(ip) * beam%p0_mc
  py_g = sl%py(ip) * beam%p0_mc
  btpar = 1 + px_g*px_g + py_g*py_g + aw*aw*awloc*awloc

  call fel_grid_weights (wf, sl%x(ip), sl%y(ip), ix, iy, wx, wy, on_grid)
  if (on_grid) then
    cpart =         wf%Ex(ix,   iy,   1) * wx * wy
    cpart = cpart + wf%Ex(ix+1, iy,   1) * (1-wx) * wy
    cpart = cpart + wf%Ex(ix,   iy+1, 1) * wx * (1-wy)
    cpart = cpart + wf%Ex(ix+1, iy+1, 1) * (1-wx) * (1-wy)
    rpart = rtmp * awloc * conjg(cpart)
  else
    rpart = 0
  endif

  call fel_runge_kutta (delz, xks, xku, btpar, rpart, gamma, theta)

  ! Back to the stored chart: pz from gamma (the subtraction is exact once p_mc is
  ! formed), z from tau = (phi0_new - theta)/ks with the updated beta.

  p_mc = sqrt(gamma**2 - 1)
  sl%pz(ip) = (p_mc - beam%p0_mc) / beam%p0_mc
  beta = p_mc / gamma
  sl%z(ip) = -beta * (phi0_new - theta) / xks
enddo

end subroutine fel_advance

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_runge_kutta (delz, xks, xku, btpar, rpart, gamma, theta)
!
! The RK4 stage bookkeeping of BeamSolver::RungeKutta (BeamSolver.cpp:89-141), verbatim.
!-

subroutine fel_runge_kutta (delz, xks, xku, btpar, rpart, gamma, theta)

real(rp) delz, xks, xku, btpar, gamma, theta
complex(rp) rpart
real(rp) k2gg, k2pp, k3gg, k3pp, stpz

! first step

k2gg = 0
k2pp = 0

call fel_ode (gamma, theta, xks, xku, btpar, rpart, k2gg, k2pp)

! second step

stpz = 0.5_rp * delz

gamma = gamma + stpz * k2gg
theta = theta + stpz * k2pp

k3gg = k2gg
k3pp = k2pp

k2gg = 0
k2pp = 0

call fel_ode (gamma, theta, xks, xku, btpar, rpart, k2gg, k2pp)

! third step

gamma = gamma + stpz * (k2gg - k3gg)
theta = theta + stpz * (k2pp - k3pp)

k3gg = k3gg / 6
k3pp = k3pp / 6

k2gg = k2gg * (-0.5_rp)
k2pp = k2pp * (-0.5_rp)

call fel_ode (gamma, theta, xks, xku, btpar, rpart, k2gg, k2pp)

! fourth step

stpz = delz

gamma = gamma + stpz * k2gg
theta = theta + stpz * k2pp

k3gg = k3gg - k2gg
k3pp = k3pp - k2pp

k2gg = k2gg * 2
k2pp = k2pp * 2

call fel_ode (gamma, theta, xks, xku, btpar, rpart, k2gg, k2pp)

gamma = gamma + stpz * (k3gg + k2gg / 6.0_rp)
theta = theta + stpz * (k3pp + k2pp / 6.0_rp)

end subroutine fel_runge_kutta

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_ode (tgam, tthet, xks, xku, btpar, rpart, k2gg, k2pp)
!
! The longitudinal equations of motion, BeamSolver::ODE (BeamSolver.cpp:144-163),
! fundamental only, ez = 0.
!-

subroutine fel_ode (tgam, tthet, xks, xku, btpar, rpart, k2gg, k2pp)

real(rp) tgam, tthet, xks, xku, btpar, k2gg, k2pp
complex(rp) rpart, ctmp
real(rp) ztemp1, btper0, btpar0

!

ztemp1 = -2.0_rp / xks
ctmp = rpart * cmplx(cos(tthet), -sin(tthet), rp)

btper0 = btpar + ztemp1 * real(ctmp, rp)
btpar0 = sqrt(1 - btper0 / (tgam * tgam))

k2pp = k2pp + xks * (1 - 1/btpar0) + xku
k2gg = k2gg + aimag(ctmp) / btpar0 / tgam       ! - ez, which is zero here.

end subroutine fel_ode

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_grid_weights (wf, x, y, ix, iy, wx, wy, on_grid)
!
! Routine to map a particle position to its lower-left grid cell corner and bilinear
! weights. Transcribed from Field::getLLGridpoint; wx is the weight of the LOWER x point.
! A particle outside |x|,|y| < gridmax neither feels the field nor radiates into it.
!-

subroutine fel_grid_weights (wf, x, y, ix, iy, wx, wy, on_grid)

type (wavefront_struct) wf
real(rp) x, y, wx, wy, gridmax
integer ix, iy, ngrid_arr(3), ngrid
logical on_grid

!

ngrid_arr = wavefront_shape(wf)
ngrid = ngrid_arr(1)
gridmax = (ngrid - 1) * wf%dx / 2       ! Genesis's gridmax; dgrid = 2*gridmax/(ngrid-1).

if (x > -gridmax .and. x < gridmax .and. y > -gridmax .and. y < gridmax) then
  wx = (x + gridmax) / wf%dx
  wy = (y + gridmax) / wf%dy
  ix = int(floor(wx))
  iy = int(floor(wy))
  wx = 1 + floor(wx) - wx
  wy = 1 + floor(wy) - wy
  ix = ix + 1                           ! To Fortran 1-based.
  iy = iy + 1
  on_grid = .true.
else
  on_grid = .false.
endif

end subroutine fel_grid_weights

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_field_step (und, beam, sl, wf, delz, err_flag)
!
! Routine to advance the field one step: build the source from the (already pushed)
! particles, propagate exp(K2 delz) in transverse Fourier space, add the source.
! Transcribed from FieldSolverFFT::advance and FFT, unfiltered path.
!
! The source scale, weighted. Genesis (FieldSolverFFT.cpp:21-22):
!
!   scl = fc(harm)*vacimp*current*ks*delz / (4*eev*npart*dgrid^2)
!
! carries current/npart per macroparticle. With per-particle weights that becomes
! c*w_j/slicelength (identical for uniform w_j = I*slicelength/(c*n)):
!
!   scl_w = fc*vacimp*ks*delz*c / (4*eev*dgrid^2*slicelength);  per particle scl_w*w_j
!
! and per particle part = sqrt(faw2(x,y))*scl_w*w_j/gamma, deposited as
! (sin theta + i cos theta)*part with the bilinear weights, added times 2 in real space
! after the transform pair (FieldSolverFFT.cpp:111).
!-

subroutine fel_field_step (und, beam, sl, wf, delz, err_flag)

type (fel_und_struct) und
type (fel_beam_struct) beam
type (fel_slice_struct) sl
type (wavefront_struct), target :: wf
real(rp) delz
logical err_flag

real(rp) xks, dgrid, scl_w, part, theta, beta, wx, wy, gam
complex(rp) cpart
integer ip, ix, iy, ngrid_arr(3), ngrid
logical on_grid, err

!

err_flag = .true.

ngrid_arr = wavefront_shape(wf)
ngrid = ngrid_arr(1)
xks = twopi / wf%wavelength
dgrid = wf%dx

call fel_field_cache_check (ngrid, dgrid, xks)

fel_crsource = 0

scl_w = fel_und_coupling(und, 1) * fel_vacimp * xks * delz * c_light
scl_w = scl_w / (4 * fel_eev * dgrid * dgrid * beam%slicelength)

if (scl_w /= 0) then
  do ip = 1, sl%n
    call fel_grid_weights (wf, sl%x(ip), sl%y(ip), ix, iy, wx, wy, on_grid)
    if (.not. on_grid) cycle

    gam = fel_gamma_of(beam%p0_mc, sl%pz(ip))
    beta = fel_beta_of(beam%p0_mc, sl%pz(ip))
    theta = beam%phi0 + xks * sl%z(ip) / beta          ! harm = 1: theta not scaled.

    part = sqrt(faw2(und, sl%x(ip), sl%y(ip))) * scl_w * sl%weight(ip) / gam
    cpart = cmplx(sin(theta), cos(theta), rp) * part

    fel_crsource(ix,   iy)   = fel_crsource(ix,   iy)   + (wx * wy) * cpart
    fel_crsource(ix+1, iy)   = fel_crsource(ix+1, iy)   + ((1-wx) * wy) * cpart
    fel_crsource(ix,   iy+1) = fel_crsource(ix,   iy+1) + (wx * (1-wy)) * cpart
    fel_crsource(ix+1, iy+1) = fel_crsource(ix+1, iy+1) + ((1-wx) * (1-wy)) * cpart
  enddo
endif

! Propagate and add: FFT, multiply exp(K2 delz), inverse FFT, normalize, plus
! 2*crsource in real space.

call wavefront_fft2 (wf%Ex(:,:,1), wf_fft_forward$, err);  if (err) return
wf%Ex(:,:,1) = wf%Ex(:,:,1) * exp(fel_k2 * delz)
call wavefront_fft2 (wf%Ex(:,:,1), wf_fft_backward$, err);  if (err) return
wf%Ex(:,:,1) = wf%Ex(:,:,1) / real(ngrid*ngrid, rp) + 2 * fel_crsource

err_flag = .false.

end subroutine fel_field_step

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_field_cache_check (ngrid, dgrid, ks)
!
! Routine to (re)build the K2 kernel and source accumulator when the grid changes,
! exactly as FieldSolverFFT::init builds them: dk = twopi/(ngrid*dgrid), integer offsets
! from the grid center, fftshift index mapping, K2 = -i*(dx^2+dy^2)*dk^2/(2*ks).
!-

subroutine fel_field_cache_check (ngrid, dgrid, ks)

integer ngrid
real(rp) dgrid, ks
real(rp) dk, shift, dx, dy
integer ix, iy, iix, iiy

!

if (ngrid == fel_cache_ngrid .and. dgrid == fel_cache_dgrid .and. ks == fel_cache_ks) return

if (allocated(fel_k2)) deallocate(fel_k2, fel_crsource)
allocate (fel_k2(ngrid, ngrid), fel_crsource(ngrid, ngrid))

dk = twopi / (ngrid * dgrid)
shift = -0.5_rp * (ngrid - 1)

do iy = 0, ngrid-1
  dy = iy + shift
  do ix = 0, ngrid-1
    dx = ix + shift
    iiy = mod(iy + (ngrid+1)/2, ngrid)
    iix = mod(ix + (ngrid+1)/2, ngrid)
    fel_k2(iix+1, iiy+1) = cmplx(0.0_rp, -(dx*dx + dy*dy) * dk * dk / 2 / ks, rp)
  enddo
enddo

fel_cache_ngrid = ngrid
fel_cache_dgrid = dgrid
fel_cache_ks = ks

end subroutine fel_field_cache_check

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_field_diag (wf, power, on_axis_intensity)
!
! Routine to compute the field diagnostics in physical units from a field in Genesis
! internal units, matching DiagField::calc: power = sum|u|^2 * (dgrid*eev/ks)^2/vacimp
! [W]; on-axis intensity = |u(center)|^2 * eev^2/(ks^2 vacimp) [W/m^2]. Accumulation
! order matches Genesis's: y outer, x inner.
!-

subroutine fel_field_diag (wf, power, on_axis_intensity)

type (wavefront_struct) wf
real(rp) power, on_axis_intensity
real(rp) ks, scl, wei
integer ix, iy, ngrid_arr(3), ngrid, ic

!

ngrid_arr = wavefront_shape(wf)
ngrid = ngrid_arr(1)
ks = twopi / wf%wavelength
scl = wf%dx * fel_eev / ks

power = 0
do iy = 1, ngrid
  do ix = 1, ngrid
    wei = real(wf%Ex(ix,iy,1), rp)**2 + aimag(wf%Ex(ix,iy,1))**2
    power = power + wei
  enddo
enddo
power = power * scl * scl / fel_vacimp

ic = ngrid/2 + 1
on_axis_intensity = (real(wf%Ex(ic,ic,1), rp)**2 + aimag(wf%Ex(ic,ic,1))**2) &
                    * fel_eev**2 / (ks**2 * fel_vacimp)

end subroutine fel_field_diag

end module fel_track_mod
