!+
! Module fel_track_mod
!
! The steady-state FEL step: period-averaged particle push, source deposition and FFT field
! solve inside an undulator segment. Transcribed from Genesis 1.3 Version 4, which GPL
! permits; the sources transcribed from are, by routine:
!
!   fel_transverse_track   <- TrackBeam::track, applyDrift, applyFQuad, applyDQuad
!                             (src/Core/TrackBeam.cpp:8-116)
!   fel_advance            <- BeamSolver::advance, RungeKutta, ODE
!                             (src/Core/BeamSolver.cpp:13-163)
!   fel_field_step         <- FieldSolverFFT::advance, FFT, init (unfiltered path)
!                             (src/Core/FieldSolverFFT.cpp)
!   fel_und_coupling       <- Undulator::fc (src/Core/Undulator.cpp:141-167)
!   faw, faw2              <- Undulator::faw, faw2 (src/Core/Undulator.cpp:170-182)
!
! and the step composition (transverse half step, longitudinal advance, transverse half
! step, then the field step) is Beam::track (src/Core/Beam.cpp:97) followed by Gencore's
! step 4 (src/Core/Gencore.cpp:279-289). Variable names are kept from Genesis wherever
! reasonable (btpar, k2gg, k2pp, rpart, crsource ...) so the transcription can be audited
! line against line.
!
! Units. The field inside this module is in Genesis's internal units, NOT volts per meter.
! The relation, from the power diagnostic (Diagnostic.cpp:801-980) and the dump scale
! (writeFieldHDF5.cpp:70), is
!
!   dfl [sqrt(W)] = u * dgrid * eev / (ks * sqrt(vacimp)),    E [V/m] = dfl * sqrt(2 Z0) / dgrid
!
! with dgrid the grid spacing. Working in internal units means every coefficient in the
! dynamics -- the source scale, the coupling, the energy exchange -- is computed from the
! same numbers Genesis computes it from, so the two implementations can be compared without
! a units conversion inside the loop. Conversion to and from the wavefront_struct's V/m
! happens once at the program boundary; see fel_field_unit_scale.
!
! What is deliberately absent, per the deliverable: no time dependence (single slice, no
! slippage), no space charge (Genesis's EFieldSolver zeroes ez when not enabled, so the
! transcription drops the term), no wakes, no incoherent synchrotron radiation, no
! harmonics (the fundamental only; fel_und_coupling still implements the planar Bessel
! factors so it is not wrong for a planar undulator), no source filter (the deck default;
! the filtered path in FieldSolverFFT needs one more forward transform and is not
! transcribed), and no orbit or field errors.
!-

module fel_track_mod

use fel_beam_mod
use wavefront_mod

implicit none

!+
! Struct fel_und_struct
!
! One undulator segment, constant parameters along it. This corresponds to one contiguous
! run of aw > 0 steps in Genesis's unrolled lattice (Lattice::unrollLattice): the segment
! length is divided into nstep equal steps of dz = l / nstep with
! nstep = round(l / delz_target), delz_target being the &setup delz.
!
! kx and ky here carry Genesis's unroll scaling (Lattice.cpp:412-413): the deck values
! (0.5 and 0.5 for a helical device by default, 0 and 1 for a planar one) multiplied by
! ku^2. faw and faw2 below use them directly, as Genesis does.
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

! Cached kernel state for fel_field_step, mirroring FieldSolverFFT's init: the propagation
! kernel K2 and the source accumulator are built once per (ngrid, dgrid, ks) and reused.
! Module state, single threaded, same caveat as the wavefront plan cache.

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
! Routine to return the factor converting a field in V/m (the wavefront_struct convention)
! to Genesis internal units. Multiply by it on entry to FEL tracking, divide on exit.
!
! Composition of the two exact relations in the module header:
!
!   u = E * dgrid / sqrt(2 Z0) * ks * sqrt(vacimp) / (dgrid * eev)
!     = E * ks * sqrt(vacimp) / (sqrt(2 Z0) * eev)
!
! Z0 here is Bmad's mu_0_vac * c_light because that is what the wavefront Genesis4 reader
! used to make V/m; vacimp and eev are Genesis's constants because that is what the
! internal units are defined by. The dgrid factors cancel exactly.
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
! Routine to return the coupling factor fc(h). Transcribed from Undulator::fc
! (src/Core/Undulator.cpp:141-167): a helical device couples the fundamental with strength
! aw and nothing else; a planar one couples odd harmonics through the Bessel factor
! JJ = J_{(h-1)/2}(xi) - J_{(h+1)/2}(xi), xi = h/2 * aw^2/(1+aw^2), with the (-1)^((h-1)/2)
! sign, and even harmonics not at all (that is a polarisation limitation, not physics; see
! the design brief section 5).
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
! gather. Undulator::faw (src/Core/Undulator.cpp:178-182).
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
! Square of the transverse dependence, used in the source deposition.
! Undulator::faw2 (src/Core/Undulator.cpp:170-174).
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
! Subroutine fel_track_und_step (und, sl, wf, gamma_ref, err_flag)
!
! Routine to advance one integration step of length und%dz inside an undulator segment.
! The composition is Genesis's, in Genesis's order:
!
!   1. Transverse half step (TrackBeam::track over dz/2, natural focusing only).
!   2. Longitudinal RK4 advance of (theta, gamma) over dz, gathering the field at the
!      particle positions of this moment and holding it fixed through the stages.
!   3. Transverse half step.
!   4. Field step: deposit the source from the pushed particles, then one exp(K2 dz)
!      propagation in transverse Fourier space, then add the source.
!
! The gather in step 2 sees the field before the solve of step 4 -- first order operator
! splitting, exactly as Genesis does it (Gencore.cpp: beam step 2 precedes field step 4).
!
! Input:
!   und         -- fel_und_struct: Segment parameters.
!   sl          -- fel_slice_struct: The slice. Steady state: one slice, periodic.
!   wf          -- wavefront_struct: The field, in Genesis internal units. One z slice.
!   gamma_ref   -- real(rp): Reference gamma of the run.
!
! Output:
!   sl, wf      -- Advanced one step.
!   err_flag    -- logical: Set True on error.
!-

subroutine fel_track_und_step (und, sl, wf, gamma_ref, err_flag)

type (fel_und_struct) und
type (fel_slice_struct) sl
type (wavefront_struct) wf
real(rp) gamma_ref
logical err_flag

!

call fel_transverse_track (und, sl, und%dz/2, gamma_ref)
call fel_advance (und, sl, wf, und%dz, gamma_ref)
call fel_transverse_track (und, sl, und%dz/2, gamma_ref)
call fel_field_step (und, sl, wf, und%dz, err_flag)

end subroutine fel_track_und_step

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_track_interlude_genesis (qf, length, sl, wf, gamma_ref, err_flag)
!
! Routine to advance one field-free interlude element -- a drift or a quadrupole -- the
! way Genesis does it, as one integration step of length `length`:
!
!   1. Transverse half step: drift or thick quad map over length/2, chromatic through the
!      per particle gammaz, focusing qquad = qf*gamma0 (TrackBeam::track with aw = 0).
!   2. theta advance over the full length with dtheta/dz = ks*(1 - 1/beta_z) + xku,
!      xku = ks*0.5/gamma0/gamma0 (the BeamSolver drift special case,
!      BeamSolver.cpp:35-38), beta_z from the px, py of this mid-element moment. The RK4
!      collapses to a single exact step because with no field and no space charge the
!      slope is theta independent; the collapse is applied here directly.
!   3. Transverse half step.
!   4. Field diffraction over the element: the FFT propagation with zero source, exactly
!      Genesis's unconditional field track.
!
! This routine exists for one purpose: it is the Genesis interlude model, transcribed, so
! that running the full lattice with it isolates what the Bmad seam changes. The seam
! (track1_bunch plus the exact theta mapping in fel_bunch_to_slice) differs from this in
! one measurable way -- Genesis samples the path length term px^2+py^2 once at mid
! element, Bmad integrates it through the quad map -- and comparing the two runs prices
! that difference; see the benchmark README. The production configuration is the seam.
!
! Input:
!   qf          -- real(rp): Quadrupole strength k1 [1/m^2], Genesis's deck convention,
!                    zero for a drift.
!   length      -- real(rp): Element length [m].
!   sl          -- fel_slice_struct: The slice.
!   wf          -- wavefront_struct: The field, Genesis internal units.
!   gamma_ref   -- real(rp): Reference gamma.
!
! Output:
!   sl, wf      -- Advanced through the element.
!   err_flag    -- logical: Set True on error.
!-

subroutine fel_track_interlude_genesis (qf, length, sl, wf, gamma_ref, err_flag)

type (fel_slice_struct) sl
type (wavefront_struct), target :: wf
real(rp) qf, length, gamma_ref
logical err_flag

type (fel_und_struct) und0
real(rp) xks, xku, qquad, gammaz, btpar, btpar0, slope
integer ip, ngrid_arr(3), ngrid
logical err

!

err_flag = .true.

xks = twopi / wf%wavelength
xku = xks * 0.5_rp / gamma_ref / gamma_ref     ! Genesis's division order, kept.
qquad = qf * gamma_ref

! Step 1: transverse half step. aw = 0, so gammaz = sqrt(gamma^2 - 1 - px^2 - py^2)
! and the only focusing is the quadrupole's.

do ip = 1, sl%n
  gammaz = sqrt(sl%gamma(ip)**2 - 1 - sl%px(ip)**2 - sl%py(ip)**2)
  call fel_apply_focus (length/2, qquad, sl%x(ip), sl%px(ip), gammaz)
  call fel_apply_focus (length/2, -qquad, sl%y(ip), sl%py(ip), gammaz)
enddo

! Step 2: theta advance, the collapsed RK4. btpar = 1 + px^2 + py^2 (aw = 0), no energy
! change (no field coupling, no space charge).

do ip = 1, sl%n
  btpar = 1 + sl%px(ip)**2 + sl%py(ip)**2
  btpar0 = sqrt(1 - btpar / (sl%gamma(ip) * sl%gamma(ip)))
  slope = xks * (1 - 1/btpar0) + xku
  sl%theta(ip) = sl%theta(ip) + length * slope
enddo

! Step 3: second transverse half step.

do ip = 1, sl%n
  gammaz = sqrt(sl%gamma(ip)**2 - 1 - sl%px(ip)**2 - sl%py(ip)**2)
  call fel_apply_focus (length/2, qquad, sl%x(ip), sl%px(ip), gammaz)
  call fel_apply_focus (length/2, -qquad, sl%y(ip), sl%py(ip), gammaz)
enddo

! Step 4: field diffraction, Genesis's kernel, zero source (aw = 0 makes the coupling and
! hence the source scale exactly zero, matching Genesis skipping the deposition).

und0%aw = 0
call fel_field_step (und0, sl, wf, length, err);  if (err) return

err_flag = .false.

end subroutine fel_track_interlude_genesis

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_transverse_track (und, sl, delz, gamma_ref)
!
! Routine to advance the transverse coordinates over delz inside an undulator: drift plus
! the undulator's natural focusing. Transcribed from TrackBeam::track with the quadrupole,
! corrector and chicane branches dropped -- in this deliverable quadrupoles never overlap
! undulator fields (they are separate elements handled by Bmad) and there are no correctors
! or chicanes.
!
! Genesis's effective focusing strengths (TrackBeam.cpp:23-31):
!
!   betpar0 = sqrt(1 - (1+aw^2)/gamma0^2)
!   qnatx   = kx*aw^2/(gamma0*betpar0),  qnaty = ky*aw^2/(gamma0*betpar0)
!
! applied per particle with gammaz = sqrt(gamma^2 - 1 - aw^2 - px^2 - py^2) through the
! exact cos/sin thick lens maps.
!-

subroutine fel_transverse_track (und, sl, delz, gamma_ref)

type (fel_und_struct) und
type (fel_slice_struct) sl
real(rp) delz, gamma_ref
real(rp) betpar0, qnatx, qnaty, qx, qy, aw2, gammaz
integer ip

!

aw2 = und%aw**2
betpar0 = sqrt(1 - (1 + aw2)/gamma_ref**2)

qnatx = und%kx * aw2 / gamma_ref / betpar0
qnaty = und%ky * aw2 / gamma_ref / betpar0
qx = qnatx      ! No quadrupole component inside undulator segments in this deliverable.
qy = qnaty

do ip = 1, sl%n
  gammaz = sqrt(sl%gamma(ip)**2 - 1 - aw2 - sl%px(ip)**2 - sl%py(ip)**2)
  call fel_apply_focus (delz, qx, sl%x(ip), sl%px(ip), gammaz)
  call fel_apply_focus (delz, qy, sl%y(ip), sl%py(ip), gammaz)
enddo

end subroutine fel_transverse_track

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_apply_focus (delz, qf, x, px, gammaz)
!
! Routine to apply one transverse plane's map over delz: TrackBeam::applyDrift for qf = 0,
! applyFQuad for qf > 0, applyDQuad for qf < 0 (TrackBeam.cpp:84-116), with the offset dx
! of those routines zero (no misaligned focusing here).
!-

subroutine fel_apply_focus (delz, qf, x, px, gammaz)

real(rp) delz, qf, x, px, gammaz
real(rp) foc, omg, a1, a2, a3, xtmp

!

if (qf == 0) then
  x = x + px * delz / gammaz
elseif (qf > 0) then
  foc = sqrt(qf/gammaz)
  omg = foc * delz
  a1 = cos(omg)
  a2 = sin(omg)/foc
  a3 = -a2 * foc * foc
  xtmp = x
  x  = a1 * xtmp + a2 * px / gammaz
  px = a3 * xtmp * gammaz + a1 * px
else
  foc = sqrt(-qf/gammaz)
  omg = foc * delz
  a1 = cosh(omg)
  a2 = sinh(omg)/foc
  a3 = a2 * foc * foc
  xtmp = x
  x  = a1 * xtmp + a2 * px / gammaz
  px = a3 * xtmp * gammaz + a1 * px
endif

end subroutine fel_apply_focus

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_advance (und, sl, wf, delz, gamma_ref)
!
! Routine to advance (theta, gamma) of every particle over delz by fourth order Runge
! Kutta. Transcribed from BeamSolver::advance, RungeKutta and ODE (BeamSolver.cpp:13-163),
! fundamental field only, no space charge (ez = 0), autophase = 0.
!
! Per particle: gather the field at (x, y) by bilinear interpolation (the same four points
! and weights as the deposition), form
!
!   rpart = (fc(1)/ks) * faw(x,y) * conj(E_gathered)
!
! and integrate
!
!   dtheta/dz = ks*(1 - 1/btpar0) + ku
!   dgamma/dz = Im(rpart * exp(-i theta)) / (btpar0 * gamma)
!
! with btper0 = 1 + px^2 + py^2 + (aw*faw)^2 - (2/ks)*Re(rpart * exp(-i theta)) and
! btpar0 = sqrt(1 - btper0/gamma^2). rpart, px, py and faw are gathered once and held fixed
! through the four stages; only theta and gamma vary. The RK4 bookkeeping is pushp's from
! the old Fortran source, kept verbatim so every add and multiply matches.
!-

subroutine fel_advance (und, sl, wf, delz, gamma_ref)

type (fel_und_struct) und
type (fel_slice_struct) sl
type (wavefront_struct) wf
real(rp) delz, gamma_ref

real(rp) xks, xku, aw, rtmp, awloc, btpar, gamma, theta, wx, wy
complex(rp) cpart, rpart
integer ip, ix, iy
logical on_grid

! Genesis: xks from the field, xku from the undulator; in a drift xku would be
! ks/(2 gamma_ref^2) (BeamSolver.cpp:35-38) but this routine is only called inside
! undulator segments.

xks = twopi / wf%wavelength
xku = und%ku
aw = und%aw
rtmp = fel_und_coupling(und, 1) / xks     ! fc(harm)/field->xks, harm = 1.

do ip = 1, sl%n
  gamma = sl%gamma(ip)
  theta = sl%theta(ip)                    ! autophase = 0: no phase shifters here.
  awloc = faw(und, sl%x(ip), sl%y(ip))
  btpar = 1 + sl%px(ip)**2 + sl%py(ip)**2 + aw*aw*awloc*awloc

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

  sl%gamma(ip) = gamma
  sl%theta(ip) = theta
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
! weights. Transcribed from Field::getLLGridpoint (src/Core/Field.cpp:160-175); the
! translation from Genesis's flat 0-based idx to Fortran (ix, iy) 1-based pairs is
! idx = (ix-1) + (iy-1)*ngrid, and Genesis's four gather/deposit points idx, idx+1,
! idx+ngrid, idx+ngrid+1 are (ix,iy), (ix+1,iy), (ix,iy+1), (ix+1,iy+1).
!
! wx is the weight of the LOWER x point (Genesis: wx = 1 + floor(w) - w).
!
! Output:
!   on_grid   -- logical: False if the particle is outside |x|,|y| < gridmax, in which
!                  case it neither feels the field nor radiates into it.
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
! Subroutine fel_field_step (und, sl, wf, delz, err_flag)
!
! Routine to advance the field one step: build the source from the (already pushed)
! particles, propagate the field exp(K2 delz) in transverse Fourier space, add the source.
! Transcribed from FieldSolverFFT::advance and FFT, unfiltered path
! (src/Core/FieldSolverFFT.cpp:9-61 and 91-113).
!
! The source scale (FieldSolverFFT.cpp:21-22):
!
!   scl = fc(harm) * vacimp * current * ks * delz / (4 * eev * npart * dgrid^2)
!
! and per particle part = sqrt(faw2(x,y)) * scl / gamma deposited as
! (sin(theta) + i cos(theta)) * part with the bilinear weights. The source is added in
! real space after the transform pair, times 2, since IFFT(FFT(s))/n^2 = s
! (FieldSolverFFT.cpp:111).
!
! The kernel K2 = -i*(kx^2+ky^2)/(2 ks) is built as Genesis builds it -- integer offsets
! from the grid center times dk, fftshift index mapping -- so the phases match Genesis's
! association exactly (FieldSolverFFT.cpp:131-145).
!-

subroutine fel_field_step (und, sl, wf, delz, err_flag)

type (fel_und_struct) und
type (fel_slice_struct) sl
type (wavefront_struct), target :: wf
real(rp) delz
logical err_flag

real(rp) xks, dgrid, scl, part, theta, wx, wy, gam
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

! Source construction. In Genesis this is gated on inUndulator && isEnabled && odd
! harmonic; this routine is only called inside undulator segments with the fundamental.

fel_crsource = 0

scl = fel_und_coupling(und, 1) * fel_vacimp * sl%current * xks * delz
scl = scl / (4 * fel_eev * real(sl%n, rp) * dgrid * dgrid)

do ip = 1, sl%n
  theta = sl%theta(ip)                 ! harm = 1: theta not scaled.
  gam = sl%gamma(ip)

  call fel_grid_weights (wf, sl%x(ip), sl%y(ip), ix, iy, wx, wy, on_grid)
  if (.not. on_grid) cycle

  part = sqrt(faw2(und, sl%x(ip), sl%y(ip))) * scl / gam
  cpart = cmplx(sin(theta), cos(theta), rp) * part

  fel_crsource(ix,   iy)   = fel_crsource(ix,   iy)   + (wx * wy) * cpart
  fel_crsource(ix+1, iy)   = fel_crsource(ix+1, iy)   + ((1-wx) * wy) * cpart
  fel_crsource(ix,   iy+1) = fel_crsource(ix,   iy+1) + (wx * (1-wy)) * cpart
  fel_crsource(ix+1, iy+1) = fel_crsource(ix+1, iy+1) + ((1-wx) * (1-wy)) * cpart
enddo

! Propagate and add the source: FFT, multiply exp(K2 delz), inverse FFT, normalize,
! plus 2*crsource in real space.

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
! Routine to (re)build the K2 kernel and source accumulator when the grid changes.
! The kernel construction is FieldSolverFFT::init (FieldSolverFFT.cpp:123-145): dk =
! twopi/(ngrid*dgrid), offsets dx, dy = index + shift with shift = -(ngrid-1)/2, fftshift
! index mapping ii = (i + (ngrid+1)/2) mod ngrid, and K2 = -i*(dx^2+dy^2)*dk^2/(2*ks).
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
! internal units, matching DiagField::calc (Diagnostic.cpp:801-986):
!
!   power [W]              = sum |u|^2 * (dgrid*eev/ks)^2 / vacimp
!   on-axis intensity      = |u(center)|^2 * eev^2/(ks^2 vacimp)   [W/m^2]
!
! The accumulation order matches Genesis's: y outer, x inner.
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

ic = ngrid/2 + 1     ! Genesis: i = (ngrid/2)*ngrid + ngrid/2, 0-based -> center for odd ngrid.
on_axis_intensity = (real(wf%Ex(ic,ic,1), rp)**2 + aimag(wf%Ex(ic,ic,1))**2) &
                    * fel_eev**2 / (ks**2 * fel_vacimp)

end subroutine fel_field_diag

end module fel_track_mod
