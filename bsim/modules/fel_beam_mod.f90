!+
! Module fel_beam_mod
!
! The packed particle representation for FEL tracking, its Genesis 1.3 Version 4 particle
! dump I/O, its conversion to and from Bmad's coord_struct at element boundaries, and the
! per-slice beam diagnostics.
!
! Coordinates: Bmad's, exactly. Each slice carries structure-of-arrays copies of
! coord_struct%vec(1:6) plus a weight:
!
!   x, y      [m]
!   px, py    transverse momentum over the reference momentum, P_xy/p0
!   z         -beta*c*(t - t_ref)  [m], the Bmad longitudinal coordinate
!   pz        (p - p0)/p0
!   weight    macroparticle charge [C], mapping to coord_struct%charge
!
! so conversion at the seam is a plain copy. (z, pz) is the conjugate pair of Bmad's
! s-based Hamiltonian. The normalization reference p0 is carried on the beam (p0c in eV
! and p0_mc = p0c/m_e c^2 = gamma0*beta0 dimensionless) and asserted against each
! element's p0c at conversion time.
!
! The ponderomotive phase is derived, not stored. Genesis's per-particle theta is the sum
! of a common reference advance (the undulator's ku term, and the drift slippage term
! ks/(2 gamma0^2)) and a particle-specific lag. The common part lives in one scalar per
! beam, phi0, maintained by the tracker; the particle part is the Bmad z:
!
!   theta_j = phi0 - ks * tau_j,     tau_j = -z_j / beta_j = c*(t_j - t_ref)
!
! This is the reference-offset split that the design brief's section 8 identifies as the
! FP32-safe formulation, and it removes the brief's 6.4 hazard outright: z does not wrap,
! so there is no theta-wrap-plus-slice-index update to get wrong at slice migration; the
! slice index is derived from z when needed. (A future single precision GPU struct still
! needs per-slice re-referencing of z, since the phase needs ~1e-6 rad across ~1e5
! wavelengths of bunch; that is a device-struct choice the brief already anticipates.)
!
! Weights are carried from day one (brief section 5): every reduction here and in
! fel_track_mod is weighted, slice current is derived as I = c * sum(w) / slicelength, and
! N_eff = (sum w)^2 / sum w^2 is a per-slice diagnostic. A uniform-weight beam reproduces
! Genesis, which the benchmark gates.
!
! Why packed arrays at all (brief section 4.2): the FEL step advances every particle every
! internal step, and coord_struct is ~224 bytes against the ~56 needed. coord_struct
! appears only at element boundaries. The arrays are allocated to a capacity that may
! exceed the fill count n, so slice migration can later move particles without per-step
! reallocation (brief section 6.4).
!-

module fel_beam_mod

use bmad
use hdf5_interface

implicit none

! Genesis's own physical constants, transcribed from src/Main/GenMain.cpp:63-67 and used
! in every transcribed formula so the arithmetic matches Genesis digit for digit. They
! deliberately differ from Bmad's constants: vacimp is the impedance of free space
! truncated to 376.73 (Bmad's mu_0_vac * c_light is 376.7303134...), and eev is the
! electron rest mass in eV. Do not "fix" these; the comparison against Genesis depends on
! them. They enter only the field normalization and source scale, never the coordinate
! conversions, which use Bmad's own mass so the Bmad side stays self-consistent.

real(rp), parameter :: fel_vacimp = 376.73_rp        ! [Ohm]     GenMain.cpp:63
real(rp), parameter :: fel_eev    = 510998.95069_rp  ! [eV]      GenMain.cpp:67

!+
! Struct fel_slice_struct
!
! One beam slice in packed structure-of-arrays form, Bmad coordinates plus weight.
! Live particles are 1:n; the arrays may be larger.
!-

type fel_slice_struct
  real(rp), allocatable :: x(:)        ! [m]
  real(rp), allocatable :: px(:)       ! P_x/p0
  real(rp), allocatable :: y(:)        ! [m]
  real(rp), allocatable :: py(:)       ! P_y/p0
  real(rp), allocatable :: z(:)        ! -beta*c*(t - t_ref) [m]
  real(rp), allocatable :: pz(:)       ! (p - p0)/p0
  real(rp), allocatable :: weight(:)   ! Macroparticle charge [C]
  integer :: n = 0                     ! Fill count.
end type

!+
! Struct fel_beam_struct
!
! The whole beam: slices, the normalization reference, the common phase, and the dump
! metadata.
!-

type fel_beam_struct
  type (fel_slice_struct), allocatable :: slice(:)
  real(rp) :: p0c = 0              ! Reference momentum times c [eV].
  real(rp) :: p0_mc = 0            ! p0/(m_e c) = gamma0*beta0, dimensionless.
  real(rp) :: gamma0 = 0           ! Genesis's reference gamma (of the run, not of p0c).
  real(rp) :: phi0 = 0             ! Common ponderomotive reference phase [rad].
  real(rp) :: reflength = 0        ! Radiation wavelength [m]. 'slicelength' in the dump.
  real(rp) :: slicelength = 0      ! Slice spacing [m]. 'slicespacing' in the dump.
  real(rp) :: s0 = 0               ! Start of the time window [m]. 'refposition' in the dump.
  integer :: nbins = 0             ! Beamlet size at generation. Carried, not used here.
  logical :: one4one = .false.
end type

!+
! Struct fel_slice_diag_struct
!
! Per-slice beam diagnostics at one output position. Genesis's DiagBeam::calc definitions,
! weighted: plain weighted averages, bunching b = |sum w exp(i theta)| / sum w. N_eff per
! the brief's 6.2.
!-

type fel_slice_diag_struct
  real(rp) :: energy = 0           ! <gamma>
  real(rp) :: energyspread = 0     ! sqrt(|<gamma^2> - <gamma>^2|)
  real(rp) :: bunching = 0         ! |b(1)|
  real(rp) :: bunchingphase = 0    ! arg b(1)
  real(rp) :: xposition = 0, yposition = 0
  real(rp) :: xsize = 0, ysize = 0
  real(rp) :: pxposition = 0, pyposition = 0   ! <gamma beta_x>, <gamma beta_y> (Genesis px)
  real(rp) :: n_eff = 0            ! (sum w)^2 / sum w^2
  real(rp) :: current = 0          ! c * sum(w) / slicelength [A]
end type

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Elemental functions for the derived kinematic quantities. All exact algebra on the
! stored (pz, p0_mc):
!
!   P     = p0_mc * (1 + pz)            total momentum / m_e c
!   gamma = sqrt(P^2 + 1)
!   beta  = P / gamma
!   tau   = -z / beta                   c*(t - t_ref) [m]
!-

elemental function fel_gamma_of (p0_mc, pz) result (gamma)
real(rp), intent(in) :: p0_mc, pz
real(rp) gamma
gamma = sqrt((p0_mc * (1 + pz))**2 + 1)
end function

elemental function fel_beta_of (p0_mc, pz) result (beta)
real(rp), intent(in) :: p0_mc, pz
real(rp) beta, p_mc
p_mc = p0_mc * (1 + pz)
beta = p_mc / sqrt(p_mc**2 + 1)
end function

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function fel_theta (beam, sl, ip, ks) result (theta)
!
! The ponderomotive phase of particle ip: theta = phi0 - ks*tau, tau = -z/beta.
! See the module header for the split.
!-

function fel_theta (beam, sl, ip, ks) result (theta)

type (fel_beam_struct) beam
type (fel_slice_struct) sl
integer ip
real(rp) ks, theta, beta

!

beta = fel_beta_of(beam%p0_mc, sl%pz(ip))
theta = beam%phi0 + ks * sl%z(ip) / beta      ! phi0 - ks*(-z/beta)

end function fel_theta

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_read_genesis4_beam (beam, file_name, gamma0, err_flag)
!
! Routine to read a Genesis 1.3 Version 4 particle dump (.par.h5) into a packed beam,
! converting Genesis coordinates (x, y in m; px, py = gamma*beta; theta; gamma) to Bmad's.
!
! The conversion, exact:
!
!   px_bmad = px_genesis / p0_mc                    (p0_mc = sqrt(gamma0^2 - 1))
!   pz      = (sqrt(gamma^2 - 1) - p0_mc) / p0_mc
!   tau     = -theta / ks,   z = -beta * tau        (phi0 starts at zero)
!   weight  = I * slicelength / (c * n)  [C]        (uniform; the dump carries no weights)
!
! ks comes from the dump's own wavelength. p0c is set from gamma0 and Bmad's electron
! mass, which is what makes the Bmad-side tracking see exactly this normalization.
!
! Input:
!   file_name   -- character(*): File to read.
!   gamma0      -- real(rp): Genesis's reference gamma for the run.
!
! Output:
!   beam        -- fel_beam_struct: Beam read from the file.
!   err_flag    -- logical: Set True on error, False otherwise.
!-

subroutine fel_read_genesis4_beam (beam, file_name, gamma0, err_flag)

type (fel_beam_struct), target :: beam
type (fel_slice_struct), pointer :: sl
integer(hid_t) f_id, g_id
integer is, ip, n_slice, np, ivec(1), h5_err
real(rp) rvec(1), gamma0, ks, current, w_uniform, gam, p_mc, beta, tau
real(rp), allocatable :: work_gamma(:), work_theta(:)
logical err_flag, err
character(*) file_name
character(*), parameter :: r_name = 'fel_read_genesis4_beam'
character(20) group_name
type (hdf5_info_struct) info

!

err_flag = .true.

if (gamma0 <= 1) then
  call out_io (s_error$, r_name, 'gamma0 MUST EXCEED 1. GOT: \es12.4\ ', r_array = [gamma0])
  return
endif

call hdf5_open_file (file_name, 'READ', f_id, err);  if (err) return

call hdf5_read_dataset_int (f_id, 'slicecount', ivec, err, 'slicecount');  if (err) return
n_slice = ivec(1)
if (n_slice < 1) then
  call out_io (s_error$, r_name, 'FILE HAS A NON-POSITIVE slicecount: \i0\ ', i_array = [n_slice])
  return
endif

call hdf5_read_dataset_int (f_id, 'beamletsize', ivec, err, 'beamletsize');  if (err) return
beam%nbins = ivec(1)

call hdf5_read_dataset_int (f_id, 'one4one', ivec, err, 'one4one');  if (err) return
beam%one4one = (ivec(1) /= 0)

call hdf5_read_dataset_real (f_id, 'slicelength', rvec, err, 'slicelength');  if (err) return
beam%reflength = rvec(1)

call hdf5_read_dataset_real (f_id, 'slicespacing', rvec, err, 'slicespacing');  if (err) return
beam%slicelength = rvec(1)

call hdf5_read_dataset_real (f_id, 'refposition', rvec, err, 'refposition');  if (err) return
beam%s0 = rvec(1)

beam%gamma0 = gamma0
beam%p0_mc = sqrt(gamma0**2 - 1)
beam%p0c = beam%p0_mc * mass_of(electron$)
beam%phi0 = 0
ks = twopi / beam%reflength

if (allocated(beam%slice)) deallocate(beam%slice)
allocate (beam%slice(n_slice))

do is = 1, n_slice
  sl => beam%slice(is)
  write (group_name, '(a, i0.6)') 'slice', is
  g_id = hdf5_open_group (f_id, trim(group_name), err, .true.);  if (err) return

  ! Particle count from the actual dataset extent (FINDINGS.md 5.1: H5LT reads have no
  ! buffer bound; extents are checked, never assumed).

  info = hdf5_object_info (g_id, 'gamma', err, .true.);  if (err) return
  np = int(info%data_dim(1))

  call fel_slice_reallocate (sl, np)
  sl%n = np

  call hdf5_read_dataset_real (g_id, 'current', rvec, err, trim(group_name) // '/current');  if (err) return
  current = rvec(1)

  allocate (work_gamma(np), work_theta(np))
  call hdf5_read_dataset_real (g_id, 'gamma', work_gamma, err, trim(group_name) // '/gamma');  if (err) return
  call hdf5_read_dataset_real (g_id, 'theta', work_theta, err, trim(group_name) // '/theta');  if (err) return
  call hdf5_read_dataset_real (g_id, 'x',  sl%x(1:np),  err, trim(group_name) // '/x');   if (err) return
  call hdf5_read_dataset_real (g_id, 'y',  sl%y(1:np),  err, trim(group_name) // '/y');   if (err) return
  call hdf5_read_dataset_real (g_id, 'px', sl%px(1:np), err, trim(group_name) // '/px');  if (err) return
  call hdf5_read_dataset_real (g_id, 'py', sl%py(1:np), err, trim(group_name) // '/py');  if (err) return

  ! Convert in place: px, py from gamma*beta to P/p0; (theta, gamma) to (z, pz).

  w_uniform = 0
  if (np > 0) w_uniform = current * beam%slicelength / (c_light * np)

  do ip = 1, np
    gam = work_gamma(ip)
    p_mc = sqrt(gam**2 - 1)
    beta = p_mc / gam

    sl%px(ip) = sl%px(ip) / beam%p0_mc
    sl%py(ip) = sl%py(ip) / beam%p0_mc
    sl%pz(ip) = (p_mc - beam%p0_mc) / beam%p0_mc

    tau = -work_theta(ip) / ks              ! theta = phi0 - ks*tau with phi0 = 0.
    sl%z(ip) = -beta * tau

    sl%weight(ip) = w_uniform
  enddo
  deallocate (work_gamma, work_theta)

  call H5Gclose_f (g_id, h5_err)
enddo

call H5Fclose_f (f_id, h5_err)

err_flag = .false.

end subroutine fel_read_genesis4_beam

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_write_genesis4_beam (beam, file_name, err_flag)
!
! Routine to write a packed beam as a Genesis 1.3 Version 4 particle dump, converting
! back from Bmad coordinates. Inverse of fel_read_genesis4_beam: theta = phi0 - ks*tau
! reconstructs Genesis's unwrapped theta including the accumulated common phase, gamma
! from pz, px py rescaled by p0_mc, current from the weights. The dump format carries no
! per-particle weight, so nonuniform weights do not survive a round trip; that is the
! format's limitation, not this representation's.
!-

subroutine fel_write_genesis4_beam (beam, file_name, err_flag)

type (fel_beam_struct), target :: beam
type (fel_slice_struct), pointer :: sl
integer(hid_t) f_id, g_id
integer is, ip, h5_err, one4one_int
real(rp) ks, gam, beta, p_mc
real(rp), allocatable :: work(:)
logical err_flag, err
character(*) file_name
character(*), parameter :: r_name = 'fel_write_genesis4_beam'
character(20) group_name

!

err_flag = .true.

if (.not. allocated(beam%slice)) then
  call out_io (s_error$, r_name, 'BEAM HAS NO SLICES.')
  return
endif

ks = twopi / beam%reflength

call hdf5_open_file (file_name, 'WRITE', f_id, err);  if (err) return

one4one_int = 0
if (beam%one4one) one4one_int = 1

call hdf5_write_dataset_real (f_id, 'slicelength',  [beam%reflength],   err);  if (err) return
call hdf5_write_dataset_real (f_id, 'slicespacing', [beam%slicelength], err);  if (err) return
call hdf5_write_dataset_real (f_id, 'refposition',  [beam%s0],          err);  if (err) return
call hdf5_write_dataset_int  (f_id, 'beamletsize',  [beam%nbins],       err);  if (err) return
call hdf5_write_dataset_int  (f_id, 'slicecount',   [size(beam%slice)], err);  if (err) return
call hdf5_write_dataset_int  (f_id, 'one4one',      [one4one_int],      err);  if (err) return

do is = 1, size(beam%slice)
  sl => beam%slice(is)
  write (group_name, '(a, i0.6)') 'slice', is
  call H5Gcreate_f (f_id, trim(group_name), g_id, h5_err)
  if (h5_err < 0) then
    call out_io (s_error$, r_name, 'CANNOT CREATE GROUP: ' // trim(group_name))
    return
  endif

  call hdf5_write_dataset_real (g_id, 'current', [c_light * sum(sl%weight(1:sl%n)) / beam%slicelength], err)
  if (err) return

  allocate (work(sl%n))

  do ip = 1, sl%n                                             ! gamma
    work(ip) = fel_gamma_of(beam%p0_mc, sl%pz(ip))
  enddo
  call hdf5_write_dataset_real (g_id, 'gamma', work, err);  if (err) return

  do ip = 1, sl%n                                             ! theta = phi0 + ks*z/beta
    work(ip) = beam%phi0 + ks * sl%z(ip) / fel_beta_of(beam%p0_mc, sl%pz(ip))
  enddo
  call hdf5_write_dataset_real (g_id, 'theta', work, err);  if (err) return

  call hdf5_write_dataset_real (g_id, 'x', sl%x(1:sl%n), err);  if (err) return
  call hdf5_write_dataset_real (g_id, 'y', sl%y(1:sl%n), err);  if (err) return

  work = sl%px(1:sl%n) * beam%p0_mc                           ! back to gamma*beta_x
  call hdf5_write_dataset_real (g_id, 'px', work, err);  if (err) return
  work = sl%py(1:sl%n) * beam%p0_mc
  call hdf5_write_dataset_real (g_id, 'py', work, err);  if (err) return

  deallocate (work)
  call H5Gclose_f (g_id, h5_err)
enddo

call H5Fclose_f (f_id, h5_err)

err_flag = .false.

end subroutine fel_write_genesis4_beam

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_slice_reallocate (sl, capacity)
!
! Routine to allocate a slice's arrays to at least the given capacity. Existing live
! particles 1:n are preserved. The fill count n is not changed.
!-

subroutine fel_slice_reallocate (sl, capacity)

type (fel_slice_struct) sl
integer capacity

!

if (allocated(sl%x)) then
  if (size(sl%x) >= capacity) return
  call re_allocate (sl%x, capacity)
  call re_allocate (sl%px, capacity)
  call re_allocate (sl%y, capacity)
  call re_allocate (sl%py, capacity)
  call re_allocate (sl%z, capacity)
  call re_allocate (sl%pz, capacity)
  call re_allocate (sl%weight, capacity)
else
  allocate (sl%x(capacity), sl%px(capacity), sl%y(capacity), sl%py(capacity), &
            sl%z(capacity), sl%pz(capacity), sl%weight(capacity))
endif

end subroutine fel_slice_reallocate

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_slice_to_bunch (beam, sl, ele, bunch, err_flag)
!
! Routine to convert a packed slice to a Bmad bunch_struct: plain copies, since the
! stored coordinates are coord_struct's. The element's p0c must match the beam's
! normalization; a mismatch is refused, not rescaled, since nothing in this deliverable
! changes the reference momentum.
!-

subroutine fel_slice_to_bunch (beam, sl, ele, bunch, err_flag)

type (fel_beam_struct) beam
type (fel_slice_struct) sl
type (ele_struct) ele
type (bunch_struct) bunch
real(rp) vec(6)
integer ip
logical err_flag
character(*), parameter :: r_name = 'fel_slice_to_bunch'

!

err_flag = .true.

if (abs(ele%value(p0c$) - beam%p0c) > 1e-10_rp * beam%p0c) then
  call out_io (s_error$, r_name, 'ELEMENT p0c \es20.12\ DOES NOT MATCH BEAM p0c \es20.12\ ', &
               'AT ELEMENT: ' // trim(ele%name), r_array = [ele%value(p0c$), beam%p0c])
  return
endif

if (allocated(bunch%particle)) then
  if (size(bunch%particle) /= sl%n) deallocate(bunch%particle)
endif
if (.not. allocated(bunch%particle)) allocate(bunch%particle(sl%n))

do ip = 1, sl%n
  vec = [sl%x(ip), sl%px(ip), sl%y(ip), sl%py(ip), sl%z(ip), sl%pz(ip)]

  ! init_coord derives beta, t and state consistently. shift_vec6 exists for elements
  ! whose reference momentum changes and must not touch vec(6) here.
  call init_coord (bunch%particle(ip), vec, ele, upstream_end$, electron$, shift_vec6 = .false.)
  bunch%particle(ip)%charge = sl%weight(ip)
enddo

bunch%n_live = sl%n
bunch%charge_live = sum(sl%weight(1:sl%n))

err_flag = .false.

end subroutine fel_slice_to_bunch

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_bunch_to_slice (bunch, ele, sl, err_flag)
!
! Routine to copy a tracked Bmad bunch back into the packed slice. Plain copies; the
! phase bookkeeping is one phi0 update per element, done by the caller
! (fel_phi0_advance), not here, because it is per beam and not per particle.
!-

subroutine fel_bunch_to_slice (bunch, ele, sl, err_flag)

type (bunch_struct) bunch
type (ele_struct) ele
type (fel_slice_struct) sl
integer ip
logical err_flag
character(*), parameter :: r_name = 'fel_bunch_to_slice'

!

err_flag = .true.

do ip = 1, sl%n
  if (bunch%particle(ip)%state /= alive$) then
    call out_io (s_error$, r_name, 'PARTICLE \i0\ LOST TRACKING THROUGH: ' // trim(ele%name), &
                                   i_array = [ip])
    return
  endif

  sl%x(ip)  = bunch%particle(ip)%vec(1)
  sl%px(ip) = bunch%particle(ip)%vec(2)
  sl%y(ip)  = bunch%particle(ip)%vec(3)
  sl%py(ip) = bunch%particle(ip)%vec(4)
  sl%z(ip)  = bunch%particle(ip)%vec(5)
  sl%pz(ip) = bunch%particle(ip)%vec(6)
  sl%weight(ip) = bunch%particle(ip)%charge
enddo

err_flag = .false.

end subroutine fel_bunch_to_slice

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function fel_phi0_rate (ks, ku_like, p0_mc) result (rate)
!
! Routine to return dphi0/ds, the advance rate of the common ponderomotive reference
! phase:
!
!   dphi0/ds = ku_like + ks*(1 - 1/beta0)
!
! where ku_like is the undulator wavenumber inside an undulator, or Genesis's drift
! surrogate ks/(2*gamma0^2) in field-free regions (BeamSolver.cpp:35-38; the caller
! supplies it), and beta0 is the reference beta from p0_mc. Combined with the
! per-particle theta_j = phi0 - ks*tau_j this reproduces Genesis's
! dtheta/ds = ks*(1 - 1/beta_z) + ku_like exactly, because
! -ks*dtau/ds = -ks*(1/beta_z - 1/beta0).
!
! The 1 - 1/beta0 factor is formed cancellation-free:
! 1 - 1/beta0 = -1/(beta0*(gamma0b^2)*(1+beta0)) with gamma0b^2 = p0_mc^2 + 1.
!-

function fel_phi0_rate (ks, ku_like, p0_mc) result (rate)

real(rp) ks, ku_like, p0_mc, rate
real(rp) gamma0b_sq, beta0

!

gamma0b_sq = p0_mc**2 + 1
beta0 = p0_mc / sqrt(gamma0b_sq)
rate = ku_like - ks / (beta0 * gamma0b_sq * (1 + beta0))

end function fel_phi0_rate

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_slice_diag (beam, sl, ks, diag)
!
! Routine to compute the per-slice beam diagnostics, weighted. Genesis's DiagBeam::calc
! (src/Core/Diagnostic.cpp:515-575) definitions with 1/N generalized to w_j/sum(w);
! uniform weights reproduce Genesis exactly up to the regrouped arithmetic. Positions and
! sizes are of x, y directly; pxposition is reported in Genesis's units (gamma*beta_x) for
! comparability. Adds N_eff and the derived current.
!-

subroutine fel_slice_diag (beam, sl, ks, diag)

type (fel_beam_struct) beam
type (fel_slice_struct) sl
type (fel_slice_diag_struct) diag
real(rp) ks
real(rp) g1, g2, x1, x2, y1, y2, px1, py1, br, bi, wsum, w2sum, w, gam, theta
integer ip

!

g1 = 0; g2 = 0; x1 = 0; x2 = 0; y1 = 0; y2 = 0; px1 = 0; py1 = 0
br = 0; bi = 0; wsum = 0; w2sum = 0

do ip = 1, sl%n
  w = sl%weight(ip)
  gam = fel_gamma_of(beam%p0_mc, sl%pz(ip))
  theta = fel_theta(beam, sl, ip, ks)

  wsum = wsum + w
  w2sum = w2sum + w*w
  x1 = x1 + w * sl%x(ip)
  x2 = x2 + w * sl%x(ip)**2
  y1 = y1 + w * sl%y(ip)
  y2 = y2 + w * sl%y(ip)**2
  g1 = g1 + w * gam
  g2 = g2 + w * gam**2
  px1 = px1 + w * sl%px(ip)
  py1 = py1 + w * sl%py(ip)
  br = br + w * cos(theta)
  bi = bi + w * sin(theta)
enddo

if (wsum <= 0) then
  diag = fel_slice_diag_struct()
  return
endif

g1 = g1/wsum;  g2 = g2/wsum
x1 = x1/wsum;  x2 = x2/wsum
y1 = y1/wsum;  y2 = y2/wsum

diag%energy = g1
diag%energyspread = sqrt(abs(g2 - g1*g1))
diag%xposition = x1
diag%xsize = sqrt(abs(x2 - x1*x1))
diag%yposition = y1
diag%ysize = sqrt(abs(y2 - y1*y1))
diag%pxposition = px1/wsum * beam%p0_mc
diag%pyposition = py1/wsum * beam%p0_mc
diag%bunching = sqrt((br/wsum)**2 + (bi/wsum)**2)
diag%bunchingphase = atan2(bi/wsum, br/wsum)
diag%n_eff = wsum*wsum / w2sum
diag%current = c_light * wsum / beam%slicelength

end subroutine fel_slice_diag

end module fel_beam_mod
