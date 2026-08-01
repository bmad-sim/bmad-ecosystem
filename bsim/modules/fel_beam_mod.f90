!+
! Module fel_beam_mod
!
! The packed particle representation for FEL tracking, its Genesis 1.3 Version 4 particle
! dump I/O, its conversion to and from Bmad's coord_struct at element boundaries, and the
! per-slice beam diagnostics.
!
! Why this exists (design brief section 4.2): the FEL step advances every particle of every
! slice every internal step, and coord_struct is roughly 224 bytes per particle against the
! 48 bytes actually needed, including a quad precision time. The FEL physics therefore
! operates on this packed structure-of-arrays form, and coord_struct appears only at element
! boundaries, converted by the routines here. That boundary is also exactly where a future
! GPU path would marshal, so the layout is chosen now (brief section 6.4): each slice has a
! capacity (the allocated size) and a fill count n, so that slice migration can later change
! populations without reallocating per step.
!
! Genesis coordinates, kept verbatim: x, y in meters; px = gamma*beta_x, py = gamma*beta_y
! (transverse momentum over m_e c, dimensionless); theta the ponderomotive phase in radians;
! gamma the Lorentz factor. theta wraps within one radiation wavelength while Bmad's z does
! not, so conversion of the longitudinal plane is NOT part of the coord_struct conversion
! here: the FEL program keeps theta in the packed array at all times and updates it across
! non-FEL elements from Bmad's change in z (see fel_theta_advance below).
!-

module fel_beam_mod

use bmad
use hdf5_interface

implicit none

! Genesis's own physical constants, transcribed from src/Main/GenMain.cpp:63-67 and used in
! every transcribed formula so that the arithmetic matches Genesis digit for digit. They
! deliberately differ from Bmad's constants: vacimp is the impedance of free space truncated
! to 376.73 (Bmad's mu_0_vac * c_light is 376.7303134...), and eev is the electron rest mass
! in eV. Do not "fix" these; the comparison against Genesis depends on them.

real(rp), parameter :: fel_vacimp = 376.73_rp        ! [Ohm]     GenMain.cpp:63
real(rp), parameter :: fel_eev    = 510998.95069_rp  ! [eV]      GenMain.cpp:67

!+
! Struct fel_slice_struct
!
! One beam slice in packed structure-of-arrays form. The arrays are allocated to a capacity
! that may exceed the fill count n; only elements 1:n are live. Genesis's Particle struct is
! the same six doubles in array-of-structures form.
!-

type fel_slice_struct
  real(rp), allocatable :: gamma(:)    ! Lorentz factor.
  real(rp), allocatable :: theta(:)    ! Ponderomotive phase [rad].
  real(rp), allocatable :: x(:)        ! [m]
  real(rp), allocatable :: y(:)        ! [m]
  real(rp), allocatable :: px(:)       ! gamma*beta_x, dimensionless.
  real(rp), allocatable :: py(:)       ! gamma*beta_y, dimensionless.
  integer :: n = 0                     ! Fill count. Live particles are 1:n.
  real(rp) :: current = 0              ! Slice current [A].
end type

!+
! Struct fel_beam_struct
!
! The whole beam: slices plus the global metadata the Genesis dump carries.
!-

type fel_beam_struct
  type (fel_slice_struct), allocatable :: slice(:)
  real(rp) :: reflength = 0        ! Reference (radiation) wavelength [m]. 'slicelength' in the dump.
  real(rp) :: slicelength = 0      ! Slice spacing [m]. 'slicespacing' in the dump.
  real(rp) :: s0 = 0               ! Start of the time window [m]. 'refposition' in the dump.
  integer :: nbins = 0             ! Beamlet size used at generation. Carried, not used here.
  logical :: one4one = .false.
end type

!+
! Struct fel_slice_diag_struct
!
! Per-slice beam diagnostics at one output position, matching Genesis's DiagBeam::calc
! definitions exactly: plain averages over particles, bunching b = |sum exp(i h theta)| / n.
!-

type fel_slice_diag_struct
  real(rp) :: energy = 0           ! <gamma>
  real(rp) :: energyspread = 0     ! sqrt(|<gamma^2> - <gamma>^2|)
  real(rp) :: bunching = 0         ! |b(1)|
  real(rp) :: bunchingphase = 0    ! arg b(1)
  real(rp) :: xposition = 0, yposition = 0    ! <x>, <y> [m]
  real(rp) :: xsize = 0, ysize = 0            ! rms sizes [m]
  real(rp) :: pxposition = 0, pyposition = 0  ! <px>, <py>
end type

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_read_genesis4_beam (beam, file_name, err_flag)
!
! Routine to read a Genesis 1.3 Version 4 particle dump (.par.h5) into a packed beam.
!
! The format, from src/IO/writeBeamHDF5.cpp: root scalars slicelength (the reference
! wavelength, note the name), slicespacing, refposition, beamletsize, slicecount, one4one;
! then one group per slice named slice000001 upward, each holding current (one value) and
! the six coordinate arrays gamma, theta, x, y, px, py.
!
! Input:
!   file_name   -- character(*): File to read.
!
! Output:
!   beam        -- fel_beam_struct: Beam read from the file.
!   err_flag    -- logical: Set True on error, False otherwise.
!-

subroutine fel_read_genesis4_beam (beam, file_name, err_flag)

type (fel_beam_struct), target :: beam
type (fel_slice_struct), pointer :: sl
integer(hid_t) f_id, g_id
integer is, n_slice, np, ivec(1), h5_err
real(rp) rvec(1)
logical err_flag, err
character(*) file_name
character(*), parameter :: r_name = 'fel_read_genesis4_beam'
character(20) group_name
type (hdf5_info_struct) info

!

err_flag = .true.

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

if (allocated(beam%slice)) deallocate(beam%slice)
allocate (beam%slice(n_slice))

do is = 1, n_slice
  sl => beam%slice(is)
  write (group_name, '(a, i0.6)') 'slice', is
  g_id = hdf5_open_group (f_id, trim(group_name), err, .true.);  if (err) return

  ! Take the particle count from the actual dataset extent, per the rule (FINDINGS.md 5.1)
  ! that H5LT reads have no buffer bound and extents are checked, never assumed.

  info = hdf5_object_info (g_id, 'gamma', err, .true.);  if (err) return
  np = int(info%data_dim(1))

  call fel_slice_reallocate (sl, np)
  sl%n = np

  call hdf5_read_dataset_real (g_id, 'current', rvec, err, trim(group_name) // '/current');  if (err) return
  sl%current = rvec(1)

  call hdf5_read_dataset_real (g_id, 'gamma', sl%gamma(1:np), err, trim(group_name) // '/gamma');  if (err) return
  call hdf5_read_dataset_real (g_id, 'theta', sl%theta(1:np), err, trim(group_name) // '/theta');  if (err) return
  call hdf5_read_dataset_real (g_id, 'x',     sl%x(1:np),     err, trim(group_name) // '/x');      if (err) return
  call hdf5_read_dataset_real (g_id, 'y',     sl%y(1:np),     err, trim(group_name) // '/y');      if (err) return
  call hdf5_read_dataset_real (g_id, 'px',    sl%px(1:np),    err, trim(group_name) // '/px');     if (err) return
  call hdf5_read_dataset_real (g_id, 'py',    sl%py(1:np),    err, trim(group_name) // '/py');     if (err) return

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
! Routine to write a packed beam as a Genesis 1.3 Version 4 particle dump, in the format
! documented at fel_read_genesis4_beam.
!
! Input:
!   beam        -- fel_beam_struct: Beam to write.
!   file_name   -- character(*): File to create.
!
! Output:
!   err_flag    -- logical: Set True on error, False otherwise.
!-

subroutine fel_write_genesis4_beam (beam, file_name, err_flag)

type (fel_beam_struct), target :: beam
type (fel_slice_struct), pointer :: sl
integer(hid_t) f_id, g_id
integer is, h5_err, one4one_int
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

  call hdf5_write_dataset_real (g_id, 'current', [sl%current],     err);  if (err) return
  call hdf5_write_dataset_real (g_id, 'gamma',   sl%gamma(1:sl%n), err);  if (err) return
  call hdf5_write_dataset_real (g_id, 'theta',   sl%theta(1:sl%n), err);  if (err) return
  call hdf5_write_dataset_real (g_id, 'x',       sl%x(1:sl%n),     err);  if (err) return
  call hdf5_write_dataset_real (g_id, 'y',       sl%y(1:sl%n),     err);  if (err) return
  call hdf5_write_dataset_real (g_id, 'px',      sl%px(1:sl%n),    err);  if (err) return
  call hdf5_write_dataset_real (g_id, 'py',      sl%py(1:sl%n),    err);  if (err) return

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

if (allocated(sl%gamma)) then
  if (size(sl%gamma) >= capacity) return
  call re_allocate (sl%gamma, capacity)
  call re_allocate (sl%theta, capacity)
  call re_allocate (sl%x, capacity)
  call re_allocate (sl%y, capacity)
  call re_allocate (sl%px, capacity)
  call re_allocate (sl%py, capacity)
else
  allocate (sl%gamma(capacity), sl%theta(capacity), sl%x(capacity), sl%y(capacity), &
            sl%px(capacity), sl%py(capacity))
endif

end subroutine fel_slice_reallocate

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_slice_to_bunch (sl, ele, bunch)
!
! Routine to convert a packed slice to a Bmad bunch_struct for tracking through a
! non-FEL element.
!
! The transverse map is exact: Bmad's px is the transverse momentum over the reference
! momentum, Genesis's is over m_e c, so the ratio is the reference gamma*beta taken from the
! element's p0c and Bmad's own electron mass. Using Bmad's mass here, not Genesis's eev,
! keeps the Bmad side self-consistent: p0c/mass_of(electron$) is exactly what Bmad's
! tracking divides by internally.
!
! The longitudinal plane: vec(5) (Bmad z) is set to zero for every particle. theta is not
! representable in a coord_struct and stays in the packed slice; after tracking, the change
! in vec(5) gives the change in theta through fel_theta_advance. vec(6) is set from gamma
! exactly.
!
! Input:
!   sl        -- fel_slice_struct: Packed slice.
!   ele       -- ele_struct: Element the bunch is about to be tracked through (start of).
!
! Output:
!   bunch     -- bunch_struct: Bunch ready for track1_bunch.
!-

subroutine fel_slice_to_bunch (sl, ele, bunch)

type (fel_slice_struct) sl
type (ele_struct) ele
type (bunch_struct) bunch
real(rp) p0_mc, mc2, p_mc, vec(6)
integer ip

!

mc2 = mass_of(electron$)
p0_mc = ele%value(p0c$) / mc2                  ! Reference gamma*beta.

if (allocated(bunch%particle)) then
  if (size(bunch%particle) /= sl%n) deallocate(bunch%particle)
endif
if (.not. allocated(bunch%particle)) allocate(bunch%particle(sl%n))

do ip = 1, sl%n
  p_mc = sqrt(sl%gamma(ip)**2 - 1)             ! Particle gamma*beta.

  vec(1) = sl%x(ip)
  vec(2) = sl%px(ip) / p0_mc
  vec(3) = sl%y(ip)
  vec(4) = sl%py(ip) / p0_mc
  vec(5) = 0
  vec(6) = (p_mc - p0_mc) / p0_mc

  ! init_coord derives beta, t and state from vec and the element consistently.
  ! shift_vec6 is explicitly false: it exists for elements whose reference momentum changes
  ! and must not touch vec(6) here.
  call init_coord (bunch%particle(ip), vec, ele, upstream_end$, electron$, shift_vec6 = .false.)
enddo

bunch%n_live = sl%n
bunch%charge_live = 0   ! Charge bookkeeping is not used in this deliverable.

end subroutine fel_slice_to_bunch

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_bunch_to_slice (bunch, ele, ks, gamma_ref, length, sl, err_flag)
!
! Routine to convert a tracked Bmad bunch back into the packed slice, advancing theta from
! the tracked change in Bmad's z. Inverse of fel_slice_to_bunch plus the theta update.
!
! The theta update, derived exactly rather than approximated. Genesis advances theta
! through a field-free region as (BeamSolver.cpp:35-37,161)
!
!   dtheta/dz = ks*(1 - 1/beta_z) + ks/(2*gamma_ref^2)
!
! Bmad's z advance through the same region is (track_a_drift.f90)
!
!   dz = L*( [beta/beta0 - 1] - [1/ps_rel - 1] )      with ps_rel = beta_z/beta
!
! so L/beta_z = L/beta0 - dz/beta, giving
!
!   dtheta = ks*L*(1 - 1/beta0 + 1/(2*gamma_ref^2)) + ks*dz/beta
!
! where beta0 is the reference beta from the element's p0c and beta is the particle's total
! beta, constant through a field-free region. The first term is particle independent; both
! terms are exact, no expansion in 1/gamma anywhere. Note the first term is a difference of
! numbers agreeing to ~4e-9, evaluated in doubles; Genesis's own evaluation of
! ks*(1-1/beta_z) has the same cancellation. This is an irreducible source of ~1e-6 rad
! level phase differences between the codes per interlude, documented in the test README.
!
! Input:
!   bunch       -- bunch_struct: Bunch after track1_bunch, with vec(5) started at 0.
!   ele         -- ele_struct: Element just tracked through (for p0c at its exit).
!   ks          -- real(rp): Radiation wavenumber twopi/lambda [1/m].
!   gamma_ref   -- real(rp): Genesis's reference gamma (gamma0 of the run).
!   length      -- real(rp): Element length [m].
!
! Output:
!   sl          -- fel_slice_struct: Updated packed slice.
!   err_flag    -- logical: Set True if a particle was lost in tracking, False otherwise.
!-

subroutine fel_bunch_to_slice (bunch, ele, ks, gamma_ref, length, sl, err_flag)

type (bunch_struct) bunch
type (ele_struct) ele
type (fel_slice_struct) sl
real(rp) ks, gamma_ref, length
real(rp) p0_mc, mc2, beta0, gamma0b, dtheta_ref, p_mc, beta
integer ip
logical err_flag
character(*), parameter :: r_name = 'fel_bunch_to_slice'

!

err_flag = .true.

mc2 = mass_of(electron$)
p0_mc = ele%value(p0c$) / mc2
gamma0b = sqrt(p0_mc**2 + 1)
beta0 = p0_mc / gamma0b

! The particle independent term, in cancellation-free form. Directly, 1 - 1/beta0 +
! 1/(2*gamma_ref^2) differences numbers agreeing to ~4e-9 and keeps ~1e-16 of absolute
! rounding, which times ks*L is microradians of phase noise. Exactly,
! 1 - 1/beta0 = -1/(beta0*gamma0b^2*(1+beta0)), so the whole coefficient is
! (1/2 - gamma_ref^2/(beta0*gamma0b^2*(1+beta0))) / gamma_ref^2, every factor well
! conditioned. gamma0b is the lattice reference gamma (from p0c); gamma_ref is Genesis's
! gamma0. They agree to ~1e-8 here but are kept distinct on principle.

dtheta_ref = ks * length * (0.5_rp - gamma_ref**2/(beta0*gamma0b**2*(1+beta0))) / gamma_ref**2

do ip = 1, sl%n
  if (bunch%particle(ip)%state /= alive$) then
    call out_io (s_error$, r_name, 'PARTICLE \i0\ LOST TRACKING THROUGH: ' // trim(ele%name), &
                                   i_array = [ip])
    return
  endif

  p_mc = p0_mc * (1 + bunch%particle(ip)%vec(6))
  sl%gamma(ip) = sqrt(p_mc**2 + 1)
  beta = p_mc / sl%gamma(ip)

  sl%x(ip)  = bunch%particle(ip)%vec(1)
  sl%px(ip) = bunch%particle(ip)%vec(2) * p0_mc
  sl%y(ip)  = bunch%particle(ip)%vec(3)
  sl%py(ip) = bunch%particle(ip)%vec(4) * p0_mc

  sl%theta(ip) = sl%theta(ip) + dtheta_ref + ks * bunch%particle(ip)%vec(5) / beta
enddo

err_flag = .false.

end subroutine fel_bunch_to_slice

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_slice_diag (sl, diag)
!
! Routine to compute the per-slice beam diagnostics, matching Genesis's DiagBeam::calc
! (src/Core/Diagnostic.cpp:515-575) definitions and accumulation order: plain sums over
! particles in storage order, normalized by n, bunching from sum of exp(i*theta).
!
! Input:
!   sl      -- fel_slice_struct: Slice.
!
! Output:
!   diag    -- fel_slice_diag_struct: Diagnostics.
!-

subroutine fel_slice_diag (sl, diag)

type (fel_slice_struct) sl
type (fel_slice_diag_struct) diag
real(rp) g1, g2, x1, x2, y1, y2, px1, py1, norm, br, bi
integer ip

!

g1 = 0; g2 = 0; x1 = 0; x2 = 0; y1 = 0; y2 = 0; px1 = 0; py1 = 0
br = 0; bi = 0

do ip = 1, sl%n
  x1 = x1 + sl%x(ip)
  x2 = x2 + sl%x(ip)**2
  y1 = y1 + sl%y(ip)
  y2 = y2 + sl%y(ip)**2
  g1 = g1 + sl%gamma(ip)
  g2 = g2 + sl%gamma(ip)**2
  px1 = px1 + sl%px(ip)
  py1 = py1 + sl%py(ip)
  br = br + cos(sl%theta(ip))
  bi = bi + sin(sl%theta(ip))
enddo

norm = 1
if (sl%n > 0) norm = 1.0_rp / sl%n

diag%energy = g1 * norm
diag%energyspread = sqrt(abs(g2*norm - (g1*norm)**2))
diag%xposition = x1 * norm
diag%xsize = sqrt(abs(x2*norm - (x1*norm)**2))
diag%yposition = y1 * norm
diag%ysize = sqrt(abs(y2*norm - (y1*norm)**2))
diag%pxposition = px1 * norm
diag%pyposition = py1 * norm
diag%bunching = sqrt((br*norm)**2 + (bi*norm)**2)
diag%bunchingphase = atan2(bi*norm, br*norm)

end subroutine fel_slice_diag

end module fel_beam_mod
