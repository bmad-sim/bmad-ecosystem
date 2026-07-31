!+
! Module wavefront_mod
!
! Representation and free-space propagation of a paraxial radiation wavefront.
!
! A wavefront_struct holds the complex transverse electric field components Ex and Ey
! sampled on a uniform 3D grid of shape (nx, ny, nz) in V/m. The first two indices are
! the transverse coordinates, the third is the longitudinal slice index. The wavefront
! propagates in the +z direction with the pulse head at max(z).
!
! This mirrors the Wavefront class of openPMD-beamphysics
! (beamphysics/wavefront/wavefront.py) closely enough that the two can be compared field
! by field. See bsim/wavefront/tests/run_validation.sh, which drifts the same input
! through both and reports the largest relative difference.
!
! Genesis 1.3 Version 4 HDF5 input and output for this structure is in
! wavefront_hdf5_mod.
!
! Index order. Fortran is column major, so Ex(:,:,iz) -- one longitudinal slice -- is
! contiguous. That is simultaneously what a per-slice 2D FFTW plan wants and what
! parallelising over slices wants. The Python class uses the same (nx, ny, nz) shape but
! in C order, so there its slices are strided.
!
! Grid convention, matching the Python class: the transverse grids are centered on zero,
! with x running from -(nx-1)*dx/2 to +(nx-1)*dx/2, and likewise for y and z. There is no
! separate origin offset. The one absolute position carried is ref_position, which locates
! the pulse along the lattice.
!-

module wavefront_mod

use, intrinsic :: iso_c_binding
use sim_utils

implicit none

! Kind of the complex field arrays.
!
! This is the single point of control for the field precision: every field array, and
! every routine argument carrying field values, is declared complex(wf_rp). Changing this
! parameter changes the representation everywhere without touching a call site.
!
! It is rp (double) now. Single precision is reserved for a future GPU path, and flipping
! this parameter is not by itself enough to get there: the plan cache below binds the
! double precision FFTW entry points, and its work buffer is declared complex(wf_rp), so a
! single precision field fails to compile at fftw_execute_dft with a type mismatch. That
! is deliberate. Transforming a single precision field in double precision would be a
! silent fallback, and a loud compile error is the correct way to be told that the fftwf_*
! entry points still have to be bound.

integer, parameter :: wf_rp = rp

! Transform directions for wavefront_fft2. Same sign convention as FFTW and as numpy.fft:
! forward carries exp(-i k x) and is unnormalized, backward carries exp(+i k x) and is
! unnormalized. Neither direction applies the 1/(nx*ny) factor; the caller does.
!
! Worth knowing when reading the validation results. Interchanging these two values does not
! change what wavefront_drift computes, and no test could be written that would notice. The
! backward transform is the forward transform with k negated, and the propagation kernel
! depends on kx^2 + ky^2 and so is even in k, so substituting k -> -k in the sum over k
! recovers the original expression term for term. The composition forward, multiply by the
! kernel, backward is therefore invariant under the interchange. The individual signs would
! matter for any operator whose kernel is not even in k.

integer, parameter :: wf_fft_forward$ = -1, wf_fft_backward$ = 1

!+
! Struct wavefront_struct
!
! Ex, Ey are allocated independently. An unallocated component means that polarisation
! component is absent, which is how the Python class's Ex = None is represented. At least
! one of the two must be allocated for the wavefront to be usable; wavefront_check tests
! that along with the rest of the Python class's __post_init__ validation.
!-

type wavefront_struct
  complex(wf_rp), allocatable :: Ex(:,:,:)   ! x polarized field [V/m], shape (nx, ny, nz).
  complex(wf_rp), allocatable :: Ey(:,:,:)   ! y polarized field [V/m], shape (nx, ny, nz).
  real(rp) :: dx = 1                         ! Transverse grid spacing in x [m].
  real(rp) :: dy = 1                         ! Transverse grid spacing in y [m].
  real(rp) :: dz = 1                         ! Slice spacing [m].
  real(rp) :: wavelength = 1                 ! Central wavelength [m].
  real(rp) :: ref_position = 0               ! Position of the pulse reference point along
                                             !   the lattice [m]. This is Genesis's
                                             !   refposition, whose job is keeping an
                                             !   imported field file and an imported beam
                                             !   file aligned in a start to end run. It is
                                             !   carried, not used, by this module.
end type

! Cached FFTW plan state for wavefront_fft2, private to the module.
!
! Following bmad/space_charge/fft_interface_mod.f90: plans are created once per transverse
! grid size and reused. Two differences from that routine. First, the work buffer here is
! allocated by fftw_alloc_complex and transforms are executed on it rather than on the
! caller's array, so the plan always sees the alignment it was created with; FFTW's
! new-array execute rule makes executing a plan on a differently aligned array undefined.
! Second, the buffer is module state, so wavefront_fft2 is not thread safe. Parallelising
! over slices, which is a later deliverable, wants one cache per thread; that is a change
! confined to these five variables and to wavefront_fft2.
!
! On the critical section inside wavefront_fft2: it guards FFTW's rule that plan creation is
! not reentrant, and nothing else. It does not make the routine callable from a parallel
! region, because wf_buf is written and read outside it. Both statements are repeated at the
! section itself, since that is where the pragma would otherwise be read as a thread safety
! claim. Note the per-thread version still needs the section, because the FFTW planner is
! globally serialised regardless of how many buffers there are.

! The transform itself is deliberately single threaded. The parallelisation axis for a
! wavefront is the slice index, not the transverse transform, so fftw_plan_with_nthreads is
! not called here: threading a small 2D transform underneath a parallel loop over slices
! would only oversubscribe.

type(C_PTR), private, save :: wf_plan_fwd = C_NULL_PTR
type(C_PTR), private, save :: wf_plan_bwd = C_NULL_PTR
type(C_PTR), private, save :: wf_buf_cptr = C_NULL_PTR
complex(wf_rp), private, pointer, save :: wf_buf(:,:) => null()
integer, private, save :: wf_cache_nx = 0, wf_cache_ny = 0

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine wavefront_init (wf, nx, ny, nz, dx, dy, dz, wavelength, polarization, ref_position)
!
! Routine to allocate the field arrays of a wavefront and set its grid parameters.
! The allocated components are zeroed.
!
! Input:
!   nx, ny, nz    -- integer: Grid dimensions.
!   dx, dy, dz    -- real(rp): Grid spacings [m]. Must be positive.
!   wavelength    -- real(rp): Central wavelength [m]. Must be positive.
!   polarization  -- character(*), optional: Which components to allocate. One of 'x',
!                      'y' or 'xy'. Default is 'x'.
!   ref_position  -- real(rp), optional: Pulse reference position along the lattice [m].
!                      Default is 0.
!
! Output:
!   wf            -- wavefront_struct: Initialized wavefront.
!-

subroutine wavefront_init (wf, nx, ny, nz, dx, dy, dz, wavelength, polarization, ref_position)

type (wavefront_struct) wf
integer nx, ny, nz
real(rp) dx, dy, dz, wavelength
real(rp), optional :: ref_position
character(*), optional :: polarization
character(*), parameter :: r_name = 'wavefront_init'
character(4) pol

!

pol = 'x'
if (present(polarization)) pol = polarization

if (pol /= 'x' .and. pol /= 'y' .and. pol /= 'xy') then
  call out_io (s_error$, r_name, 'POLARIZATION MUST BE "x", "y" OR "xy". GOT: ' // pol)
  return
endif

if (allocated(wf%Ex)) deallocate(wf%Ex)
if (allocated(wf%Ey)) deallocate(wf%Ey)

if (index(pol, 'x') /= 0) then
  allocate (wf%Ex(nx, ny, nz))
  wf%Ex = 0
endif

if (index(pol, 'y') /= 0) then
  allocate (wf%Ey(nx, ny, nz))
  wf%Ey = 0
endif

wf%dx = dx
wf%dy = dy
wf%dz = dz
wf%wavelength = wavelength
wf%ref_position = real_option(0.0_rp, ref_position)

end subroutine wavefront_init

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine wavefront_check (wf, err_flag)
!
! Routine to validate a wavefront. This is the Fortran counterpart of the Python class's
! __post_init__: at least one field component must be present, the two components must
! have the same shape if both are present, and dx, dy, dz and wavelength must be positive.
!
! Input:
!   wf        -- wavefront_struct: Wavefront to check.
!
! Output:
!   err_flag  -- logical: Set True if the wavefront is not valid, False otherwise.
!-

subroutine wavefront_check (wf, err_flag)

type (wavefront_struct) wf
logical err_flag
character(*), parameter :: r_name = 'wavefront_check'

!

err_flag = .true.

if (.not. allocated(wf%Ex) .and. .not. allocated(wf%Ey)) then
  call out_io (s_error$, r_name, 'AT LEAST ONE OF Ex OR Ey MUST BE PRESENT IN A WAVEFRONT.')
  return
endif

if (allocated(wf%Ex) .and. allocated(wf%Ey)) then
  if (any(shape(wf%Ex) /= shape(wf%Ey))) then
    call out_io (s_error$, r_name, 'Ex SHAPE \3i6\ DOES NOT MATCH Ey SHAPE \3i6\.', &
                                                 i_array = [shape(wf%Ex), shape(wf%Ey)])
    return
  endif
endif

if (wf%dx <= 0 .or. wf%dy <= 0 .or. wf%dz <= 0 .or. wf%wavelength <= 0) then
  call out_io (s_error$, r_name, 'dx, dy, dz AND wavelength MUST ALL BE POSITIVE. GOT: \4es14.6\', &
                                                 r_array = [wf%dx, wf%dy, wf%dz, wf%wavelength])
  return
endif

err_flag = .false.

end subroutine wavefront_check

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function wavefront_shape (wf) result (n_grid)
!
! Routine to return the shape of the field arrays of a wavefront.
!
! Input:
!   wf          -- wavefront_struct: Wavefront.
!
! Output:
!   n_grid(3)   -- integer: [nx, ny, nz]. Set to zero if neither component is present.
!-

function wavefront_shape (wf) result (n_grid)

type (wavefront_struct) wf
integer n_grid(3)

!

if (allocated(wf%Ex)) then
  n_grid = shape(wf%Ex)
elseif (allocated(wf%Ey)) then
  n_grid = shape(wf%Ey)
else
  n_grid = 0
endif

end function wavefront_shape

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function wavefront_k0 (wf) result (k0)
!
! Routine to return the central wavenumber twopi / wavelength [rad/m].
!-

function wavefront_k0 (wf) result (k0)

type (wavefront_struct) wf
real(rp) k0

!

k0 = twopi / wf%wavelength

end function wavefront_k0

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function wavefront_photon_energy (wf) result (e_photon)
!
! Routine to return the central photon energy h_bar * c * k0 [eV].
!-

function wavefront_photon_energy (wf) result (e_photon)

type (wavefront_struct) wf
real(rp) e_photon

!

e_photon = wavefront_k0(wf) * h_bar_planck * c_light

end function wavefront_photon_energy

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function wavefront_coord_vec (n_pt, d_pt) result (vec)
!
! Routine to return the coordinate vector of a zero centered uniform grid:
! n_pt points spaced d_pt apart running from -(n_pt-1)*d_pt/2 to +(n_pt-1)*d_pt/2.
!
! wavefront_xvec, wavefront_yvec and wavefront_zvec are this applied to the three axes.
!
! Input:
!   n_pt      -- integer: Number of points.
!   d_pt      -- real(rp): Point spacing.
!
! Output:
!   vec(n_pt) -- real(rp), allocatable: Coordinates.
!-

function wavefront_coord_vec (n_pt, d_pt) result (vec)

integer n_pt, i
real(rp) d_pt
real(rp), allocatable :: vec(:)

!

allocate (vec(n_pt))
do i = 1, n_pt
  vec(i) = (i - 1 - (n_pt - 1) / 2.0_rp) * d_pt
enddo

end function wavefront_coord_vec

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function wavefront_xvec (wf) result (vec)
!
! Routine to return the x coordinates of the grid points [m].
!-

function wavefront_xvec (wf) result (vec)

type (wavefront_struct) wf
real(rp), allocatable :: vec(:)
integer n_grid(3)

!

n_grid = wavefront_shape(wf)
vec = wavefront_coord_vec(n_grid(1), wf%dx)

end function wavefront_xvec

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function wavefront_yvec (wf) result (vec)
!
! Routine to return the y coordinates of the grid points [m].
!-

function wavefront_yvec (wf) result (vec)

type (wavefront_struct) wf
real(rp), allocatable :: vec(:)
integer n_grid(3)

!

n_grid = wavefront_shape(wf)
vec = wavefront_coord_vec(n_grid(2), wf%dy)

end function wavefront_yvec

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function wavefront_zvec (wf) result (vec)
!
! Routine to return the z coordinates of the slices [m].
!-

function wavefront_zvec (wf) result (vec)

type (wavefront_struct) wf
real(rp), allocatable :: vec(:)
integer n_grid(3)

!

n_grid = wavefront_shape(wf)
vec = wavefront_coord_vec(n_grid(3), wf%dz)

end function wavefront_zvec

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function wavefront_fft_wavenumber_vec (n_pt, d_pt) result (k_vec)
!
! Routine to return the transverse wavenumbers in FFT storage order, that is, in the order
! the values come out of a forward transform: zero frequency first, ascending positive
! frequencies, then the negative frequencies ascending to just below zero. This is
! twopi * numpy.fft.fftfreq(n_pt, d = d_pt), and no fftshift is applied or wanted since a
! propagation kernel built from these is used between a forward and a backward transform.
!
! Input:
!   n_pt        -- integer: Number of points.
!   d_pt        -- real(rp): Point spacing [m].
!
! Output:
!   k_vec(n_pt) -- real(rp), allocatable: Wavenumbers [rad/m].
!-

function wavefront_fft_wavenumber_vec (n_pt, d_pt) result (k_vec)

integer n_pt, i, m
real(rp) d_pt
real(rp), allocatable :: k_vec(:)

!

allocate (k_vec(n_pt))

do i = 1, n_pt
  m = i - 1
  if (2 * m > n_pt - 1) m = m - n_pt   ! The second half of the array holds negative wavenumbers.
  k_vec(i) = twopi * m / (n_pt * d_pt)
enddo

end function wavefront_fft_wavenumber_vec

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function wavefront_intensity (wf) result (intensity)
!
! Routine to return the field intensity c * eps_0 / 2 * (|Ex|^2 + |Ey|^2) [W/m^2].
!
! Output:
!   intensity(nx,ny,nz)  -- real(rp), allocatable: Intensity.
!-

function wavefront_intensity (wf) result (intensity)

type (wavefront_struct) wf
real(rp), allocatable :: intensity(:,:,:)
integer n_grid(3)

!

n_grid = wavefront_shape(wf)
allocate (intensity(n_grid(1), n_grid(2), n_grid(3)))
intensity = 0

if (allocated(wf%Ex)) intensity = intensity + abs(wf%Ex)**2
if (allocated(wf%Ey)) intensity = intensity + abs(wf%Ey)**2

intensity = c_light * eps_0_vac / 2 * intensity

end function wavefront_intensity

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function wavefront_energy (wf) result (energy)
!
! Routine to return the total field energy eps_0/2 * Integral |E|^2 dx dy dz [J].
!
! Note this is conserved exactly by wavefront_drift, since the propagation kernel has unit
! modulus everywhere and the transform pair is unitary up to the 1/(nx*ny) factor. It is
! therefore a useful invariant to watch.
!-

function wavefront_energy (wf) result (energy)

type (wavefront_struct) wf
real(rp) energy

!

energy = sum(wavefront_intensity(wf)) / c_light * wf%dx * wf%dy * wf%dz

end function wavefront_energy

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function wavefront_fluence (wf) result (fluence)
!
! Routine to return the transverse fluence eps_0/2 * Integral |E(x,y,z)|^2 dz [J/m^2].
!
! Output:
!   fluence(nx,ny) -- real(rp), allocatable: Fluence.
!-

function wavefront_fluence (wf) result (fluence)

type (wavefront_struct) wf
real(rp), allocatable :: fluence(:,:)

!

fluence = sum(wavefront_intensity(wf), dim = 3) / c_light * wf%dz

end function wavefront_fluence

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function wavefront_power (wf) result (power)
!
! Routine to return the longitudinal power profile Integral I(x,y,z) dx dy [W].
!
! Output:
!   power(nz)  -- real(rp), allocatable: Power per slice.
!-

function wavefront_power (wf) result (power)

type (wavefront_struct) wf
real(rp), allocatable :: power(:)

!

power = sum(sum(wavefront_intensity(wf), dim = 1), dim = 1) * wf%dx * wf%dy

end function wavefront_power

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine wavefront_transverse_moments (wf, mean_x, mean_y, sigma_x, sigma_y)
!
! Routine to return the intensity weighted transverse centroid and rms size of a
! wavefront, summed over all slices. Mirrors the Python class's mean_x, mean_y, sigma_x
! and sigma_y.
!
! These give an independent analytic check on wavefront_drift: a Gaussian of waist size
! sigma0 at the waist has sigma(z) = sqrt(sigma0^2 + (z * wavelength / (fourpi * sigma0))^2)
! after drifting a distance z.
!
! Input:
!   wf        -- wavefront_struct: Wavefront.
!
! Output:
!   mean_x    -- real(rp): Intensity weighted <x> [m].
!   mean_y    -- real(rp): Intensity weighted <y> [m].
!   sigma_x   -- real(rp): Intensity weighted rms x [m].
!   sigma_y   -- real(rp): Intensity weighted rms y [m].
!-

subroutine wavefront_transverse_moments (wf, mean_x, mean_y, sigma_x, sigma_y)

type (wavefront_struct) wf
real(rp) mean_x, mean_y, sigma_x, sigma_y
real(rp) wt_sum
real(rp), allocatable :: intensity(:,:,:), wt_x(:), wt_y(:), xvec(:), yvec(:)

!

intensity = wavefront_intensity(wf)
wt_x = sum(sum(intensity, dim = 3), dim = 2)
wt_y = sum(sum(intensity, dim = 3), dim = 1)
xvec = wavefront_xvec(wf)
yvec = wavefront_yvec(wf)

wt_sum = sum(wt_x)
if (wt_sum <= 0) then
  mean_x = 0; mean_y = 0; sigma_x = 0; sigma_y = 0
  return
endif

mean_x = sum(wt_x * xvec) / wt_sum
mean_y = sum(wt_y * yvec) / wt_sum

! Two pass variance rather than <x^2> - <x>^2. The one pass form differences two numbers of
! similar size whenever the centroid is large against the width, which is the normal
! situation for an off-axis beam, and loses most of its significant figures doing it. This
! is also the form mean_variance_calc uses in openPMD-beamphysics.

sigma_x = sqrt(sum(wt_x * (xvec - mean_x)**2) / wt_sum)
sigma_y = sqrt(sum(wt_y * (yvec - mean_y)**2) / wt_sum)

end subroutine wavefront_transverse_moments

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine wavefront_drift (wf, z_drift, err_flag, curvature)
!
! Routine to propagate a wavefront a distance z_drift through free space.
!
! Each slice is propagated independently with the exact transfer function of the paraxial
! wave equation, applied in transverse k-space:
!
!   E(kx, ky) -> E(kx, ky) * exp(-i * z_drift * wavelength * (kx^2 + ky^2) / fourpi)
!
! The exponent is -i * z_drift * (kx^2 + ky^2) / (2 * k0) written in terms of the
! wavelength, wavelength / fourpi being 1 / (2 * k0). This is spectrally exact and
! unconditionally stable: it is not a finite difference step and there is no accuracy to
! be gained by subdividing z_drift. It is the same kernel as Genesis's FieldSolverFFT and
! as drift_wavefront_basic in openPMD-beamphysics.
!
! The transverse boundary is periodic, since that is what a discrete Fourier transform
! imposes and no apodisation is applied. Field reaching the transverse edge of the grid
! wraps around and re-enters. Over a long drift the grid must therefore be large enough to
! contain the field, or padded until it is.
!
! Input:
!   wf          -- wavefront_struct: Wavefront to propagate.
!   z_drift     -- real(rp): Drift distance [m]. May be negative.
!   curvature   -- real(rp), optional: Accepted only as zero. A nonzero curvature selects
!                    quadratic phase rescaling with an expanding grid, which is not
!                    implemented here; passing one is an error rather than being silently
!                    ignored.
!
! Output:
!   wf          -- wavefront_struct: Propagated wavefront.
!   err_flag    -- logical, optional: Set True on error, False otherwise.
!-

subroutine wavefront_drift (wf, z_drift, err_flag, curvature)

type (wavefront_struct), target :: wf
real(rp) z_drift
real(rp), optional :: curvature
logical, optional :: err_flag
logical err

integer n_grid(3), nx, ny, nz, iz, i_pol
complex(wf_rp), allocatable :: kernel(:,:)
complex(wf_rp), pointer :: field(:,:,:)
character(*), parameter :: r_name = 'wavefront_drift'

!

if (present(err_flag)) err_flag = .true.

call wavefront_check (wf, err);  if (err) return

if (present(curvature)) then
  if (curvature /= 0) then
    call out_io (s_error$, r_name, 'NONZERO CURVATURE IS NOT IMPLEMENTED BY THIS PROPAGATOR.', &
                 'Curvature corrected propagation with grid rescaling, as in', &
                 'drift_wavefront_advanced of openPMD-beamphysics, is a separate propagator.')
    return
  endif
endif

! Note err_flag stays true from here until the transforms have actually finished. Clearing it
! early and then returning out of the loop below on a failed transform would report success on
! failure.

if (z_drift == 0) then
  if (present(err_flag)) err_flag = .false.
  return
endif

n_grid = wavefront_shape(wf)
nx = n_grid(1); ny = n_grid(2); nz = n_grid(3)

kernel = wavefront_drift_kernel(wf, z_drift)

! Apply slice by slice. Both polarisation components see the same kernel: the paraxial
! operator does not couple Ex and Ey in free space.

do i_pol = 1, 2
  if (i_pol == 1) then
    if (.not. allocated(wf%Ex)) cycle
    field => wf%Ex
  else
    if (.not. allocated(wf%Ey)) cycle
    field => wf%Ey
  endif

  do iz = 1, nz
    call wavefront_fft2 (field(:,:,iz), wf_fft_forward$, err);  if (err) return
    field(:,:,iz) = kernel * field(:,:,iz)
    call wavefront_fft2 (field(:,:,iz), wf_fft_backward$, err);  if (err) return
    field(:,:,iz) = field(:,:,iz) / real(nx * ny, wf_rp)
  enddo
enddo

if (present(err_flag)) err_flag = .false.

end subroutine wavefront_drift

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function wavefront_drift_kernel (wf, z_drift) result (kernel)
!
! Routine to return the transverse k-space propagation kernel for a drift of z_drift:
!
!   exp(-i * z_drift * wavelength * (kx^2 + ky^2) / fourpi)
!
! in FFT storage order, so that it can be applied between a forward and a backward
! transform with no fftshift.
!
! The arithmetic here is deliberately written in the same association order as
! drift_wavefront_basic in openPMD-beamphysics, so that the two implementations differ by
! round-off rather than by grouping.
!
! Input:
!   wf              -- wavefront_struct: Supplies nx, ny, dx, dy and wavelength.
!   z_drift         -- real(rp): Drift distance [m].
!
! Output:
!   kernel(nx,ny)   -- complex(wf_rp), allocatable: Propagation kernel.
!-

function wavefront_drift_kernel (wf, z_drift) result (kernel)

type (wavefront_struct) wf
real(rp) z_drift
complex(wf_rp), allocatable :: kernel(:,:)
integer n_grid(3), nx, ny, ix, iy
real(rp), allocatable :: kx_vec(:), ky_vec(:)
real(rp) phase, k2_scale

!

n_grid = wavefront_shape(wf)
nx = n_grid(1); ny = n_grid(2)

kx_vec = wavefront_fft_wavenumber_vec(nx, wf%dx)
ky_vec = wavefront_fft_wavenumber_vec(ny, wf%dy)
k2_scale = -z_drift * wf%wavelength

allocate (kernel(nx, ny))
do iy = 1, ny
  do ix = 1, nx
    phase = k2_scale * (kx_vec(ix)**2 + ky_vec(iy)**2) / (4 * pi)
    kernel(ix, iy) = cmplx(cos(phase), sin(phase), wf_rp)
  enddo
enddo

end function wavefront_drift_kernel

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine wavefront_drift_reference (wf, z_drift, err_flag)
!
! Reference implementation of wavefront_drift that shares no code with it.
!
! This exists as a validation instrument, not as a propagator anyone should use: it costs
! O(nx*ny*(nx+ny)) per slice against the FFT's O(nx*ny*log(nx*ny)). It builds its own
! propagation kernel and performs its own transform, so it is independent both of
! wavefront_drift_kernel and of every FFTW convention -- the dimension order handed to
! fftw_plan_dft_2d, the sign attached to each direction, the normalization factor, and the
! alignment rules for executing a cached plan on a differently allocated array.
!
! Why that independence is worth the code. The comparison against openPMD-beamphysics goes
! through the Genesis4 file format, which requires nx = ny and dx = dy. On a square grid
! with equal spacings the propagation kernel is symmetric under interchanging the two
! transverse axes, so that comparison cannot see a transposed transform, nor an
! interchanged dx and dy, however asymmetric the test field is. Both were confirmed by
! mutation to pass the Python comparison at round-off. Only a rectangular grid with unequal
! spacings can see them, and only this routine can supply a reference on one.
!
! The kernel below is deliberately written in a different but equivalent form from
! wavefront_drift_kernel: the wavenumbers are formed inline rather than through
! wavefront_fft_wavenumber_vec, and the coefficient is 1/(2*k0) rather than
! wavelength/fourpi. That is what makes the comparison a check on the arithmetic and not
! merely on the transform. The two therefore differ by round-off in the phase as well as in
! the transform.
!
! Being a direct sum, this routine also accumulates round-off as O(n) rather than the FFT's
! O(log n), so it is the less accurate of the two. Expect agreement near 1e-14 rather than
! 1e-16, and do not read the residual as an error in wavefront_drift.
!
! Input:
!   wf          -- wavefront_struct: Wavefront to propagate.
!   z_drift     -- real(rp): Drift distance [m].
!
! Output:
!   wf          -- wavefront_struct: Propagated wavefront.
!   err_flag    -- logical, optional: Set True on error, False otherwise.
!-

subroutine wavefront_drift_reference (wf, z_drift, err_flag)

type (wavefront_struct), target :: wf
real(rp) z_drift
logical, optional :: err_flag
logical err

integer n_grid(3), nx, ny, nz, ix, iy, iz, i_pol, mx, my
complex(wf_rp), allocatable :: kernel(:,:), slice(:,:)
complex(wf_rp), pointer :: field(:,:,:)
real(rp) k0, kx, ky, phase

!

if (present(err_flag)) err_flag = .true.
call wavefront_check (wf, err);  if (err) return

! As in wavefront_drift, err_flag is cleared only once the work is done. Nothing below can
! fail today, since the direct sum has no fallible call in it, but clearing it up here is the
! shape that reports success on failure the moment one is added.

if (z_drift == 0) then
  if (present(err_flag)) err_flag = .false.
  return
endif

n_grid = wavefront_shape(wf)
nx = n_grid(1); ny = n_grid(2); nz = n_grid(3)

! Kernel, built here rather than shared. See the note above on why.

k0 = twopi / wf%wavelength
allocate (kernel(nx, ny), slice(nx, ny))

do iy = 1, ny
  my = iy - 1
  if (my > (ny - 1) / 2) my = my - ny
  ky = my * twopi / (ny * wf%dy)

  do ix = 1, nx
    mx = ix - 1
    if (mx > (nx - 1) / 2) mx = mx - nx
    kx = mx * twopi / (nx * wf%dx)

    phase = -z_drift * (kx * kx + ky * ky) / (2 * k0)
    kernel(ix, iy) = cmplx(cos(phase), sin(phase), wf_rp)
  enddo
enddo

do i_pol = 1, 2
  if (i_pol == 1) then
    if (.not. allocated(wf%Ex)) cycle
    field => wf%Ex
  else
    if (.not. allocated(wf%Ey)) cycle
    field => wf%Ey
  endif

  do iz = 1, nz
    ! Forward: transform along x, then along y.
    do iy = 1, ny
      slice(:,iy) = wavefront_dft_1d(field(:,iy,iz), wf_fft_forward$)
    enddo
    do ix = 1, nx
      slice(ix,:) = wavefront_dft_1d(slice(ix,:), wf_fft_forward$)
    enddo

    slice = kernel * slice

    ! Backward, same two passes with the opposite sign, then the 1/(nx*ny) normalization.
    do iy = 1, ny
      slice(:,iy) = wavefront_dft_1d(slice(:,iy), wf_fft_backward$)
    enddo
    do ix = 1, nx
      slice(ix,:) = wavefront_dft_1d(slice(ix,:), wf_fft_backward$)
    enddo

    field(:,:,iz) = slice / real(nx * ny, wf_rp)
  enddo
enddo

if (present(err_flag)) err_flag = .false.

end subroutine wavefront_drift_reference

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function wavefront_dft_1d (dat, direction) result (out)
!
! Routine to return the unnormalized discrete Fourier transform of a vector, computed as a
! direct sum:
!
!   out(k) = Sum_j dat(j) * exp(direction * i * twopi * (j-1) * (k-1) / n)
!
! with direction = wf_fft_forward$ giving the negative exponent, matching FFTW and
! numpy.fft. Used only by wavefront_drift_reference; see the discussion there.
!
! Input:
!   dat(:)      -- complex(wf_rp): Data to transform.
!   direction   -- integer: wf_fft_forward$ or wf_fft_backward$.
!
! Output:
!   out(size(dat)) -- complex(wf_rp), allocatable: Transform.
!-

function wavefront_dft_1d (dat, direction) result (out)

complex(wf_rp) dat(:)
complex(wf_rp), allocatable :: out(:)
integer direction, n_pt, j, k
integer(8) n8
real(rp) angle

!

n_pt = size(dat)
n8 = n_pt
allocate (out(n_pt))

! The index product is reduced modulo n before being turned into an angle. That keeps the
! angle inside one turn, which both avoids a needless argument reduction in cos and sin and
! keeps the product from overflowing a default integer at large n.

do k = 1, n_pt
  out(k) = 0
  do j = 1, n_pt
    angle = direction * twopi * modulo(int(j-1, 8) * int(k-1, 8), n8) / n_pt
    out(k) = out(k) + dat(j) * cmplx(cos(angle), sin(angle), wf_rp)
  enddo
enddo

end function wavefront_dft_1d

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine wavefront_fft2 (dat, direction, err_flag)
!
! Routine to apply an unnormalized in-place 2D complex FFT to dat using a cached FFTW plan.
!
! Neither direction applies a normalization factor, so a forward followed by a backward
! transform multiplies the data by nx*ny. This matches FFTW and matches numpy.fft.fft2;
! numpy's ifft2 differs from wf_fft_backward$ only by that factor.
!
! Plans are created once per grid size and cached, following
! bmad/space_charge/fft_interface_mod.f90. The transform is executed on an internal
! fftw_alloc_complex buffer rather than on dat, because a plan may only be executed on an
! array with the alignment it was created with. dat is copied in and out.
!
! Not thread safe: the plan cache and its work buffer are module state. See the note at
! the cache declaration.
!
! Input:
!   dat(:,:)    -- complex(wf_rp): Data to transform.
!   direction   -- integer: wf_fft_forward$ or wf_fft_backward$.
!
! Output:
!   dat(:,:)    -- complex(wf_rp): Transformed data.
!   err_flag    -- logical, optional: Set True on error, False otherwise.
!-

subroutine wavefront_fft2 (dat, direction, err_flag)

include 'fftw3.f03'

complex(wf_rp) dat(:,:)
integer direction
logical, optional :: err_flag
logical cache_ok
integer nx, ny
character(*), parameter :: r_name = 'wavefront_fft2'

!

if (present(err_flag)) err_flag = .true.

if (direction /= wf_fft_forward$ .and. direction /= wf_fft_backward$) then
  call out_io (s_error$, r_name, 'BAD TRANSFORM DIRECTION: \i0\ ', i_array = [direction])
  return
endif

nx = size(dat, 1)
ny = size(dat, 2)

! The critical section below guards one specific thing: FFTW's rule that plan creation is not
! reentrant. It does NOT make this routine thread safe, and must not be read as doing so. The
! work buffer wf_buf is module state, written and read outside the section, so two threads in
! here at once corrupt each other's data whatever the planner does. See the note at the cache
! declaration for what parallelising over slices actually requires.
!
! Failures are recorded in wf_cache_nx and reported after the section rather than returned
! from inside it, since branching out of a critical construct is not allowed.

!$OMP CRITICAL (wavefront_fft_plan_lock)

if (nx /= wf_cache_nx .or. ny /= wf_cache_ny) then
  if (C_ASSOCIATED(wf_plan_fwd)) call fftw_destroy_plan(wf_plan_fwd)
  if (C_ASSOCIATED(wf_plan_bwd)) call fftw_destroy_plan(wf_plan_bwd)
  if (C_ASSOCIATED(wf_buf_cptr)) call fftw_free(wf_buf_cptr)
  wf_plan_fwd = C_NULL_PTR
  wf_plan_bwd = C_NULL_PTR
  wf_buf_cptr = C_NULL_PTR
  wf_buf => null()

  ! Mark the cache invalid up front, so that any failure below leaves it invalid rather than
  ! half built. It is set to the real size only once every piece has been obtained.
  wf_cache_nx = 0; wf_cache_ny = 0

  ! fftw_alloc_complex and fftw_plan_dft_2d both return null on failure, and a null plan
  ! reaching fftw_execute_dft is a crash rather than a diagnosable error.
  wf_buf_cptr = fftw_alloc_complex(int(nx, C_SIZE_T) * int(ny, C_SIZE_T))

  if (C_ASSOCIATED(wf_buf_cptr)) then
    call C_F_POINTER (wf_buf_cptr, wf_buf, [nx, ny])

    ! Note the reversed dimension order: fftw_plan_dft_2d takes the slowest varying
    ! dimension first, and in a Fortran (nx, ny) array that is ny. Getting this backwards is
    ! invisible on a square grid, which is why wavefront_drift_reference exists.
    wf_plan_fwd = fftw_plan_dft_2d(ny, nx, wf_buf, wf_buf, FFTW_FORWARD,  FFTW_MEASURE)
    wf_plan_bwd = fftw_plan_dft_2d(ny, nx, wf_buf, wf_buf, FFTW_BACKWARD, FFTW_MEASURE)

    if (C_ASSOCIATED(wf_plan_fwd) .and. C_ASSOCIATED(wf_plan_bwd)) then
      wf_cache_nx = nx; wf_cache_ny = ny
    endif
  endif
endif

cache_ok = (wf_cache_nx == nx .and. wf_cache_ny == ny)

!$OMP END CRITICAL (wavefront_fft_plan_lock)

if (.not. cache_ok) then
  call out_io (s_error$, r_name, &
        'FFTW COULD NOT ALLOCATE A WORK BUFFER OR CREATE A PLAN FOR A \i0\ BY \i0\ TRANSFORM.', &
        i_array = [nx, ny])
  return
endif

wf_buf = dat

! Passing the same array as both in and out is FFTW's documented way to ask for an in-place
! transform, and is what bmad/space_charge/fft_interface_mod.f90 does. gfortran with -Wall
! warns about it, because the fftw3.f03 interface declares both dummies intent(out) and
! aliasing those would be nonconforming for a normal Fortran procedure. It is bind(C) with
! assumed-size dummies, so the arguments are passed by address with no copy in or copy out,
! and there is nothing to go wrong. Bmad's own build does not enable that warning.

if (direction == wf_fft_forward$) then
  call fftw_execute_dft (wf_plan_fwd, wf_buf, wf_buf)
else
  call fftw_execute_dft (wf_plan_bwd, wf_buf, wf_buf)
endif

dat = wf_buf

if (present(err_flag)) err_flag = .false.

end subroutine wavefront_fft2

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine wavefront_fft_free ()
!
! Routine to destroy the cached FFTW plans and free the work buffer. Not needed for
! correctness; useful for making a leak check clean.
!-

subroutine wavefront_fft_free ()

include 'fftw3.f03'

!

!$OMP CRITICAL (wavefront_fft_plan_lock)
if (C_ASSOCIATED(wf_plan_fwd)) call fftw_destroy_plan(wf_plan_fwd)
if (C_ASSOCIATED(wf_plan_bwd)) call fftw_destroy_plan(wf_plan_bwd)
if (C_ASSOCIATED(wf_buf_cptr)) call fftw_free(wf_buf_cptr)
wf_plan_fwd = C_NULL_PTR
wf_plan_bwd = C_NULL_PTR
wf_buf_cptr = C_NULL_PTR
wf_buf => null()
wf_cache_nx = 0; wf_cache_ny = 0
!$OMP END CRITICAL (wavefront_fft_plan_lock)

end subroutine wavefront_fft_free

end module wavefront_mod
