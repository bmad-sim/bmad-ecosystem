!+
! Program wavefront_drift_test
!
! Test driver for the Fortran free-space wavefront propagator. Two modes.
!
!   wavefront_drift_test <input_file> <output_file> <z_drift>
!
! Reads a Genesis 1.3 Version 4 field file into a wavefront_struct, drifts it, and writes
! the result back out in the same format. This is the Fortran half of the Fortran against
! Python comparison; tests/run_validation.sh drives both halves and reports the largest
! relative difference.
!
! The input carries a single polarisation component, which the reader puts in Ex. This mode
! copies it to Ey as well, so that both components are propagated, and then checks that they
! came out identical. That is a cheap test of the claim that free-space propagation does not
! couple the two: the components go through separate arrays and separate transforms, so a
! mistake that mixed them up, or treated one differently, would show. Only Ex is written to
! the output file, since the Genesis4 format holds one component.
!
!   wavefront_drift_test -self_check [tolerance]
!
! Compares wavefront_drift against wavefront_drift_reference, the longhand DFT
! implementation, on rectangular grids with unequal transverse spacings.
!
! This mode exists because the comparison against openPMD-beamphysics has a blind spot it
! cannot close. That comparison passes through the Genesis4 file format, which requires
! nx = ny and dx = dy, and on a square grid with equal spacings the propagation kernel is
! symmetric under interchanging the two transverse axes. A transposed transform is therefore
! invisible to it however asymmetric the test field is -- confirmed by mutation: swapping
! the two arguments of fftw_plan_dft_2d passes the Python comparison at 6e-16. Only a
! rectangular grid can see it, and the reference propagator is what supplies a trustworthy
! answer on one.
!
! Invariants printed in both modes: the total field energy, which the propagator conserves
! exactly because its kernel has unit modulus, and the intensity weighted transverse rms
! sizes.
!-

program wavefront_drift_test

use wavefront_hdf5_mod

implicit none

integer n_arg
character(200) arg1

!

n_arg = command_argument_count()
if (n_arg > 0) call get_command_argument (1, arg1)

if (n_arg >= 1 .and. arg1 == '-self_check') then
  call self_check_mode()
elseif (n_arg == 3) then
  call file_drift_mode()
else
  print '(a)', 'Usage: wavefront_drift_test <input_file> <output_file> <z_drift>'
  print '(a)', '       wavefront_drift_test -self_check [tolerance]'
  stop 1
endif

contains

!------------------------------------------------------------------------------
!+
! Read a Genesis4 field file, drift it, write it back out.
!-

subroutine file_drift_mode ()

type (wavefront_struct) wf, wf_pol
integer n_grid(3)
real(rp) z_drift, energy0, energy1, sigma_x0, sigma_y0, sigma_x1, sigma_y1
real(rp) mean_x0, mean_y0, mean_x1, mean_y1, d_pol, d_ex
logical err
character(200) in_file, out_file
character(40) z_str

!

call get_command_argument (1, in_file)
call get_command_argument (2, out_file)
call get_command_argument (3, z_str)
read (z_str, *) z_drift

call wavefront_read_genesis4 (wf, in_file, err)
if (err) stop 1

n_grid = wavefront_shape(wf)

print '(a)',            'Fortran wavefront_drift_test'
print '(a, a)',         '  Input file:            ', trim(in_file)
print '(a, 3i7)',       '  Grid (nx, ny, nz):     ', n_grid
print '(a, 3es16.8)',   '  dx, dy, dz [m]:        ', wf%dx, wf%dy, wf%dz
print '(a, es16.8)',    '  Wavelength [m]:        ', wf%wavelength
print '(a, es16.8)',    '  Ref position [m]:      ', wf%ref_position
print '(a, es16.8)',    '  Drift [m]:             ', z_drift

call wavefront_transverse_moments (wf, mean_x0, mean_y0, sigma_x0, sigma_y0)
energy0 = wavefront_energy(wf)

! A second copy carrying both polarisation components, to check two things at once: that
! free-space propagation leaves Ex and Ey equal if they started equal, and that adding a
! second component does not perturb the first. The second is the stronger statement, since
! it would catch a shared work buffer being reused across components.

wf_pol = wf
wf_pol%Ey = wf_pol%Ex

call wavefront_drift (wf, z_drift, err);      if (err) stop 1
call wavefront_drift (wf_pol, z_drift, err);  if (err) stop 1

call wavefront_transverse_moments (wf, mean_x1, mean_y1, sigma_x1, sigma_y1)
energy1 = wavefront_energy(wf)

! Peak normalized, so that these do not read as large errors when everything involved is
! small.

d_pol = maxval(abs(wf_pol%Ey - wf_pol%Ex)) / maxval(abs(wf_pol%Ex))
d_ex  = maxval(abs(wf_pol%Ex - wf%Ex))     / maxval(abs(wf%Ex))

print '(a, 2es16.8)',   '  Energy before, after [J]:  ', energy0, energy1
print '(a, es16.8)',    '  Energy relative change:    ', (energy1 - energy0) / energy0
print '(a, 2es16.8)',   '  sigma_x before, after [m]: ', sigma_x0, sigma_x1
print '(a, 2es16.8)',   '  sigma_y before, after [m]: ', sigma_y0, sigma_y1
print '(a, es16.8)',    '  Max relative |Ey - Ex| with both components: ', d_pol
print '(a, es16.8)',    '  Max relative Ex change from adding Ey:       ', d_ex

if (d_pol /= 0) then
  print '(a)', '  FAIL: Ex and Ey diverged, but free-space propagation must not couple them.'
  stop 1
endif

if (d_ex /= 0) then
  print '(a)', '  FAIL: adding an Ey component changed Ex.'
  stop 1
endif

call wavefront_write_genesis4 (wf, out_file, err, 'x')
if (err) stop 1

print '(a, a)',         '  Output file:           ', trim(out_file)

! Machine readable lines for validate_drift.py, which recomputes each of these from the
! output file and checks it. Without this the Fortran side of the mirrored quantities --
! wavefront_energy, wavefront_transverse_moments and the coordinate vectors they rest on --
! would be printed but never verified against anything.

print '(a)', ''
print '(a, es24.16)', 'CHECK photon_energy ', wavefront_photon_energy(wf)
print '(a, es24.16)', 'CHECK energy        ', energy1
print '(a, es24.16)', 'CHECK mean_x        ', mean_x1
print '(a, es24.16)', 'CHECK mean_y        ', mean_y1
print '(a, es24.16)', 'CHECK sigma_x       ', sigma_x1
print '(a, es24.16)', 'CHECK sigma_y       ', sigma_y1

call wavefront_fft_free()

end subroutine file_drift_mode

!------------------------------------------------------------------------------
!+
! Compare wavefront_drift against wavefront_drift_reference on rectangular grids.
!-

subroutine self_check_mode ()

real(rp) tolerance, worst, d_rel
integer i_case
logical err, failed
character(40) tol_str

! Grid sizes chosen so that nx /= ny -- which is the whole point -- and so that the four
! combinations of even and odd in the two transverse directions are all covered, since the
! negative wavenumber half of the FFT ordering is indexed differently for even and odd n.
! Sizes are small because the reference propagator is a direct sum.

integer, parameter :: n_case = 4
integer, parameter :: nx_case(n_case) = [24, 15, 18, 21]
integer, parameter :: ny_case(n_case) = [40, 22, 27, 16]

!

tolerance = 1.0e-13_rp
if (command_argument_count() >= 2) then
  call get_command_argument (2, tol_str)
  read (tol_str, *) tolerance
endif

print '(a)',         'Fortran wavefront_drift against wavefront_drift_reference'
print '(a, es12.4)', '  Tolerance: ', tolerance
print '(a)',         ''

worst = 0
failed = .false.

do i_case = 1, n_case
  call one_self_check (nx_case(i_case), ny_case(i_case), d_rel, err)
  if (err) stop 1
  worst = max(worst, d_rel)
enddo

print '(a)',         ''
print '(a, es16.8)', '  LARGEST RELATIVE DIFFERENCE (FFT against direct DFT): ', worst

if (worst > tolerance) then
  print '(a)', '  FAIL'
  stop 1
endif

print '(a)', '  PASS'

call wavefront_fft_free()

end subroutine self_check_mode

!------------------------------------------------------------------------------
!+
! One rectangular grid case: build a field, drift it both ways, compare.
!-

subroutine one_self_check (nx, ny, d_rel, err)

type (wavefront_struct) wf_fft, wf_ref
integer nx, ny
integer, parameter :: nz = 3
real(rp) d_rel, dx, dy, dz, wavelength, z_drift, scale
real(rp) energy0, energy1, e_rel
logical err

!

! Unequal transverse spacings as well as unequal grid sizes, so that a swapped dx and dy is
! caught alongside a swapped nx and ny.

dx = 2.0e-6_rp
dy = 3.5e-6_rp
dz = 1.0e-9_rp
wavelength = 1.0e-9_rp
z_drift = 0.02_rp

call wavefront_init (wf_fft, nx, ny, nz, dx, dy, dz, wavelength, 'xy')
call fill_test_field (wf_fft)

energy0 = wavefront_energy(wf_fft)

wf_ref = wf_fft

call wavefront_drift (wf_fft, z_drift, err);            if (err) return
call wavefront_drift_reference (wf_ref, z_drift, err);  if (err) return

energy1 = wavefront_energy(wf_fft)
e_rel = (energy1 - energy0) / energy0

! Peak normalized, per the rule that a difference is only meaningful against the largest
! value in the transform rather than against a local one that may be near zero.

scale = max(maxval(abs(wf_ref%Ex)), maxval(abs(wf_ref%Ey)))
d_rel = max(maxval(abs(wf_fft%Ex - wf_ref%Ex)), maxval(abs(wf_fft%Ey - wf_ref%Ey))) / scale

print '(a, i4, a, i4, a, es14.6, a, es14.6)', &
      '  nx, ny = ', nx, ',', ny, '   max relative difference = ', d_rel, &
      '   energy relative change = ', e_rel

end subroutine one_self_check

!------------------------------------------------------------------------------
!+
! Fill a wavefront with a deterministic asymmetric test field.
!
! Asymmetry is the point, and every feature below breaks a symmetry that would otherwise
! let a bug through. The centroid is off axis in both planes and differently in each, so an
! x and y interchange shows. The amplitude carries a cross term, so the ellipse is tilted
! and a transpose shows even where the two widths happen to agree. The linear phase gives a
! nonzero mean transverse wavenumber, so the centroid moves during the drift and the sign of
! the kernel matters. The quadratic phase has opposite signs in x and y, so the two axes
! evolve differently. Everything varies with slice index, so a mistake on the slice axis
! cannot cancel. Ey is a different field from Ex, not a copy.
!-

subroutine fill_test_field (wf)

type (wavefront_struct) wf
integer n_grid(3), nx, ny, nz, ix, iy, iz
real(rp) k0, w_ref, f, x0, y0, sx, sy, u, v, amp, phase
real(rp), allocatable :: xvec(:), yvec(:)

!

n_grid = wavefront_shape(wf)
nx = n_grid(1); ny = n_grid(2); nz = n_grid(3)

xvec = wavefront_xvec(wf)
yvec = wavefront_yvec(wf)
k0 = wavefront_k0(wf)
w_ref = 4 * wf%dx    ! Transverse scale: resolved by the grid and well inside it.

do iz = 1, nz
  f = (iz - 1.0_rp) / max(nz - 1, 1)

  x0 = (0.9_rp + 0.4_rp * f) * w_ref
  y0 = -1.3_rp * w_ref
  sx = (1.0_rp + 0.5_rp * f) * w_ref
  sy = 1.7_rp * w_ref

  do iy = 1, ny
    do ix = 1, nx
      u = (xvec(ix) - x0) / sx
      v = (yvec(iy) - y0) / sy
      amp = (0.4_rp + f) * exp(-(u**2 + v**2 + 0.6_rp * u * v) / 2)
      phase = k0 * (2.5e-5_rp * (xvec(ix) - x0) - 1.1e-5_rp * (yvec(iy) - y0)) &
            + k0 * ((xvec(ix) - x0)**2 / 8 - (yvec(iy) - y0)**2 / 14) + 0.37_rp * iz
      wf%Ex(ix, iy, iz) = amp * cmplx(cos(phase), sin(phase), wf_rp)

      ! Ey: a different field, so that the two components cannot be confused for each other.
      wf%Ey(ix, iy, iz) = 0.6_rp * amp * cmplx(cos(2 * phase + 1), sin(2 * phase + 1), wf_rp)
    enddo
  enddo
enddo

end subroutine fill_test_field

end program wavefront_drift_test
