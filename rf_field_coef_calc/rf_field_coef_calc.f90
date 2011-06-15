!+
! Program to calculate the coefficients to specify the electric and magnetic field
! in an RF cavity for a given mode given field data on the surface of a cylinder 
! running though the cavity.
!
! These coefficients can then be used for particle tracking.
! See:
!    Dan Abell, PRST-AB 9, 052001 (2006)
!   "Numerical computation of high-order transfer maps for rf cavities."
!-

program rf_field_coef_calc

use bmad
use nr

implicit none

integer, parameter :: m_max = 4  ! Maximum mode m value considered.

type (rf_field_mode_struct), target :: mode(-m_max:m_max)
type (rf_field_mode_struct), pointer :: md


type field_cylinder_values
  complex(rp) :: rho = 0, phi = 0, z = 0
end type

!type cylinder_coords
! real(rp) :: rho = 0, phi = 0, z = 0
!end type

type (field_cylinder_values), allocatable :: E_dat(:,:)  ! phi, z
type (field_cylinder_values) E_here

complex(rp) expi, ez_coef, erho_coef, ephi_coef, E_cos, E_sin, E_zm, kappa_n
complex(rp), allocatable :: Ez_fft(:), Ephi_fft(:), Erho_fft(:), Ez_resid(:,:)

real(rp) e_mag, e_phase(3), e_re(3), e_im(3), b_mag, b_phase(3), b_re(3), b_im(3), rdummy4(4)
real(rp) radius, rad_in, x, y, z, k_z, k_t, dz, k_zz, cos_amp, sin_amp, z_here, amp_E_z_dat
real(rp) kappa2_n, kap_rho, freq_t, r_hat(2), amp_E_zm(-m_max:m_max)
real(rp), allocatable :: z_pos(:), phi_pos(:), cos_term(:)
real(rp) ez_dat_rms, ez_fit_rms, length


integer i, j, m, n, ios, n_phi, n_z, n2_z

logical significant_mode(-m_max:m_max)

character(100) file_name, field_file, line

namelist / params / field_file, freq_t

! Read in the parameters

if (cesr_iargc() > 1) then
  print *, 'Usage:'
  print *, '  rf_field_coef_calc {<input_file>}'
  stop
endif

file_name = 'rf_field_coef_calc.init'
if (cesr_iargc() == 1) call cesr_getarg(1, file_name)
print *, 'Input file: ', trim(file_name)

open (1, file = file_name, status = 'old')
read (1, nml = params)
close (1)

mode(:)%freq = freq_t
k_t = twopi * freq_t / c_light

! Read in the field values...
! First count the number of points

open (1, file = field_file, status = 'old')
read (1, '(a)') line  ! Skip header

i = 0
radius = -1
do 
  read (1, *, iostat = ios) e_mag, e_phase, e_re, e_im, b_mag, b_phase, b_re, b_im, rdummy4, x, y, z 
  if (ios > 0) then
    print *, 'Read error in field file on line', i+3
    call err_exit
  endif
  if (ios < 0) exit
  i = i + 1
  rad_in = sqrt(x**2 + y**2)
  if (radius < 0) radius = rad_in
  if (abs(radius - rad_in) > 1e-10) then
    print *, 'Radius mismatch', radius, rad_in
    call err_exit
  endif
enddo

print *, 'Sample_radius:', radius

! Now read in the data

n_phi = 1
n_z = i
n2_z = 2**ceiling(0.99999 * log(real(n_z)) / log(2.0))

allocate (E_dat(n_phi,n2_z), z_pos(n_z), phi_pos(n_phi), Ez_resid(n_phi, n2_z), cos_term(n_phi))
allocate (Ez_fft(n2_z), Erho_fft(n2_z), Ephi_fft(n2_z))

do m = -m_max, m_max
  allocate(mode(m)%term(n2_z))
  mode(m)%m = abs(m)
enddo

print *, 'Number of data points in z:', n_z
print *, 'Array size for FFT:        ', n2_z

rewind(1)
read (1, '(a)') line  ! Skip header

do i = 1, n_z
  read (1, *, iostat = ios) e_mag, e_phase, e_re, e_im, b_mag, b_phase, b_re, b_im, rdummy4, x, y, z 
  r_hat = [x, y] / sqrt(x**2 + y**2)
  z_pos(i) = z
  
  E_dat(1,i)%rho = cmplx(dot_product(e_re(1:2), r_hat), dot_product(e_im(1:2), r_hat))
  E_dat(1,i)%phi = cmplx(e_re(2) * r_hat(1) - e_re(1) * r_hat(2), e_im(2) * r_hat(1) - e_im(1) * r_hat(2))
  E_dat(1,i)%z   = cmplx(e_re(3), e_im(3))
enddo

close(1)

! Error check

length = z_pos(n_z) - z_pos(1)
dz = length / (n_z - 1)
mode%dz = dz

do i = 2, n_z
  if (abs(z_pos(i) - z_pos(1) - (i - 1) * dz) > 2e-5) then
    print *, 'Z POSITIONS NOT EVENLY SPACED!', i, z_pos(i) - z_pos(1), (i - 1) * dz
    call err_exit
  endif
enddo

! Analyze the phi dependence to get the best phi0.
! mode(-m) is out-of-phase with mode(m) and if there is only a single mode present, it
! will be assocaited with mode(m) and mode(-m) will have zero amplitude.

do m = 1, m_max
  if (modulo(n_phi, max(1, 4*m)) /= 0) then
    print '(a, i0)', 'Note: Due to number of points in phi, will not consider modes with m = ', m
    cycle
  endif

  cos_amp = 0; sin_amp = 0
  do i = 1, n2_z

    E_cos = 0; E_sin = 0
    do j = 1, n_phi
      E_cos = E_cos + E_dat(j, i)%z * cos(m * (j - 1) * twopi / n_phi)
      E_sin = E_sin + E_dat(j, i)%z * sin(m * (j - 1) * twopi / n_phi)
    enddo

    if (real(E_cos) > 0) then
      cos_amp = cos_amp + real(E_cos)
      sin_amp = sin_amp + real(E_sin)
    else
      cos_amp = cos_amp - real(E_cos)
      sin_amp = sin_amp - real(E_sin)
    endif

    if (aimag(E_cos) > 0) then
      cos_amp = cos_amp + aimag(E_cos)
      sin_amp = sin_amp + aimag(E_sin)
    else
      cos_amp = cos_amp - aimag(E_cos)
      sin_amp = sin_amp - aimag(E_sin)
    endif

  enddo

  mode(m)%phi_0 = modulo2 (atan2(sin_amp, cos_amp), pi/(2*m))
  mode(-m)%phi_0 = mode(m)%phi_0 + pi / (2*m) ! Out-of-phase with mode(m)

enddo

! Compute mode amplitudes and choose largest modes

Ez_resid = E_dat%z
amp_E_z_dat = sum(abs(E_dat%z)) / (n2_z * n_phi)
significant_mode = .false.
amp_E_zm = 0

do m = -m_max, m_max

  if (modulo(n_phi, max(1, 4*abs(m))) /= 0) cycle

  do j = 1, n_phi
    cos_term(j) = cos(m * (j - 1) * twopi / n_phi - mode(m)%phi_0)
  enddo

  do i = 1, n2_z
    if (m == 0) then
      E_zm = sum(E_dat(:, i)%z * cos_term) / n_phi
    else
      E_zm = sum(E_dat(:, i)%z * cos_term) / (2 * n_phi)
    endif
    amp_E_zm(m) = amp_E_zm(m) + abs(E_zm) * sum(abs(cos_term)) / n2_z
    Ez_resid(:, i) = Ez_resid(:, i) - E_zm * cos_term
  enddo

  significant_mode(m) = (amp_E_zm(m) > 0.01 * Amp_E_z_dat) 

enddo

print *
print *, '  m      <|E_zm|>'
do m = -m_max, m_max
  print '(i4, f14.4, a, l)', m, amp_E_zm(m), '  ', significant_mode(m)
enddo

print *
print *, '<|E_z|>           ', amp_E_z_dat
print *, '<|E_z - E_z(fit)|>', sum(abs(Ez_resid)) / (n2_z * n_phi)
print *, '<|E_phi|>         ', sum(abs(E_dat%phi)) / (n2_z * n_phi)
print *, '<|E_rho|>         ', sum(abs(E_dat%rho)) / (n2_z * n_phi)

!--------------------------------------------------
! Find coefficients for m == 0

if (significant_mode(0)) then

  Ez_fft(:) = sum(E_dat(:,:)%z, 1) / n_phi
  call four1(Ez_fft, -1)

  Ephi_fft(:) = sum(E_dat(:,:)%phi, 1) / n_phi
  call four1(Ephi_fft, -1)

  do i = 1, n2_z
    k_z = twopi * (i - 1) / (n2_z * dz)  
    if (2 * i > n2_z) k_z = k_z - twopi / dz
    kappa2_n = k_z**2 - k_t**2
    kappa_n = sqrt(abs(kappa2_n))
    kap_rho = kappa_n * radius
    if (kappa2_n < 0) then
      kappa_n = -i_imaginary * kappa_n
      kap_rho = -kap_rho
    endif
  
    mode(0)%term(i)%e = Ez_fft(i) / (I_bessel(0, kap_rho) * n2_z)
    mode(0)%term(i)%b = Ephi_fft(i) * kappa_n / (I_bessel(1, kap_rho) * n2_z)
  enddo
endif

! Find coefficients for m /= 0

do m = -m_max, m_max

  if (m == 0) cycle
  if (.not. significant_mode(m)) cycle

  do j = 1, n_phi
    cos_term(j) = cos(m * (j - 1) * twopi / n_phi - mode(m)%phi_0)
  enddo

  do i = 1, n2_z
    Ez_fft(i)   = sum(E_dat(:, i)%z * cos_term) / (2 * n_phi)
    Erho_fft(i) = sum(E_dat(:, i)%rho * cos_term) / (2 * n_phi)
  enddo

  call four1(Ez_fft, -1)
  call four1(Erho_fft, -1)

  do i = 1, n2_z
    k_z = twopi * (i - 1) / (n2_z * dz)  
    if (2 * i > n2_z) k_z = k_z - twopi / dz
    kappa2_n = k_z**2 - k_t**2
    kappa_n = sqrt(abs(kappa2_n))
    kap_rho = kappa_n * radius
    if (kappa2_n < 0) then
      kappa_n = -i_imaginary * kappa_n
      kap_rho = -kap_rho
    endif
  
    mode(m)%term(i)%e = Ez_fft(i) * kappa_n**m / (I_bessel(m, kap_rho) * n2_z)
    mode(m)%term(i)%b = (radius * kappa_n**m / I_bessel(m, kap_rho)) * &
          (i_imaginary * Erho_fft(i) / n2_z - &
           k_z * mode(m)%term(i)%e * I_bessel(m+1, kap_rho)/kappa_n**(m+1))
  enddo

enddo

!--------------------------------------------------------------------
! Write coefs
! Assume just one m = 0 mode for now.

open (1, file = 'rf_field_coef.dat')

write (1, '(a, f10.5, a)') 'l = ', length, ','
write (1, '(a)') 'rf_field = {'

do m = -m_max, m_max
  if (.not. significant_mode(m)) cycle
  md => mode(m)

  write (1, '(2x, a)') 'mode = {' 
  write (1, '(4x, a, i0, a)')     'm             =', abs(m),           ','
  write (1, '(4x, a, es11.4, a)') 'freq          =', md%freq,          ','
  write (1, '(4x, a, f0.1, a)')   'f_damp        =', 0.0_rp,           ','
  write (1, '(4x, a, f0.1, a)')   'theta_t0      =', 0.0_rp,           ','
  write (1, '(4x, a, f0.4, a)')   'phi_0         =', md%phi_0,         ','
  write (1, '(4x, a, f0.6, a)')   'dz            =', md%dz,            ','
  write (1, '(4x, a, f0.1, a)')   'stored_energy =', 0.0_rp,           ','
  
  n = size(md%term)

  write (1, '(4x, a)') 'e_re = ('
  do i = 1, (n - 1) / 10
    write (1, '(6x, 10(es12.4, a))') (real(md%term(j)%e), ',', j = (i-1)*10+1, i*10)
  enddo
  write (1, '(6x, 10(es12.4, a))') (real(md%term(j)%e), ',', j = (i-1)*10+1, n-1), real(md%term(n)%e), '),'

  write (1, '(4x, a)') 'e_im = ('
  do i = 1, (n - 1) / 10
    write (1, '(6x, 10(es12.4, a))') (aimag(md%term(j)%e), ',', j = (i-1)*10+1, i*10)
  enddo
  write (1, '(6x, 10(es12.4, a))') (aimag(md%term(j)%e), ',', j = (i-1)*10+1, n-1), aimag(md%term(n)%e), ')}}'

enddo

close(1)

!--------------------------------------------------------------------
! Write field files as a check.

Ez_fft = 0
Erho_fft = 0
Ephi_fft = 0

open (1, file = 'check30_fft', recl = 200)

do i = 1, n2_z
  z_here = (i-1) * dz
  E_here = e_field_calc (radius, 0.0_rp, z_here, mode, significant_mode)
  Ez_fft(i)   = E_here%z
  Erho_fft(i) = E_here%rho
  Ephi_fft(i) = E_here%phi

  write (1, '(f7.4, 6(2es13.3, 2x), i4)') z_here, real(Ez_fft(i)), real(E_dat(1,i)%z), aimag(Ez_fft(i)), aimag(E_dat(1,i)%z), &
                         real(Erho_fft(i)), real(E_dat(1,i)%rho), aimag(Erho_fft(i)), aimag(E_dat(1,i)%rho), &
                         real(Ephi_fft(i)), real(E_dat(1,i)%phi), aimag(Ephi_fft(i)), aimag(E_dat(1,i)%phi), i
enddo
close (1)

!----------------------------------
! Check at full radius

do i = 1, n2_z
  k_z = twopi * (i - 1) / (n2_z * dz)  
  if (2 * i > n2_z) k_z = k_z - twopi / dz
  kappa2_n = k_z**2 - k_t**2
  kappa_n = sqrt(abs(kappa2_n))
  kap_rho = kappa_n * radius
  if (kappa2_n < 0) then
    kappa_n = -i_imaginary * kappa_n
    kap_rho = -kap_rho
  endif

  Ez_fft(i) = mode(0)%term(i)%e * I_bessel(0, kap_rho)
  Erho_fft(i) = -i_imaginary * k_z * mode(0)%term(i)%e * (I_bessel(1, kap_rho) / kappa_n)
  Ephi_fft(i) = mode(0)%term(i)%b * (I_bessel(1, kap_rho) / kappa_n)
enddo

call four1(Ez_fft, 1)
call four1(Erho_fft, 1)
call four1(Ephi_fft, 1)

open (1, file = 'check30', recl = 200)
open (2, file = 'fundamental_r30.csv', status = 'old')
read (2, *) line

ez_dat_rms = 0
ez_fit_rms = 0

do i = 1, n2_z
  if (i > n_z) then
    e_re = 0
    e_im = 0  
  else
    read (2, *, iostat = ios) e_mag, e_phase, e_re, e_im, b_mag, b_phase, b_re, b_im, rdummy4, x, y, z 
  endif
  z_here = (i-1) * dz
  write (1, '(f7.4, 6(2es13.3, 2x), i4)') z_here, real(Ez_fft(i)), e_re(3), aimag(Ez_fft(i)), e_im(3), &
                         real(Erho_fft(i)), e_re(2), aimag(Erho_fft(i)), e_im(2), &
                         real(Ephi_fft(i)), -e_re(1), aimag(Ephi_fft(i)), -e_im(1), i
  ez_dat_rms = ez_dat_rms + e_re(3)**2
  ez_fit_rms = ez_fit_rms + real(Ez_fft(i))**2
enddo

print *
print *, 'RMS Dat:', sqrt(ez_dat_rms / n2_z)
print *, 'RMS Fit:', sqrt(ez_fit_rms / n2_z)

close (1)
close (2)

!----------------------------------
! Check at half radius

do i = 1, n2_z
  k_z = twopi * (i - 1) / (n2_z * dz)  
  if (2 * i > n2_z) k_z = k_z - twopi / dz
  kappa2_n = k_z**2 - k_t**2
  kappa_n = sqrt(abs(kappa2_n))
  kap_rho = kappa_n * radius
  if (kappa2_n < 0) then
    kappa_n = -i_imaginary * kappa_n
    kap_rho = -kap_rho
  endif

  Ez_fft(i) = mode(0)%term(i)%e * I_bessel(0, kap_rho/2)
  Erho_fft(i) = -i_imaginary * k_z * mode(0)%term(i)%e * I_bessel(1, kap_rho/2) / kappa_n
  Ephi_fft(i) = mode(0)%term(i)%b * I_bessel(1, kap_rho/2) / kappa_n
enddo

call four1(Ez_fft, 1)
call four1(Erho_fft, 1)
call four1(Ephi_fft, 1)

open (1, file = 'check15', recl = 200)
open (2, file = 'fundamental_r15.csv', status = 'old')
read (2, *) line

do i = 1, n2_z
  if (i > n_z) then
    e_re = 0
    e_im = 0  
  else
    read (2, *, iostat = ios) e_mag, e_phase, e_re, e_im, b_mag, b_phase, b_re, b_im, rdummy4, x, y, z 
  endif
  z_here = (i-1) * dz
  write (1, '(f7.4, 6(2es13.3, 2x), i4)') z_here, real(Ez_fft(i)), e_re(3), aimag(Ez_fft(i)), e_im(3), &
                         real(Erho_fft(i)), e_re(2), aimag(Erho_fft(i)), e_im(2), &
                         real(Ephi_fft(i)), -e_re(1), aimag(Ephi_fft(i)), -e_im(1), i
enddo
close (1)
close (2)

!----------------------------------
! Check on centerline

do i = 1, n2_z
  k_z = twopi * (i - 1) / (n2_z * dz)  
  if (2 * i > n2_z) k_z = k_z - twopi / dz
  kappa2_n = k_z**2 - k_t**2
  kappa_n = sqrt(abs(kappa2_n))
  kap_rho = kappa_n * radius
  if (kappa2_n < 0) then
    kappa_n = -i_imaginary * kappa_n
    kap_rho = -kap_rho
  endif

  Ez_fft(i) = mode(0)%term(i)%e * I_bessel(0, 0.0_rp)
  Erho_fft(i) = -i_imaginary * k_z * mode(0)%term(i)%e * I_bessel(1, 0.0_rp) / kappa_n
  Ephi_fft(i) = mode(0)%term(i)%b * I_bessel(1, 0.0_rp) / kappa_n
enddo

call four1(Ez_fft, 1)
call four1(Erho_fft, 1)
call four1(Ephi_fft, 1)

open (1, file = 'check00', recl = 200)
open (2, file = 'fundamental_r00.csv', status = 'old')
read (2, *) line

do i = 1, n2_z
  if (i > n_z) then
    e_re = 0
    e_im = 0  
  else
    read (2, *, iostat = ios) e_mag, e_phase, e_re, e_im, b_mag, b_phase, b_re, b_im, rdummy4, x, y, z 
  endif
  z_here = (i-1) * dz
  write (1, '(f7.4, 6(2es13.3, 2x), i4)') z_here, real(Ez_fft(i)), e_re(3), aimag(Ez_fft(i)), e_im(3), &
                         real(Erho_fft(i)), e_re(2), aimag(Erho_fft(i)), e_im(2), &
                         real(Ephi_fft(i)), -e_re(1), aimag(Ephi_fft(i)), -e_im(1), i
enddo
close (1)
close (2)

!------------------------------------------------------------------
contains

function e_field_calc (rho, phi, z, modes, use_mode) result (E)

type (rf_field_mode_struct), target :: modes(:)
type (rf_field_mode_struct), pointer :: mode
type (field_cylinder_values) E

complex(rp) ikz, expi
complex(rp) Im_minus, Im0, Im_plus, Im_norm, kappa_n

real(rp) rho, phi, z, k_t, kappa2_n, kap_rho, c, s

integer i, imode

logical :: use_mode(:)

!

E%rho = 0; E%phi = 0; E%z = 0

do imode = 1, size(modes)

  if (.not. use_mode(imode)) cycle
  mode => modes(imode)

  k_t = twopi * mode%freq / c_light

  do i = 1, size(mode%term)
    k_z = twopi * (i - 1) / (size(mode%term) * mode%dz)
    if (2 * i > n2_z) k_z = k_z - twopi / dz
    kappa2_n = k_z**2 - k_t**2
    kappa_n = sqrt(abs(kappa2_n))
    kap_rho = kappa_n * radius
    if (kappa2_n < 0) then
      kappa_n = -i_imaginary * kappa_n
      kap_rho = -kap_rho
    endif
    expi = cmplx(cos(k_z * z), sin(k_z * z))

    if (mode%m == 0) then
      Im0     = I_bessel(0, kap_rho)
      Im_plus = I_bessel(1, kap_rho) / kappa_n
      E%rho = E%rho + -i_imaginary * k_z * mode%term(i)%e * Im_plus * expi
      E%phi = E%phi + mode%term(i)%b * Im_plus * expi
      E%z   = E%z   + mode%term(i)%e * Im0 * expi

    else
      c = cos(mode%m * phi - mode%phi_0)
      s = sin(mode%m * phi - mode%phi_0)
      Im_plus  = I_bessel(mode%m+1, kap_rho)
      Im_minus = I_bessel(mode%m-1, kap_rho) 

      Im_norm  = kappa_n * (Im_minus - Im_plus) / (2 * mode%m) / kappa_n**m ! I_m(kap_rho) / rho
      Im0      = rho * Im_norm

      Im_plus = Im_plus / kappa_n**(m+1)
      Im_minus = Im_minus / kappa_n**(m-1)

      E%rho = E%rho - i_imaginary * (k_z * mode%term(i)%e * Im_plus + mode%term(i)%b * Im_norm) * c * expi
      E%phi = E%phi - i_imaginary * (k_z * mode%term(i)%e * Im_plus + &
                                             mode%term(i)%b * (Im_norm - Im_minus / mode%m)) * s * expi
      E%z   = E%z   + mode%term(i)%e * Im0 * c * expi

    endif

  enddo

enddo

end function

end program
