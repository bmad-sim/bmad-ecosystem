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

use rf_field_coef_calc_mod

implicit none

integer, parameter :: m_max = 4  ! Maximum mode m value considered.

type (rf_field_mode_struct), target :: mode(-m_max:m_max), abell_mode
type (rf_field_mode_struct), pointer :: md

!type cylinder_coords
! real(rp) :: rho = 0, phi = 0, z = 0
!end type

type (field_cylinder_values), allocatable :: E_dat(:,:)  ! phi, z
type (field_cylinder_values) E_here

complex(rp) expi, ez_coef, erho_coef, ephi_coef, E_cos, E_sin, E_zm, kappa_n
complex(rp), allocatable :: Ez_fft(:), Ephi_fft(:), Erho_fft(:), Ez_resid(:,:)

real(rp) e_mag, e_phase(3), e_re(3), e_im(3), b_mag, b_phase(3), b_re(3), b_im(3), rdummy4(4)
real(rp) radius, rad_in, x, y, z, k_z, k_t, dz, k_zz, cos_amp, sin_amp, z_here, amp_E_z_dat
real(rp) kappa2_n, kap_rho, freq_t, r_hat(2), amp_E_zm(-m_max:m_max), phi
real(rp), allocatable :: z_pos(:), phi_pos(:), cos_term(:)
real(rp) ez_dat_rms, ez_fit_rms, length, ez_in(4), ephi_in(4), erho_in(4)


integer i, j, im, n, ios, n_phi, n_z, n2_z, file_type, nc

logical significant_mode(-m_max:m_max)

character(100) file_name, field_file, line

namelist / params / field_file, freq_t, file_type, radius

! Read in the parameters

if (command_argument_count() > 1) then
  print *, 'Usage:'
  print *, '  rf_field_coef_calc {<input_file>}'
  stop
endif

file_name = 'rf_field_coef_calc.init'
if (command_argument_count() == 1) call get_command_argument(1, file_name)
print *, 'Input file: ', trim(file_name)

open (1, file = file_name, status = 'old')
read (1, nml = params)
close (1)

mode(:)%freq = freq_t
k_t = twopi * freq_t / c_light

! Read in the field values...
! First count the number of points

open (1, file = field_file, status = 'old')
i = 0

if (file_type == 2) then
  do
    read (1, *, iostat = ios) line
    if (ios > 0) call err_exit
    if (ios < 0) exit
    i = i + 1
  enddo

  n_phi = 4
  n_z = i    
  n2_z = 2**ceiling(0.99999 * log(real(n_z)) / log(2.0))
  
  allocate (E_dat(n_phi,n2_z), z_pos(n_z), phi_pos(n_phi), Ez_resid(n_phi, n2_z), cos_term(n_phi))
  allocate (Ez_fft(n2_z), Erho_fft(n2_z), Ephi_fft(n2_z))
  allocate(abell_mode%term(n2_z))

  do im = -m_max, m_max
    allocate(mode(im)%term(n2_z))
    mode(im)%m = abs(im)
  enddo

  rewind(1)

  do i = 1, n_z
    read (1, *) z, ez_in, ephi_in, erho_in
    z_pos(i) = z
    
    E_dat(:,i)%rho = cmplx(erho_in)
    E_dat(:,i)%phi = cmplx(ephi_in)
    E_dat(:,i)%z   = cmplx(ez_in)
  enddo

else

  read (1, '(a)') line  ! Skip header

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

  do im = -m_max, m_max
    allocate(mode(im)%term(n2_z))
    mode(im)%m = abs(im)
  enddo

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

endif

print *, 'Number of data points in z:', n_z
print *, 'Array size for FFT:        ', n2_z

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

do im = 1, m_max
  if (modulo(n_phi, max(1, 4*im)) /= 0) then
    print '(a, i0)', 'Note: Due to number of points in phi, will not consider modes with m = ', im
    cycle
  endif

  cos_amp = 0; sin_amp = 0
  do i = 1, n2_z

    E_cos = 0; E_sin = 0
    do j = 1, n_phi
      E_cos = E_cos + E_dat(j, i)%z * cos(im * (j - 1) * twopi / n_phi)
      E_sin = E_sin + E_dat(j, i)%z * sin(im * (j - 1) * twopi / n_phi)
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

  mode(im)%phi_0 = modulo2 (atan2(sin_amp, cos_amp), pi/(2*im))
  mode(-im)%phi_0 = mode(im)%phi_0 + pi / (2*im) ! Out-of-phase with mode(im)

enddo

! Compute mode amplitudes and choose largest modes

Ez_resid = E_dat%z
amp_E_z_dat = sum(abs(E_dat%z)) / (n2_z * n_phi)
significant_mode = .false.
amp_E_zm = 0

do im = -m_max, m_max

  md => mode(im)

  if (modulo(n_phi, max(1, 4*md%m)) /= 0) cycle

  do j = 1, n_phi
    cos_term(j) = cos(md%m * (j - 1) * twopi / n_phi - mode(im)%phi_0)
  enddo

  do i = 1, n2_z
    if (im == 0) then
      E_zm = sum(E_dat(:, i)%z * cos_term) / n_phi
    else
      E_zm = 2 * sum(E_dat(:, i)%z * cos_term) / n_phi
    endif
    amp_E_zm(im) = amp_E_zm(im) + abs(E_zm) * sum(abs(cos_term)) / n2_z
    Ez_resid(:, i) = Ez_resid(:, i) - E_zm * cos_term
  enddo

  significant_mode(im) = (amp_E_zm(im) > 0.01 * Amp_E_z_dat) 

enddo

print *
print *, '  m      <|E_zm|>'
do im = -m_max, m_max
  print '(i4, f14.4, a, l)', im, amp_E_zm(im), '  ', significant_mode(im)
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
  
    mode(0)%term(i)%e = Ez_fft(i) / (I_bessel_extended(0, kap_rho) * n2_z)
    mode(0)%term(i)%b = Ephi_fft(i) * kappa_n / (I_bessel_extended(1, kap_rho) * n2_z)
  enddo
endif

! Find coefficients for m /= 0

do im = -m_max, m_max

  md => mode(im)
  if (im == 0) cycle
  if (.not. significant_mode(im)) cycle

  do j = 1, n_phi
    cos_term(j) = cos(md%m * (j - 1) * twopi / n_phi - mode(im)%phi_0)
  enddo

  do i = 1, n2_z
    Ez_fft(i)   = 2 * sum(E_dat(:, i)%z * cos_term) / n_phi
    Erho_fft(i) = 2 * sum(E_dat(:, i)%rho * cos_term) / n_phi
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
  
    mode(im)%term(i)%e = Ez_fft(i) * kappa_n**md%m / (I_bessel_extended(md%m, kap_rho) * n2_z)
    mode(im)%term(i)%b = (radius * kappa_n**md%m / I_bessel_extended(md%m, kap_rho)) * &
          (i_imaginary * Erho_fft(i) / n2_z - &
           k_z * mode(im)%term(i)%e * I_bessel_extended(md%m+1, kap_rho)/kappa_n**(md%m+1))

    abell_mode%term(i)%e = mode(im)%term(i)%e * dz / kappa_n**md%m
    abell_mode%term(i)%b = mode(im)%term(i)%b * dz / kappa_n**md%m


  enddo

enddo

!--------------------------------------------------------------------
! Write coefs

open (1, file = 'rf_field_coef.dat')

write (1, '(a, f10.5, a)') 'l = ', length, ','
write (1, '(a)') 'rf_field = {'

do im = -m_max, m_max
  if (.not. significant_mode(im)) cycle
  md => mode(im)

  write (1, '(2x, a)') 'mode = {' 
  write (1, '(4x, a, i0, a)')     'm             = ', abs(im),         ','
  write (1, '(4x, a, es11.4, a)') 'freq          =', md%freq,          ','
  write (1, '(4x, a, f0.1, a)')   'f_damp        =', 0.0_rp,           ','
  write (1, '(4x, a, f0.1, a)')   'theta_t0      =', 0.0_rp,           ','
  write (1, '(4x, a, f0.4, a)')   'phi_0         =', md%phi_0,         ','
  write (1, '(4x, a, f0.6, a)')   'dz            =', md%dz,            ','
  write (1, '(4x, a, f0.1, a)')   'stored_energy =', 0.0_rp,           ','
  
  n = size(md%term)

  nc = 6

  write (1, '(4x, a)') 'e_re = ('
  do i = 1, (n - 1) / nc
    write (1, '(6x, 10(es17.9, a))') (real(md%term(j)%e), ',', j = (i-1)*nc+1, i*nc)
  enddo
  write (1, '(6x, 10(es17.9, a))') (real(md%term(j)%e), ',', j = (i-1)*nc+1, n-1), real(md%term(n)%e), '),'

  write (1, '(4x, a)') 'e_im = ('
  do i = 1, (n - 1) / nc
    write (1, '(6x, 10(es17.9, a))') (aimag(md%term(j)%e), ',', j = (i-1)*nc+1, i*nc)
  enddo
  write (1, '(6x, 10(es17.9, a))') (aimag(md%term(j)%e), ',', j = (i-1)*nc+1, n-1), aimag(md%term(n)%e), '),'

  write (1, '(4x, a)') 'b_re = ('
  do i = 1, (n - 1) / nc
    write (1, '(6x, 10(es17.9, a))') (real(md%term(j)%b), ',', j = (i-1)*nc+1, i*nc)
  enddo
  write (1, '(6x, 10(es17.9, a))') (real(md%term(j)%b), ',', j = (i-1)*nc+1, n-1), real(md%term(n)%b), '),'

  write (1, '(4x, a)') 'b_im = ('
  do i = 1, (n - 1) / nc
    write (1, '(6x, 10(es17.9, a))') (aimag(md%term(j)%b), ',', j = (i-1)*nc+1, i*nc)
  enddo
  write (1, '(6x, 10(es17.9, a))') (aimag(md%term(j)%b), ',', j = (i-1)*nc+1, n-1), aimag(md%term(n)%b), ')}}'

enddo

close(1)

!--------------------------------------------------------------------
! Write field files as a check.
! Negative sign to correct wrong sign of E_phi in input file.

if (file_type == 2) then

  open (1, file = 'e_check.out', recl = 200)
  open (2, file = 'e_check.in', recl = 200)

  do i = 1, n2_z
    do j = 1, n_phi
      phi = (j - 1) * twopi / n_phi
      z_here = (i-1) * dz
      E_here = e_field_calc (radius, phi, z_here, mode, significant_mode)
      ez_in(j)   = E_here%z 
      Erho_in(j) = E_here%rho
      Ephi_in(j) = E_here%phi
    enddo
    write (1, '(f10.4, 12es12.4)') z_pos(1) + z_here, ez_in, ephi_in, erho_in
    write (2, '(f10.4, 12es12.4)') z_pos(1) + z_here, real(E_dat(:,i)%z), -real(E_dat(:,i)%phi), real(E_dat(:,i)%rho)
  enddo
  close (1)
  close (2)

  call exit
endif


!--------------------------------------------------------------------
! Write field files as a check.

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

  Ez_fft(i) = mode(0)%term(i)%e * I_bessel_extended(0, kap_rho)
  Erho_fft(i) = -i_imaginary * k_z * mode(0)%term(i)%e * (I_bessel_extended(1, kap_rho) / kappa_n)
  Ephi_fft(i) = mode(0)%term(i)%b * (I_bessel_extended(1, kap_rho) / kappa_n)
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

  Ez_fft(i) = mode(0)%term(i)%e * I_bessel_extended(0, kap_rho/2)
  Erho_fft(i) = -i_imaginary * k_z * mode(0)%term(i)%e * I_bessel_extended(1, kap_rho/2) / kappa_n
  Ephi_fft(i) = mode(0)%term(i)%b * I_bessel_extended(1, kap_rho/2) / kappa_n
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

  Ez_fft(i) = mode(0)%term(i)%e * I_bessel_extended(0, 0.0_rp)
  Erho_fft(i) = -i_imaginary * k_z * mode(0)%term(i)%e * I_bessel_extended(1, 0.0_rp) / kappa_n
  Ephi_fft(i) = mode(0)%term(i)%b * I_bessel_extended(1, 0.0_rp) / kappa_n
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

end program
