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

type cylinder_coords
  complex(rp) :: rho = 0, phi = 0, z = 0
end type

type (cylinder_coords), allocatable :: E(:,:)  ! phi, z

complex(rp) expi, ez_coef, erho_coef, ephi_coef
complex(rp), allocatable :: Ez_fft(:), Ephi_fft(:), Erho_fft(:), e0(:), f0(:)

real(rp) radius, e_mag, e_phase(3), e_re(3), e_im(3), b_mag, b_phase(3), b_re(3), b_im(3), rdummy4(4)
real(rp) rad_in, x, y, z, k_z, k_l, dz, k_zz
real(rp) kappa2_l, kappa_l, freq_l, r_hat(2), m_cos_max(0:m_max), m_sin_max(0:m_max)
real(rp), allocatable :: z_pos(:), phi_pos(:)

integer i, j, m, ios, n_phi, n_z, n2_z

logical m_cos_mode(0:m_max), m_sin_mode(0:m_max) 

character(100) file_name, field_file, line

namelist / params / field_file, freq_l

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

k_l = twopi * freq_l / c_light

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

! Now read in the data

n_phi = 1
n_z = i
n2_z = 2**ceiling(0.99999 * log(real(n_z)) / log(2.0))
allocate (E(n_phi,n2_z), z_pos(n_z), phi_pos(n_phi), e0(n2_z), f0(n2_z))
allocate (Ez_fft(n2_z), Erho_fft(n2_z), Ephi_fft(n2_z))

print *, 'Number of data points in z:', n_z
print *, 'Array size for FFT:        ', n2_z

rewind(1)
read (1, '(a)') line  ! Skip header

do i = 1, n_z
  read (1, *, iostat = ios) e_mag, e_phase, e_re, e_im, b_mag, b_phase, b_re, b_im, rdummy4, x, y, z 
  r_hat = [x, y] / sqrt(x**2 + y**2)
  z_pos(i) = z
  
  E(1,i)%rho = cmplx(dot_product(e_re(1:2), r_hat), dot_product(e_im(1:2), r_hat))
  E(1,i)%phi = cmplx(e_re(2) * r_hat(1) - e_re(1) * r_hat(2), e_im(2) * r_hat(1) - e_im(1) * r_hat(2))
  E(1,i)%z   = cmplx(e_re(3), e_im(3))
enddo

close(1)

! Error check

dz = (z_pos(n_z) - z_pos(1)) / (n_z - 1)
do i = 2, n_z
  if (abs(z_pos(i) - z_pos(1) - (i - 1) * dz) > 2e-5) then
    print *, 'Z POSITIONS NOT EVENLY SPACED!', i, z_pos(i) - z_pos(1), (i - 1) * dz
    call err_exit
  endif
enddo

! Analyze the phi dependence

do m = 1, m_max
  if (modulo(n_phi, max(1, 4*m)) /= 0) then
    print *, 'Note: Not considering modes with m =', m
    cycle
  endif

  E_cos_tot = 0; E_sin_tot = 0
  do i = 1, n2_z

    c = 0; s = 0
    do j = 1, n_phi
      c = c + E(j, i)%z * cos(m * (j - 1) * twopi / n_phi)
      s = s + E(j, i)%z * sin(m * (j - 1) * twopi / n_phi)
    enddo

    if (c > 0) then
      E_cos_tot = E_cos_tot + c / (2 * n_phi)
      E_sin_tot = E_sin_tot + s / (2 * n_phi)
    else
      E_cos_tot = E_cos_tot - c / (2 * n_phi)
      E_sin_tot = E_sin_tot - s / (2 * n_phi)
    endif

  enddo

  phase(m) = modulo2 (atan2(E_sin_tot, E_cos_tot), pi/(2*m))
  phase(-m) = phase(m) + pi / (2*m)
  do i = 1, n2_z
    c = 0, c2 = 0
    do j = 1, n_phi
      c = c + E(j, i)%z * cos(m * (j - 1) * twopi / n_phi - phase(m))
      c2 = c2 + E(j, i)%z * cos(m * (j - 1) * twopi / n_phi - phase(-m))
    enddo

    E_zm( m,i) = c / (2 * n_phi)
    E_zm(-m,i) = c2 / (2 * n_phi)
  enddo

enddo

E_zm = 0

E_zm(0, :) = sum(E(:,:)%z, 1) / n_phi

! Choose largest modes

Ez_resid = E%z
do m = -m_max, m_max
  E_zm_amp = 0
  do i = 1, n2_z
    do j = 1, n_phi
      dE = E_zm(m, i) * cos(m * (j - 1) * twopi / n_phi - phase(m))
      E_zm_amp(m) = E_zm_amp(m) + abs(dE)
      Ez_resid = Ez_resid - dE
    enddo
  enddo
  E_zm_amp(m) = E_zm_amp(m) / (n2_z * n_phi)
enddo

print *
print *, '  m, <|E_zm|>'
do m = -m_max, m_max
  print *, m, E_zm_amp(m)
enddo

print *
print *, '<|E_z|>           ', sum(abs(E%z)) / (n2_z * n_phi)
print *, '<|E_z - E_z(fit)|>', sum(abs(Ez_resid)) / (n2_z * n_phi)

! Find coefficients

do m = -m_max, m_max

  

  Ez_fft(i) = E_zm(m)sum(E(:,:)%z, 1) / n_phi
  call four1(Ez_fft, -1)

  Ephi_fft(i) = sum(E(:,:)%phi, 1) / n_phi
  call four1(Ephi_fft, -1)

  do i = 1, n2_z
    k_z = twopi * (i - 1) / (n2_z * dz)  
    if (2 * i > n2_z) k_z = k_z - twopi / dz
    kappa2_l = k_z**2 - k_l**2
    kappa_l = sqrt(abs(kappa2_l))
  
    e0(i) = Ez_fft(i) / (R(0, kappa2_l, kappa_l, radius) * n2_z)
    f0(i) = Ephi_fft(i) / (R(1, kappa2_l, kappa_l, radius) * n2_z)
  enddo
endif

! m > 0

do m = 1, m_max
  Ez_fft(i) = sum(E(:,:)%z, 1) / n_phi
  call four1(Ez_fft, -1)

  Ephi_fft(i) = sum(E(:,:)%phi, 1) / n_phi
  call four1(Ephi_fft, -1)

  do i = 1, n2_z
    k_z = twopi * (i - 1) / (n2_z * dz)  
    if (2 * i > n2_z) k_z = k_z - twopi / dz
    kappa2_l = k_z**2 - k_l**2
    kappa_l = sqrt(abs(kappa2_l))
  
    e0(i) = Ez_fft(i) / (R(0, kappa2_l, kappa_l, radius) * n2_z)
    f0(i) = Ephi_fft(i) / (R(1, kappa2_l, kappa_l, radius) * n2_z)
  enddo

  call err_exit
enddo

!--------------------------------------------------------------------
! Check at full radius with full FFT

Ez_fft = 0
Erho_fft = 0
Ephi_fft = 0

do i = 1, n2_z
  k_z = twopi * (i - 1) / (n2_z * dz)
  if (2 * i > n2_z) k_z = k_z - twopi / dz
  kappa2_l = k_z**2 - k_l**2
  kappa_l = sqrt(abs(kappa2_l))

  Ez_coef = e0(i) * R(0, kappa2_l, kappa_l, radius)
  Erho_coef = (-i_imaginary * k_z / kappa_l) * e0(i) * R(1, kappa2_l, kappa_l, radius)
  Ephi_coef = f0(i) * R(1, kappa2_l, kappa_l, radius) 

  do j = 1, n2_z
    expi = cmplx(cos(dz*k_z*(j-1)), sin(dz*k_z*(j-1)))
    Ez_fft(j) = Ez_fft(j) + Ez_coef * expi
    Erho_fft(j) = Erho_fft(j) + Erho_coef * expi
    Ephi_fft(j) = Ephi_fft(j) + Ephi_coef * expi
  enddo
enddo

open (1, file = 'check30_fft', recl = 200)
do i = 1, n2_z
  write (1, '(i6, 6(2es13.3, 2x))') i, real(Ez_fft(i)), real(E(1,i)%z), aimag(Ez_fft(i)), aimag(E(1,i)%z), &
                         real(Erho_fft(i)), real(E(1,i)%rho), aimag(Erho_fft(i)), aimag(E(1,i)%rho), &
                         real(Ephi_fft(i)), real(E(1,i)%phi), aimag(Ephi_fft(i)), aimag(E(1,i)%phi) 
enddo
close (1)

! Check at full radius

do i = 1, n2_z
  k_z = twopi * (i - 1) / (n2_z * dz)  
  if (2 * i > n2_z) k_z = k_z - twopi / dz
  kappa2_l = k_z**2 - k_l**2
  kappa_l = sqrt(abs(kappa2_l))

  Ez_fft(i) = e0(i) * R(0, kappa2_l, kappa_l, radius)
  Erho_fft(i) = (-i_imaginary * k_z / kappa_l) * e0(i) * R(1, kappa2_l, kappa_l, radius)
  Ephi_fft(i) = f0(i) * R(1, kappa2_l, kappa_l, radius) 
enddo

call four1(Ez_fft, 1)
call four1(Erho_fft, 1)
call four1(Ephi_fft, 1)

open (1, file = 'check30', recl = 200)
open (2, file = 'fundamental_r30.csv', status = 'old')
read (2, *) line

do i = 1, n2_z
  if (i > n_z) then
    e_re = 0
    e_im = 0  
  else
    read (2, *, iostat = ios) e_mag, e_phase, e_re, e_im, b_mag, b_phase, b_re, b_im, rdummy4, x, y, z 
  endif
  write (1, '(i6, 6(2es13.3, 2x))') i, real(Ez_fft(i)), e_re(3), aimag(Ez_fft(i)), e_im(3), &
                         real(Erho_fft(i)), e_re(2), aimag(Erho_fft(i)), e_im(2), &
                         real(Ephi_fft(i)), -e_re(1), aimag(Ephi_fft(i)), -e_im(1)
enddo
close (1)
close (2)

! Check at half radius

do i = 1, n2_z
  k_z = twopi * (i - 1) / (n2_z * dz)  
  if (2 * i > n2_z) k_z = k_z - twopi / dz
  kappa2_l = k_z**2 - k_l**2
  kappa_l = sqrt(abs(kappa2_l))

  Ez_fft(i) = e0(i) * R(0, kappa2_l, kappa_l, radius/2)
  Erho_fft(i) = (-i_imaginary * k_z / kappa_l) * e0(i) * R(1, kappa2_l, kappa_l, radius/2)
  Ephi_fft(i) = f0(i) * R(1, kappa2_l, kappa_l, radius/2) 
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
  write (1, '(i6, 6(2es13.3, 2x))') i, real(Ez_fft(i)), e_re(3), aimag(Ez_fft(i)), e_im(3), &
                         real(Erho_fft(i)), e_re(2), aimag(Erho_fft(i)), e_im(2), &
                         real(Ephi_fft(i)), -e_re(1), aimag(Ephi_fft(i)), -e_im(1)
enddo
close (1)
close (2)

! Check on centerline

do i = 1, n2_z
  k_z = twopi * (i - 1) / (n2_z * dz)  
  if (2 * i > n2_z) k_z = k_z - twopi / dz
  kappa2_l = k_z**2 - k_l**2
  kappa_l = sqrt(abs(kappa2_l))

  Ez_fft(i) = e0(i) * R(0, kappa2_l, kappa_l, 0.0_rp)
  Erho_fft(i) = (-i_imaginary * k_z / kappa_l) * e0(i) * R(1, kappa2_l, kappa_l, 0.0_rp)
  Ephi_fft(i) = f0(i) * R(1, kappa2_l, kappa_l, 0.0_rp) 
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
  write (1, '(i6, 6(2es13.3, 2x))') i, real(Ez_fft(i)), e_re(3), aimag(Ez_fft(i)), e_im(3), &
                         real(Erho_fft(i)), e_re(2), aimag(Erho_fft(i)), e_im(2), &
                         real(Ephi_fft(i)), -e_re(1), aimag(Ephi_fft(i)), -e_im(1)
enddo
close (1)
close (2)

!------------------------------------------------------------------
contains

function R(m, kappa2_l, kappa_l, rho) result (r_out)

integer m
real(rp) rho, r_out
real(rp) kappa2_l, kappa_l

!

if (m == 0) then
  if (kappa2_l > 0) then
    r_out = bessi0(kappa_l * rho)
  else
    r_out = bessj0(kappa_l * rho)
  endif
elseif (m == 1) then
  if (kappa2_l > 0) then
    r_out = bessi1(kappa_l * rho)
  else
    r_out = bessj1(kappa_l * rho)
  endif
else
  if (kappa2_l > 0) then
    r_out = bessi(m, kappa_l * rho)
  else
    r_out = bessj(m, kappa_l * rho)
  endif
endif

end function

end program
