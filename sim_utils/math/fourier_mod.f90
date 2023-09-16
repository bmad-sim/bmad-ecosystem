module fourier_mod

use fgsl
use precision_def
use physical_constants

contains

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function coarse_frequency_estimate(data, error) result(frequency)
!  
! Simple function to take periodic data and estimate 
! the most dominant frequency by FFT. 
! 
! Input: 
!   data(:)    -- real(rp): data to analyze. Preferably size(data) is a power of 2
!                           Otherwise the data is padded with zeros. 
!
! Output:
!   frequency  -- real(rp): Frequency corresponding to the largest FFT amplitude
!   err        -- logical : Error: not enough data. Frequency is near 0 or 0.5
!-  

function coarse_frequency_estimate(data, error) result(frequency)

implicit none

!integer(fgsl_size_t), parameter :: nsize = 128_fgsl_size_t
real(rp) :: data(:) 
real(rp), allocatable :: fdata(:) 
real(rp) :: amp, maxamp, frequency
integer :: status, i, imax, pow
integer(fgsl_size_t) :: nsize

logical :: err
logical, optional ::  error
!

nsize = size(data)

if (.not. (iand(nsize, nsize-1) == 0) ) then
  !print *, 'Error: data should be power of 2, is:', nsize, 'padding...'
  pow = ceiling(log(1.0_rp*size(data))/log(2.0_rp))
  nsize = 2**pow
  !print *, 'new size: ', nsize
  allocate(fdata(nsize))
  fdata(1:size(data)) = data
  fdata(size(data)+1:nsize) = 0.0_rp
else
  allocate(fdata(nsize))
  fdata = data
endif

status = fgsl_fft_real_radix2_transform(fdata, 1_fgsl_size_t, nsize)

err = .false.

! First case
maxamp = fdata(1)**2
imax = 1

do i = 2, nsize/2
  amp = fdata(i)**2 + fdata(nsize+2-i)**2
  if (amp > maxamp) then
    maxamp = amp 
    imax = i
  endif
end do

! Last case
i = nsize/2+1
amp = fdata(i)**2
if (amp > maxamp) then
  ! near 0.5
  maxamp = amp 
  imax = i
  err = .true.
endif

! Near zero frequency
if (imax == 1) err = .true.

!print *, 'Estimated frequency: ', (imax-1)/(1.0_rp*nsize)

if(present(error)) error = err 
frequency =  (imax-1)*1.0_rp/(nsize)

deallocate(fdata)

end function

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function fine_frequency_estimate(data) result(frequency)
!  
! Uses Laskar's method to accurately find the most dominant frequency
! A coarse estimate is first made by FFT.
! 
! Input: 
!   data(:)    -- real(rp): data to analyze
!
! Output:
!   frequency  -- real(rp): Frequency corresponding to the largest FFT amplitude
!-

function fine_frequency_estimate(data) result(frequency)

use super_recipes_mod, only: super_dbrent

implicit none

real(rp) :: data(:)
real(rp) :: frequency
real(rp) :: ftry, fbin, fmin, fmax, festimate, df, f1, f2
real(rp) :: cos_amp, sin_amp, dcos_amp, dsin_amp, omega
real(rp) :: amp, max_amp
integer :: i, imax, status
integer, parameter :: nscan = 4 !Scan points
logical :: err

!

! FFT 
ftry = coarse_frequency_estimate(data, err)
if (err) then
  frequency = ftry
  return
endif

!print *, 'FFT coarse ftry: ', ftry

fbin = 1.0_rp/size(data)
fmin = ftry - fbin
fmax = ftry + fbin

max_amp =  -super_dbrent(fmin,ftry,fmax, negative_ampsquared, negative_dampsquared, &
                                                                 1e-12_rp, fbin*1e-15_rp, frequency, status)

!---------------------------------------------------------------------------
contains

function negative_ampsquared(frequency, status) result(amp)

implicit none
real(rp) :: cos_amp, sin_amp
real(rp) :: frequency
real(rp) amp
integer, optional :: status

!

call fourier_amplitude(data, frequency, cos_amp, sin_amp)
amp = -cos_amp**2 - sin_amp**2

end function negative_ampsquared

!---------------------------------------------------------------------------
! contains

function negative_dampsquared(frequency, status) result(damp)
implicit none
real(rp) :: cos_amp, sin_amp, dcos_amp, dsin_amp
real(rp) :: frequency
real(rp) damp
integer, optional :: status

!

call fourier_amplitude(data, frequency, cos_amp, sin_amp, dcos_amp, dsin_amp)
damp = -cos_amp*dcos_amp - sin_amp*dsin_amp

end function negative_dampsquared

end function fine_frequency_estimate


!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+ 
! Subroutine fourier_amplitude(data, frequency, cos_amp, sin_amp, dcos_amp, dsin_amp)
!
! Computes cos_amp = (1/N) * sum_n=0^{N-1} data(n-1) cos(twopi*frequency*n)
!     and  sin_amp = (1/N) * sum_n=0^{N-1} data(n-1) sin(twopi*frequency*n)
!     and optionally dcos_amp = d/dfrequency cos_amp
!                    dsin_amp = d/dfrequency sin_amp
! Input: 
!   data(:)    -- real(rp): data to analyze
!   frequency  -- real(rp): frequency
! Output:
!   cos_amp    -- real(rp): cosine amplitude
!   sin_amp    -- real(rp): sine amplitude
!   dcos_amp   -- real(rp), optional: cosine amplitude derivative
!   dsin_amp   -- real(rp), optional: sine amplitude derivative
!-

subroutine fourier_amplitude(data, frequency, cos_amp, sin_amp, dcos_amp, dsin_amp)

implicit none

real(rp) :: data(:) 
real(rp) :: cos_amp, sin_amp, frequency, omega
real(rp), optional :: dcos_amp, dsin_amp
real(rp) :: d
integer :: i, n 

!

cos_amp = 0.0_rp
sin_amp = 0.0_rp

omega = frequency*twopi

n = 0

if (present(dcos_amp) .and. present(dsin_amp) ) then
  dcos_amp = 0.0_rp
  dsin_amp = 0.0_rp
  do i = lbound(data,1), ubound(data,1)
    d = data(i)*cos(omega*n)
    cos_amp = cos_amp + d
    dsin_amp = dsin_amp + n*d  
    d = data(i)*sin(omega*n)
    sin_amp = sin_amp + d
    dcos_amp = dcos_amp - n*d
    n = n+1
  end do
  dcos_amp = twopi*dcos_amp/n
  dsin_amp = twopi*dsin_amp/n
else
  do i = lbound(data,1), ubound(data,1)
    cos_amp = cos_amp + data(i)*cos(omega*n)
    sin_amp = sin_amp + data(i)*sin(omega*n)
    n = n+1
  end do
endif

cos_amp = cos_amp/n
sin_amp = sin_amp/n

end subroutine

end module
