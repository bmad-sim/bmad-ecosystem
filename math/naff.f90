!+
! module naff_mod
!
! This module implements the NAFF algorithm for calculating the spectra
! of periodic data.
!
! This module is useful when you need an accurate frequency decomposition from only a small number of samples.  
!
! Precision of the determined freqiencies goes as 1/N^4, where N is the number of data points.
!
! Decomposes complex spectrum data of the form D(:) = D1(:) + i D2(:).
!
! freqs contains the frequencies found in the data.
! amps contains the complex amplitudes of these frequencies.
!
! If opt_dump_spectra=<some integer> the FFT spectra will be dumped to a fort.<some integer> file.
! This will also cause other debug information to be printed to stdout.
!
! If opt_zero_first is present and .true., then the first component returned will be frequency=0.
!
! The steps of NAFF are:
! 1) Estimate peak omega_1 in frequency spectrum using interpolated FFT.
! 2) Refine estimate by using optimizer to maximize <data|e^{-i omega_1}>, and also return amplitude of this component.
! 3) Remove e_1 = amp*e^{i omega_i} component from the data.
! 4) Repeat step 1 to estimate the new frequency component.
! 5) Repeat step 2 to refine the new frequency component.
! 6) Use Gram-Schmidt to build a basis function by orthogonalizing the new frequency component relative to the existing frequency components.
! 7) Delete the orthogonalized basis function from the data.
! 8) Repeat at step 4 until new frequency components are no longer significant.
!-

module naff_mod

use sim_utils

implicit none

contains

!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
!+
! subroutine naff(cdata,freqs,amps,opt_dump_spectra,opt_zero_first)
!
! This subroutine implements the NAFF algorithm for calculating the spectra
! of periodic data.
!
! See naff_mod documentation for details.
!
! Frequencies returned are in units of 2pi. That is, freqs ranges from 0 to 1.
!
! freqs and amps must be allocated before hand.  This subroutine will repeat the
! decomposition loop until all elements of freqs and amps are populated.
!
! Input:
!   cdata(:)         - complex(rp), complex signal data. Size must be power of 2.
!   opt_dump_spectra - integer, optional, If present write FFT spectra to file named "fort.N"  
!                       Also debug info printed on the terminal
!   opt_zero_first   - logical, optional: If present and true, then the first component
!                       returned will be freq = 0.
!
! Output:
!   freqs(:)         - real(rp), frequency components found in units of 0 to 1.
!                       Will set to -1 if there is an error [Can happen, for example, if all data is zero 
!                       except for one point.]
!   amps(:)          - complex(rp), amplitudes of frequency components.
!-

subroutine naff(cdata, freqs, amps, opt_dump_spectra, opt_zero_first)

implicit none

complex(rp) cdata(:)
real(rp) freqs(:)
complex(rp) amps(:)
integer, optional :: opt_dump_spectra
logical, optional :: opt_zero_first

complex(rp), allocatable :: u(:,:)

integer n_comp, size_data
integer i, j
logical calc_ok
real(rp) throw_away

character(*), parameter :: r_name = 'naff'

!

amps = 0

if (all(cdata == cdata(1))) then
  freqs = 0
  return
endif

size_data = size(cdata)
n_comp = size(freqs)

allocate(u(n_comp, size_data))

freqs = -1

do i=1,n_comp
  if(i==1 .and. logic_option(.false., opt_zero_first)) then
    !throw_away = interpolated_fft(cdata, calc_ok, opt_dump_spectra, 0) !call to get spectrum dump
    throw_away = interpolated_fft_gsl(cdata, calc_ok, opt_dump_spectra, 0) !call to get spectrum dump
    freqs(i) = 0.0d0
  else
    !freqs(i) = interpolated_fft(cdata, calc_ok, opt_dump_spectra, 0)  !estimate location of spectrum peak using FFT
    freqs(i) = interpolated_fft_gsl(cdata, calc_ok, opt_dump_spectra, 0)  !estimate location of spectrum peak using FFT
    if (.not. calc_ok) return
    freqs(i) = maximize_projection(freqs(i), cdata)  !refine location of frequency peak
  endif
  u(i,:) = ed(freqs(i),size_data)
  do j=2,i
    u(i,:) = u(i,:) - projdd(u(j-1,:),u(i,:))*u(j-1,:)
  enddo
  amps(i) = projdd(u(i,:),cdata)
  cdata = cdata - amps(i)*u(i,:)
enddo

end subroutine naff

!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------

function projdd(a,b)

complex(rp) a(:), b(:) 
complex(rp) projdd
projdd = sum(conjg(a)*b)/size(a)

end function

!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------

function ed(freq,N)

real(rp) freq
integer N,t
complex(rp) ed(1:N)
ed = (/ (exp((0.0d0,-1.0d0)*twopi*freq*t),t=0,N-1) /)

end function

!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
!+
! function maximize_projection
!
! Optimizer that uses Numerical Recipes brent to find a local maximum, 
! which is the frequency that maximizes the projection.
!
!-

function maximize_projection(seed, cdata)

use super_recipes_mod, only: super_mnbrak, super_brent

implicit none

real(rp) seed
complex(rp) cdata(:)
real(rp) maximize_projection

integer i, N, status
real(rp) small_step
real(rp) :: tol = 1.0d-8, abs_tol = 1.0d-15
real(rp) fmin !throw away
real(rp) ax,bx,cx
real(rp) fa,fb,fc

real(rp) r, window, hsamp

complex(rp) wcdata(size(cdata))

N = size(cdata)
small_step = 0.1d0/N

!apply window
hsamp = (N-1)/2.0d0
do i=1, N
  window = 0.5d0*(1.0d0 + cos(twopi*(i-hsamp-1)/(N-1)))  !hanning
  ! r = 8.0d0
  ! window = exp( -0.5d0*(r*(1.0d0*i-hsamp-1)/(N-1))**2)  !gaussian
  wcdata(i)= cdata(i) * window
enddo

ax = seed
bx = seed+small_step
cx = 0.0d0
call super_mnbrak(ax, bx, cx, fa, fb, fc, special_projection, status)
fmin = super_brent(ax, bx, cx, special_projection, tol, abs_tol, maximize_projection, status)

!--------------------------------------------

contains
!+
! function special_projection
!
! Calculates <cdata | exp(i theta)>
! 
! Used only by maximize projection. Uses data global to the function to accomodate stock NR routine.
!-

function special_projection (f, status)
  integer, optional :: status
  real(rp), intent(in) :: f
  real(rp) special_projection
  special_projection = -1.0d0*abs(projdd(wcdata,ed(f,N)))
end function special_projection

end function maximize_projection

!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
!+
!  function interpolated_fft_gsl
!
!  Windows the complex data and uses a mixed-radix GSL routine to find the peak in the spectrum.
!  The result is interpolated to improve the accuracy.  Hanning and Gaussian windowing are
!  available.
!-

function interpolated_fft_gsl (cdata, calc_ok, opt_dump_spectrum, opt_dump_index) result (this_fft)

use fgsl
use sim_utils

complex(rp) cdata(:)
integer, optional :: opt_dump_spectrum, opt_dump_index

integer dump_spectrum, dump_index

complex(rp) wcdata(size(cdata))
real(rp) fft_amp(size(cdata))
real(rp) this_fft
real(rp) window, r, hsamp, denom
real(rp) lk, lkm, lkp, A

type(fgsl_fft_complex_wavetable) :: wavetable
type(fgsl_fft_complex_workspace) :: work
integer(fgsl_int) :: status

integer n_samples
integer i, max_ix

logical calc_ok

!

calc_ok = .false.
this_fft = -1

dump_spectrum = integer_option (0, opt_dump_spectrum)
dump_index = integer_option (0, opt_dump_index)

n_samples = size(cdata)
hsamp = (n_samples-1)/2.0d0

!apply window
do i=1, n_samples
  !window = 0.5d0*(1.0d0 + cos(twopi*(i-hsamp-1)/(n_samples-1)))  !hanning, also adjust interpolation below
  r = 8.0
  window = exp( -0.5d0*(r*(1.0d0*i-hsamp-1.0d0)/(n_samples-1.0d0))**2)  !gaussian, also adjust interpolation below
  ! window = 1
  wcdata(i)= cdata(i) * window
enddo

wavetable = fgsl_fft_complex_wavetable_alloc(int(n_samples,fgsl_size_t))
work = fgsl_fft_complex_workspace_alloc(int(n_samples,fgsl_size_t))
status = fgsl_fft_complex_backward(wcdata, 1_fgsl_size_t, int(n_samples,fgsl_size_t), wavetable, work)
fft_amp(:)=sqrt(wcdata(:)*conjg(wcdata(:)))

if( dump_spectrum > 10 ) then
  do i=1,n_samples
    write(dump_spectrum,*) dump_index, (i-1.0d0)/n_samples, fft_amp(i)
  enddo
  write(dump_spectrum,*)
  write(dump_spectrum,*)
endif

max_ix = maxloc(fft_amp(2:n_samples-1), 1) + 1
if (fft_amp(max_ix) == 0) return  ! Error

!Gaussian Interpolation (use with gaussian window)
lk = log(fft_amp(max_ix))
lkm = log(fft_amp(max_ix-1))
lkp = log(fft_amp(max_ix+1))
denom = 2.0d0*lk - lkp - lkm
if (denom == 0) return   ! Error

A = (lkp-lkm) / 2.0d0 / denom
this_fft = ( 1.0d0*(max_ix-1) + A ) / n_samples
calc_ok = .true.

! Parabolic Interpolation (use with hanning window)
! lk = fft_amp(max_ix)
! lkm = fft_amp(max_ix-1)
! lkp = fft_amp(max_ix+1)
! A = (lkp-lkm) / 2.0 / (2.0*lk - lkp - lkm)
! this_fft = 1.0d0*max_ix/n_samples + A/n_samples

call fgsl_fft_complex_workspace_free(work)
call fgsl_fft_complex_wavetable_free(wavetable)

end function interpolated_fft_gsl

!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
!+
!  Function interpolated_fft (cdata, calc_ok, opt_dump_spectrum, opt_dump_index) result (this_fft)
!
!  Windows the complex data and used Numerical Recipes four1 to find the peak in the spectrum.
!  The result is interpolated to improve the accuracy.  Hanning and Gaussian windowing are
!  available.
!-

function interpolated_fft (cdata, calc_ok, opt_dump_spectrum, opt_dump_index) result (this_fft)

use sim_utils

implicit none

complex(rp) cdata(:)
integer, optional :: opt_dump_spectrum, opt_dump_index

integer dump_spectrum, dump_index

complex(rp) wcdata(size(cdata))
real(rp) fft_amp(size(cdata))
real(rp) this_fft
real(rp) window, r, hsamp, denom
real(rp) lk, lkm, lkp, A

integer n_samples, i, max_ix
integer(8) plan

logical calc_ok

!

calc_ok = .false.
this_fft = -1

dump_spectrum = integer_option (0, opt_dump_spectrum)
dump_index = integer_option (0, opt_dump_index)

n_samples = size(cdata)
hsamp = (n_samples-1)/2.0d0

!apply window
do i=1, n_samples
  !window = 0.5d0*(1.0d0 + cos(twopi*(i-hsamp-1)/(n_samples-1)))  !hanning, also adjust interpolation below
  r = 8.0
  window = exp(-0.5d0*(r*(1.0d0*i-hsamp-1.0d0)/(n_samples-1.0d0))**2)  !gaussian, also adjust interpolation below
  ! window = 1
  wcdata(i)= cdata(i) * window
enddo

call fft_1d (wcdata, -1)
fft_amp(:)=sqrt(wcdata(:)*conjg(wcdata(:)))

if( dump_spectrum > 10 ) then
  do i=1,n_samples
    write(dump_spectrum,*) dump_index, (i-1.0d0)/n_samples, fft_amp(i)
  enddo
  write(dump_spectrum,*)
  write(dump_spectrum,*)
endif

max_ix = maxloc(fft_amp(2:n_samples-1), 1) + 1
if (fft_amp(max_ix) == 0) return  ! Error

!Gaussian Interpolation (use with gaussian window)
lk = log(fft_amp(max_ix))
lkm = log(fft_amp(max_ix-1))
lkp = log(fft_amp(max_ix+1))
denom = 2.0d0*lk - lkp - lkm
if (denom == 0) return   ! Error

A = (lkp-lkm) / 2.0d0 / denom
this_fft = ( 1.0d0*(max_ix-1) + A ) / n_samples
calc_ok = .true.

! Parabolic Interpolation (use with hanning window)
! lk = fft_amp(max_ix)
! lkm = fft_amp(max_ix-1)
! lkp = fft_amp(max_ix+1)
! A = (lkp-lkm) / 2.0 / (2.0*lk - lkp - lkm)
! this_fft = 1.0d0*max_ix/n_samples + A/n_samples

end function interpolated_fft

end module







