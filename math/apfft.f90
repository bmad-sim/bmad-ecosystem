!+
! module apfft_mod
!
! This module implements the All Phase FFT method for obtaining accurate phase from signal data.
!-

module apfft_mod

use physical_constants

implicit none

contains

!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
!+
! subroutine apfft_corr(rdata_in, bounds, window, phase, amp, freq, diag)
!
! For real signal rdata_in, computes phase, frequency, and amplitude
! of peak found within bounds.  Algorithm is corrected all-phase FFT and should.
!
! This routine finds only one peak:  the largest amplitude within the bound.  Signals with multiple
! components can be investigated by varying bounds appropriately.
! 
! Input:
!   rdata_in(:)               -- real(rp): signal data.
!   bounds(1:2)               -- real(rp): range within which to search for peak.
!   window                    -- character(3): 'rec' or 'han' for rectangular or Hann window.
!   diag                      -- integer, optional:  causes low-level routine apfft_ext to produce
!                                                    a fort.X file where X=9000+fid containing diag information.
! Output:
!   phase                     -- real(rp): phase of peak found in signal.
!   freq                      -- real(rp): frequency of peak
!   amp                       -- real(rp): amplitude of peak
!-

subroutine apfft_corr(rdata_in, bounds, window, phase, amp, freq, diag)

use, intrinsic :: ieee_arithmetic

implicit none

real(rp) rdata_in(:)
real(rp) phase, amp, freq
real(rp), optional :: bounds(2)
character(3) window
integer, optional :: diag

real(rp), allocatable :: rdata(:)
real(rp), allocatable :: u1(:), u2(:)
real(rp) phase_u1, amp_u1, freq_u1
real(rp) phase_u2, amp_u2, freq_u2
real(rp) d

integer i,N, ndata

character(*), parameter :: r_name = 'apfft_corr'

ndata = size(rdata_in)
N = floor((ndata+1.0)/3)
ndata = 3*N-1

allocate(rdata(ndata))

!remove zero offset
rdata = rdata_in(1:ndata) - sum(rdata_in(1:ndata))/size(rdata_in(1:ndata))

allocate(u1(2*N-1))
allocate(u2(2*N-1))
u1 = rdata(1:2*N-1)
u2 = rdata(N+1:3*N-1)

call apfft_ext(u1, bounds, window, phase_u1, amp_u1, freq_u1, diag)
call apfft_ext(u2, bounds, window, phase_u2, amp_u2, freq_u2, diag)

d = (phase_u2-phase_u1)/twopi
if(d .gt. 0.5) then
  d = d - 1.0d0
elseif(d .le. -0.5) then
  d = d + 1.0d0
endif

freq = freq_u1 + (1.0d0*d)/N

phase = 2.0d0*phase_u1 - phase_u2
if(phase .lt. 0) then
  phase = phase + twopi
elseif(phase .gt. twopi) then
  phase = phase - twopi
endif

if(window == 'rec') then
  if( ieee_is_finite(1.0d0/sin(pi*d)) ) then
    amp = (pi*d/sin(pi*d))**2 * amp_u1 * 2.0
  else
    amp = 2.0d0 * amp_u1
  endif
elseif(window == 'han') then
  if( ieee_is_finite(1.0d0/sin(pi*d)) ) then
    amp = (pi*d*(1.0d0-d*d)/sin(pi*d))**2 * amp_u1  * 2.0
  else
    amp = 2.0d0 * amp_u1
  endif
endif

deallocate(rdata)
deallocate(u1)
deallocate(u2)

end subroutine apfft_corr

!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
!+
! subroutine apfft(rdata_in, bounds, window, phase, diag)
!
! Implements the All Phase FFT method for obtaining accurate phase from signal data.
!
! The signal data is truncated to an odd length, and the phase is relative to the central point.
!
! Input:
!   rdata_in(:)      - real(rp): signal data.
!   bounds(2)        - real(rp): Fractional.  Boundaries within which peak is searched.
!   window           - character(3): 'rec' or 'han'
!   diag             - integer, optional :: If present and positive, produced diag file containing spectrum.
!
! Output:
!   
!   phase            - real(rp), Phase of found peak, as determined by ApFFT.
!-

subroutine apfft(rdata_in, bounds, window, phase, diag)

use fgsl

implicit none

real(rp) rdata_in(:)
real(rp) bounds(2)
character(3) window
real(rp) phase

real(rp) amp, freq
real(rp) rdata(size(rdata_in))
integer, optional :: diag

rdata = rdata_in - 1.0d0*sum(rdata_in)/size(rdata_in)

call apfft_ext(rdata, bounds, window, phase, amp, freq, diag)

end subroutine apfft

!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
!+
! subroutine apfft_ext(rdata,bounds, window, phase, amp, freq, diag)
!
! Implements the All Phase FFT method for obtaining accurate phase from signal data.
!
! This "extended" apfft subroutine returns the amplitudes and frequency as well, for use
! by the corrected apfft subroutine in this module.
!
! Input:
!   rdata_in(:)      - real(rp): signal data.
!
! Output:
!   phase            - real(rp), Phase dominant peak, as determined by ApFFT.
!-

subroutine apfft_ext(rdata, bounds, window, phase, amp, freq, diag)
use fgsl

implicit none

real(rp) rdata(:)
real(rp) bounds(2)
character(3) window
real(rp) phase, amp, freq
integer, optional :: diag

real(rp) wc(size(rdata))

integer fid
integer ndata, N, max_ix
integer i, j, lb, ub

real(rp) alpha, beta, gamma, p

complex(rp), allocatable :: c_Xap(:)

type(fgsl_fft_complex_wavetable) :: wavetable
type(fgsl_fft_complex_workspace) :: work
integer(fgsl_int) :: status
real(rp), allocatable :: apwindow(:)

character(4) fid_str
character(*), parameter :: r_name = 'apfft_ext'

ndata = size(rdata)
N = floor((ndata+1.0d0)/2.0d0)

allocate(c_Xap(N))

!Produce the all-phase vector from signal data.

if(window == 'rec') then
  !rec-rec window
  c_Xap(1) = cmplx(N*rdata(N),0.0d0)
  do i=2,N
    c_Xap(i) = cmplx((N-i+1)*rdata(N+i-1) + (i-1)*rdata(i-1),0.0d0)
  enddo
  c_Xap = c_Xap / N / N
elseif( window == 'han') then
  allocate(apwindow(ndata))
  call hanhan(N,apwindow)
  rdata = apwindow*rdata
  deallocate(apwindow)
  c_Xap(1) = cmplx(rdata(N),0.0d0)
  do i=2,N
    c_Xap(i) = cmplx(rdata(N+i-1) + rdata(i-1),0.0d0)
  enddo
  c_Xap = c_Xap 
else
  write(*,*) "Unknown window passed to apfft_ext: ", window
  call err_exit()
endif

wavetable = fgsl_fft_complex_wavetable_alloc(int(N,fgsl_size_t))
work = fgsl_fft_complex_workspace_alloc(int(N,fgsl_size_t))
status = fgsl_fft_complex_backward(c_Xap, 1_fgsl_size_t, int(N,fgsl_size_t), wavetable, work)
call fgsl_fft_complex_workspace_free(work)
call fgsl_fft_complex_wavetable_free(wavetable)

if(present(diag)) then
  if( diag .gt. 0) then
    fid = 9000 + diag
    write(fid_str,'(i4)') fid
    open(fid, file = 'apfft_'//fid_str//'.diag')
    write(fid,*) '#apfft_ext diag file'
    write(fid,*) '#freq   abs   angle'
    do i=1,N
      write(fid,'(3f15.6)') (i-1.0d0)/N, abs(c_Xap(i)), twopi_atan2(aimag(c_Xap(i)),real(c_Xap(i)))
    enddo
    write(fid,*)
    write(fid,*)
    close(fid)
  endif
endif

lb = nint(N*bounds(1)) + 1
ub = nint(N*bounds(2)) + 1
max_ix = maxloc(abs(c_Xap(lb:ub)),1) + lb - 1

phase = twopi_atan2(aimag(c_Xap(max_ix)),real(c_Xap(max_ix)))
amp = abs(c_Xap(max_ix))
freq = (max_ix-1.0d0)/N

deallocate(c_Xap)

!----------------------------------------------------------------------

contains

function twopi_atan2(Y,X) result(arg)
  use physical_constants
  real(rp) Y, X, arg
  arg = atan2(Y,X)
  if(arg < 0) arg = arg+twopi
  arg = twopi - arg
end function twopi_atan2

end subroutine apfft_ext

!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------

subroutine hanhan(N, hh)
integer N
real(rp) hh(:) !length must be 2*N-1

real(rp) h(N)
integer i

!

h = han(N)

!convolution of two Hann windows
do i=1,N
  hh(i) = sum( h(N-i+1:N) * h(1:i) )
  hh(2*N-i) = hh(i)
enddo
hh = hh 
end subroutine hanhan

!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------

function han(N) 
implicit none
integer i, N
real(rp) han(N)

!

do i=1,N
  han(i) = sin(pi*i/(N-1.0d0))**2 / (N-1.0d0) * 2.0d0
enddo
end function han

end module







