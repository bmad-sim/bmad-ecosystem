!+
! module apfft_mod
!
! This module implements the All Phase FFT method for obtaining accurate phase from signal data.
!-

module apfft_mod

use sim_utils

implicit none

contains

!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
!+
! subroutine apfft(cdata,freqs,amps,opt_dump_spectra,opt_zero_first)
!
! Output:
!   freqs(:)         - real(rp), frequency components found in units of 0 to 1.
!                       Will set to -1 if there is an error [Can happen, for example, if all data is zero 
!                       except for one point.]
!   amps(:)          - complex(rp), amplitudes of frequency components.
!-

subroutine apfft(rdata_in, phase)
  use fgsl

  implicit none

  real(rp) rdata_in(:)
  real(rp) rdata(size(rdata_in))
  real(rp) wc(size(rdata_in))
  real(rp) phase

  integer ndata, N, max_ix
  integer i, j

  complex(rp) c_Xap(size(rdata_in))
  type(fgsl_fft_complex_wavetable) :: wavetable
  type(fgsl_fft_complex_workspace) :: work
  integer(fgsl_int) :: status

  character(*), parameter :: r_name = 'apfft'

  ndata = size(rdata_in)
  N = floor((ndata+1.0)/2)

  !remove zero offset
  rdata = rdata_in - sum(rdata_in)/size(rdata_in)

  c_Xap(1) = cmplx(N*rdata(N),0.0d0)
  do i=2,N
    c_Xap(i) = cmplx((N-i+1)*rdata(N+i-1) + (i-1)*rdata(i-1),0.0d0)
  enddo
  c_Xap = c_Xap / N

  wavetable = fgsl_fft_complex_wavetable_alloc(int(N,fgsl_size_t))
  work = fgsl_fft_complex_workspace_alloc(int(N,fgsl_size_t))
  status = fgsl_fft_complex_backward(c_Xap, 1_fgsl_size_t, int(N,fgsl_size_t), wavetable, work)
  call fgsl_fft_complex_workspace_free(work)
  call fgsl_fft_complex_wavetable_free(wavetable)

  max_ix = maxloc(abs(c_Xap(1:floor(N/2.0))),1)

  phase = twopi_atan2(aimag(c_Xap(max_ix)),real(c_Xap(max_ix)))

  contains

    function twopi_atan2(Y,X) result(arg)
      use physical_constants
      real(rp) Y, X, arg
      arg = atan2(Y,X)
      if(arg .lt. 0) arg = arg+twopi
      arg = twopi - arg
    end function
end subroutine

function han(x,N) result(hanning)
  implicit none
  integer x
  real(rp) hanning
  integer N
  if(abs(x) .le. N) then
    hanning = sin(pi*x/(N-1.0d0))**2
  else
    hanning = 0.0d0
  endif
end function

end module







