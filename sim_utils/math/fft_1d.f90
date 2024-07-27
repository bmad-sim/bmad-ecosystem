!+
! Subroutine fft_1d (arr, isign)
!
! 1D discrete FFT of a complex array.
!   fft(k) = sum(j = 0, N-1): A(j) * Exp(isign * 2 * pi * j * k * i_imag / N)
! Where N is the number of data points. 
! The algorithm will be most efficient when N is a power of 2.
!
! Input:
!   arr(:)    -- complex(rp): Input array.
!   isign     -- integer: -1 => "Forward" transform, +1 => "Backwards" transform.
!
! Output:
!   arr(:)    -- complex(rp): FFT of array.
!-

subroutine fft_1d (arr, isign)

use output_mod, only: out_io, rp, global_com, s_error$, real_garbage$
use, intrinsic :: iso_c_binding

implicit none
include 'fftw3.f03'

complex(rp) arr(:)
integer isign
integer(8) plan
character(*), parameter :: r_name = 'fft_1d'

!

select case (isign)
case (1)
  call dfftw_plan_dft_1d(plan, size(arr), arr, arr, FFTW_BACKWARD, FFTW_ESTIMATE)
case (-1)
  call dfftw_plan_dft_1d(plan, size(arr), arr, arr, FFTW_FORWARD, FFTW_ESTIMATE)
case default
  call out_io(s_error$, r_name, 'BAD ISIGN ARGUMENT.')
  if (global_com%exit_on_error) call err_exit
  arr = real_garbage$
  return
end select

call dfftw_execute_dft(plan, arr, arr)
call dfftw_destroy_plan(plan)

end subroutine
