!+
! Module da2_mod
!
! Module for doing two variable Differential Algebra.
!
! The coefficients for a Taylor series are stored in a matrix rather than an array.
! This is very efficient if the number of zero terms is not large.
!
! For a given matrix ta(0:N,0:N), the Taylor series is considered to be of order N.
! That is, terms ta(i,j) with i+j > N are ignored.
!-

module da2_mod

use precision_def
use output_mod

implicit none

contains

!-----------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------
!+
! Function da2_mult (ta, tb) result (tc)
!
! Routine to multiply two da2 series.
! The size of ta and tb must be the same.
!
! Input:
!   ta(0:N,0:N) -- real(rp): Input Taylor series
!   tb(0:N,0:N) -- real(rp): Input Taylor series
!
! Output:
!   tc(0:N,0:N) -- real(rp): tc = ta * tb.
!-

function da2_mult(ta, tb) result (tc)

real(rp) ta(0:,0:), tb(0:,0:), tc(0:ubound(ta,1),0:ubound(ta,2))
integer n, i1, j1, i2, j2

!

n = ubound(ta, 1)
tc = 0

do i1 = 0, n
do j1 = 0, n-i1

  do i2 = 0, n-i1-j1
  do j2 = 0, n-i1-j1-i2

    tc(i1+i2,j1+j2) = tc(i1+i2,j1+j2) + ta(i1,j1) * tb(i2,j2)

  enddo
  enddo

enddo
enddo

end function da2_mult

!-----------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------
!+
! Function da2_div (ta, tb) result (tc)
!
! Routine to divide two da2 series.
! The size of ta, and tb must be the same.
!
! Input:
!   ta(0:N,0:N) -- real(rp): Input Taylor series
!   tb(0:N,0:N) -- real(rp): Input Taylor series
!
! Output:
!   tc(0:N,0:N) -- real(rp): tc = ta / tb.
!-

function da2_div(ta, tb) result (tc)

real(rp) ta(0:,0:), tb(0:,0:), tc(0:ubound(ta,1),0:ubound(ta,2))

!

tc = da2_mult(ta, da2_inverse (tb))

end function da2_div

!-----------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------
!+
! Function da2_inverse (ta) result (ta_inv)
!
! Routine to take the inverse of a da2 series.
!
! Input:
!   ta(0:N,0:N) -- real(rp): Input Taylor series
!
! Output:
!   ta_inv(0:N,0:N) -- real(rp): ta_inv = 1 / ta
!-

function da2_inverse(ta) result (ta_inv)

real(rp) ta(0:,0:), ta_inv(0:ubound(ta,1),0:ubound(ta,2))
real(rp) eps(0:ubound(ta,1),0:ubound(ta,2)), epsn(0:ubound(ta,1),0:ubound(ta,2))
integer n, i

character(*), parameter :: r_name = 'da2_inverse'

!

n = ubound(ta, 1)

if (ta(0,0) == 0) then
  call out_io (s_error$, r_name, 'INVERSE OF TAYLOR SERIES WITH CONSTANT TERM ZERO IS NOT POSSIBLE')
  if (global_com%exit_on_error) call err_exit
  return
endif

eps = ta / ta(0,0)
eps(0,0) = 0

ta_inv = -eps
ta_inv(0,0) = 1

epsn = -eps
do i = 2, n
  epsn = -da2_mult(eps, epsn)
  ta_inv = ta_inv + epsn
enddo

ta_inv = ta_inv / ta(0,0)

end function da2_inverse

!-----------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------
!+
! Function da2_evaluate (ta, x, y) result (value)
!
! Routine to evaluate a da2 series.
!
! Input:
!   ta(0:N,0:N) -- real(rp): Input Taylor series.
!   x, y        -- real(rp): Point to evaluate at.
!
! Output:
!   value       -- real(rp): Value of series at (x, y).
!-

function da2_evaluate (ta, x, y) result (value)

real(rp) ta(0:,0:), x, y, value
real(rp) xn(0:ubound(ta,1)), yn(0:ubound(ta,2))
integer n, i, j

n = ubound(ta, 1)
xn(0:1) = [1.0_rp, x]
yn(0:1) = [1.0_rp, y]
do i = 2, n
  xn(i) = xn(i-1) * x
  yn(i) = yn(i-1) * y
enddo

value = 0
do i = 0, n
do j = 0, n-i
  value = value + ta(i,j) * xn(i) * yn(j)
enddo
enddo

end function da2_evaluate

end module
