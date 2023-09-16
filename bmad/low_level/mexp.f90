!+
! Function mexp (x, m) result (this_exp)
!
! Returns x^m with 0^0 = 1.
!
! Input:
!   x -- Real(rp): Number.
!   m -- Integer: Exponent.
!
! Output:
!   this_exp -- Real(rp): Result.
!-

function mexp (x, m) result (this_exp)

use sim_utils

implicit none

real(rp) x, this_exp
integer m

!

if (m < 0) then
  this_exp = 0
elseif (m == 0) then
  this_exp = 1
else
  this_exp = x**m
endif

end function mexp
