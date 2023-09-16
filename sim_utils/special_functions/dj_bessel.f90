!+
! Function dJ_bessel(n, arg) result (dj_bes)
!
! Function to evaluate the derivative of the bessel function J_n.
!
! Input:
!   n      -- Integer: Bessel order.
!   arg    -- Real(rp): Bessel argument.
!
! Output:
!   dj_bes -- Real: Bessel value.
!-

function dJ_bessel(n, arg) result (dj_bes)

use sim_utils_interface, only: rp, j_bessel

implicit none

integer n
real(rp) arg, dj_bes

!

if (n == 0) then
  dj_bes = -J_bessel(1, arg)
else
  dj_bes = 0.5_rp * (J_bessel(n-1, arg) - J_bessel(n+1, arg))
endif

end function dJ_bessel

