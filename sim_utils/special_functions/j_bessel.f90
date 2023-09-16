!+
! Function J_bessel(n, arg) result (j_bes)
!
! Function to evaluate the bessel function J_n.
!
! Input:
!   n    -- Integer: Bessel order.
!   arg  -- Real(rp): Bessel argument.
!
! Output:
!   j_bes -- Real: Bessel value.
!-

function J_bessel(n, arg) result (j_bes)

use physical_constants
use fgsl

implicit none

integer n
real(rp) arg, j_bes
real(rp), parameter :: arg_min(0:50) = &
            [0.0_rp, 0.0_rp, 10.0_rp**(-153), 10.0_rp**(-101), 10.0_rp**(-76), 10.0_rp**(-60), &  !  0 -  5
               10.0_rp**(-50), 10.0_rp**(-43), 10.0_rp**(-37), 10.0_rp**(-33), 10.0_rp**(-29), &      !  6 - 10
               10.0_rp**(-26), 10.0_rp**(-24), 10.0_rp**(-22), 10.0_rp**(-20), 10.0_rp**(-19), &      ! 11 - 15
               10.0_rp**(-18), 10.0_rp**(-16), 10.0_rp**(-15), 10.0_rp**(-14), 10.0_rp**(-14), &      ! 16 - 20
               10.0_rp**(-13), 10.0_rp**(-12), 10.0_rp**(-12), 10.0_rp**(-11), 10.0_rp**(-10), &      ! 21 - 25
               10.0_rp**(-10), 10.0_rp**(-10), 10.0_rp**(-09), 10.0_rp**(-09), 10.0_rp**(-08), &      ! 26 - 30
               10.0_rp**(-08), 10.0_rp**(-08), 10.0_rp**(-07), 10.0_rp**(-07), 10.0_rp**(-07), &      ! 31 - 35
               10.0_rp**(-07), 10.0_rp**(-06), 10.0_rp**(-06), 10.0_rp**(-06), 10.0_rp**(-06), &      ! 36 - 40
               10.0_rp**(-05), 10.0_rp**(-05), 10.0_rp**(-05), 10.0_rp**(-05), 10.0_rp**(-05), &      ! 41 - 45
               10.0_rp**(-05), 10.0_rp**(-04), 10.0_rp**(-04), 10.0_rp**(-04), 10.0_rp**(-04)]        ! 46 - 50

! Note: The GSL j_bessel function does not properly trap underflow!

if (n <= ubound(arg_min, 1)) then
  if (abs(arg) < arg_min(n)) then
    j_bes = 0
    return
  endif
endif

select case(n)
case (0);     j_bes = fgsl_sf_bessel_jc0(arg)
case (1);     j_bes = fgsl_sf_bessel_jc1(arg)
case default; j_bes = fgsl_sf_bessel_jcn(n, arg)
end select

end function J_bessel

