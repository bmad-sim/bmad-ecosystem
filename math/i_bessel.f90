!+
! Function I_bessel(m, arg) result (i_bes)
!
! Function to evaluate the modified bessel function of the 
! first kind I.
!
! Modules needed:
!   use sim_utils
!
! Input:
!   m    -- Integer: Bessel order.
!   arg  -- Real(rp): Bessel argument.
!
! Output:
!   i_bes -- Real: Bessel value.
!-

function I_bessel(m, arg) result (i_bes)

use physical_constants
use fgsl

implicit none

integer m
real(rp) arg, i_bes

!

select case(m)
case (0);     i_bes = fgsl_sf_bessel_ic0(arg)
case (1);     i_bes = fgsl_sf_bessel_ic1(arg)
case default; i_bes = fgsl_sf_bessel_icn(m, arg)
end select

end function I_bessel

