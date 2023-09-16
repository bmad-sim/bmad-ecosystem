!+
! Function I_bessel_extended(m, arg) result (i_bes)
!
! Function to evaluate the "extended" bessel function which is a 
! combination of the modified bessel function I and the bessel function J:
!   I_bes_ext = I_m(arg)                         For arg > 0
!   I_bes_ext = I_m(-i*arg) = i^{-m} * J(arg)    For arg < 0
!
! Also see:
!   I_bessel
!
! Input:
!   m    -- Integer: Bessel order.
!   arg  -- Real(rp): Bessel argument.
!
! Output:
!   i_bes -- Complex(rp): Bessel value.
!-

function I_bessel_extended(m, arg) result (i_bes)

use sim_utils_struct
use fgsl

implicit none

integer m
real(rp) arg
complex(rp) i_bes

!

select case(m)
case (0)
  if (arg > 0) then
    i_bes = fgsl_sf_bessel_ic0(arg)
  else
    i_bes = fgsl_sf_bessel_jc0(-arg)
  endif

case (1)
  if (arg > 0) then
    i_bes = fgsl_sf_bessel_ic1(arg)
  else
    i_bes = -i_imaginary * fgsl_sf_bessel_jc1(-arg)
  endif

case default
  if (arg > 0) then
    i_bes = fgsl_sf_bessel_icn(m, arg)
  else
    i_bes = (-i_imaginary)**m * fgsl_sf_bessel_jcn(m, -arg)
  endif
end select

end function I_bessel_extended

