!+
! Function I_bessel(m, arg) result (i_bes)
!
! Function to evaluate the modified bessel function I:
!   I_bes = I_m(arg)                         For arg > 0
!   I_bes = I_m(-i*arg) = i^{-m} * J(arg)    For arg < 0
!
! Modules needed:
!   use sim_utils
!
! Input:
!   m    -- Integer: Bessel order.
!   arg  -- Real(rp): Bessel argument.
!
! Output:
!   i_bes -- Complex(rp): Bessel value.
!-

function I_bessel(m, arg) result (i_bes)

use physical_constants
use nr

integer m
real(rp) arg
complex(rp) i_bes

!

select case(m)
case (0)
  if (arg > 0) then
    i_bes = bessi0(arg)
  else
    i_bes = bessj0(-arg)
  endif

case (1)
  if (arg > 0) then
    i_bes = bessi1(arg)
  else
    i_bes = -i_imaginary * bessj1(-arg)
  endif

case default
  if (arg > 0) then
    i_bes = bessi(m, arg)
  else
    i_bes = (-i_imaginary)**m * bessj(m, -arg)
  endif
end select

end function I_bessel

