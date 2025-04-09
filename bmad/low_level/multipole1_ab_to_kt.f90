!+
! Subroutine multipole1_ab_to_kt (an, bn, n, knl, tn)
!
! Subroutine to convert ab type multipole to kt (MAD standard) multipole.
! Also see: multipole_ab_to_kt.
!
! Input:
!   an -- Real(rp): Skew multipole component.
!   bn -- Real(rp): Normal multipole component.
!   n  -- Integer: Order of multipole. 
!
! Output:
!   knl -- Real(rp): Multitude magnatude.
!   tn  -- Real(rp): Multipole angle.
!-

subroutine multipole1_ab_to_kt (an, bn, n, knl, tn)

use sim_utils

implicit none

real(rp) an, bn, knl, tn
integer n

!

real(rp) a, b

if (an == 0 .and. bn == 0) then
  knl = 0
  tn = 0
else
  ! Use temp a, b to avoid problems when actual (knl, tn) args are the same as (an, bn).
  a = an
  b = bn
  knl  = factorial(n) * sqrt(a**2 + b**2)
  tn = -atan2(a, b) / (n + 1)
  ! In case the user looks at this, make tn to be in the range [-pi, pi]/(n+1)
  if (2 * (n + 1) * abs(tn) > pi) then
    knl = -knl
    tn = sign_of(tn) * (abs(tn) - pi/(n+1))
  endif

endif

end subroutine multipole1_ab_to_kt

