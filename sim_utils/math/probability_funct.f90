!+
! Function probability_funct (x) result (p_func)
!
! Routine compute the normal probability function probability_funct.
!
! Defining relation:
!     probability_funct(x) = integral_{0, x} dt exp(-t^2/2) / sqrt(2*pi)
!
! Properties:
!     probability_funct(0)      = 0
!     probability_funct(\infty) = 0.5
!     probability_funct(-x)     = -probability_funct(x)
!
! In terms of the error function:
!     probability_funct(x) = erf(x/sqrt(2)) / 2
!
! Input:
!   x         -- real(rp): Function argument.
!
! Output:
!   p_func    -- real(rp): Function value.
!-

function probability_funct(x) result (p_func)

use precision_def

implicit none

real(rp) p_func, x

p_func = erf(x/sqrt(2.0)) / 2

end function
