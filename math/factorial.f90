!+
! Function factorial(n) result (fact)
!
! Function to return the Factorial n!
!
! Modules needed:
!   use sim_utils
!
! Input:
!   n -- Integer: Must be non-negative
!
! Output:
!   fact -- Real(rp): n!
!-

function factorial(n) result (fact)

use physical_constants

implicit none

real(rp) fact, lnn
real(rp), save :: f(0:30)

integer n, i

logical, save :: init_needed = .true.

!

if (n < 0) then
  print *, 'ERROR IN FACTORIAL(N). N < 0!'
  call err_exit
endif

if (init_needed) then
  f(0) = 1
  do i = 1, ubound(f, 1)
    f(i) = i * f(i-1)
  enddo
endif

! Use Sterling's formula if n is very large

if (n > ubound(f, 1)) then
  fact = exp(n * log(real(n, rp)) - n + log(twopi * n) / 2)
else
  fact = f(n)
endif

end function
