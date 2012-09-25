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
real(rp), parameter :: f(0:30) = [ &
      1.000000000000D+00, 1.000000000000D+00, 2.000000000000D+00, 6.000000000000D+00, &
      2.400000000000D+01, 1.200000000000D+02, 7.200000000000D+02, 5.040000000000D+03, &
      4.032000000000D+04, 3.628800000000D+05, 3.628800000000D+06, 3.991680000000D+07, &
      4.790016000000D+08, 6.227020800000D+09, 8.717829120000D+10, 1.307674368000D+12, &
      2.092278988800D+13, 3.556874280960D+14, 6.402373705728D+15, 1.216451004088D+17, &
      2.432902008177D+18, 5.109094217171D+19, 1.124000727778D+21, 2.585201673888D+22, &
      6.204484017332D+23, 1.551121004333D+25, 4.032914611266D+26, 1.088886945042D+28, &
      3.048883446117D+29, 8.841761993740D+30, 2.652528598122D+32]

integer n, i

! Error check

if (n < 0) then
  print *, 'ERROR IN FACTORIAL(N). N < 0!'
  if (global_com%exit_on_error) call err_exit
endif

! Use Sterling's formula if n is very large

if (n > ubound(f, 1)) then
  fact = exp(n * log(real(n, rp)) - n + log(twopi * n) / 2)
else
  fact = f(n)
endif

end function
