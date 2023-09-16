!+
! Function n_choose_k(n, k) result (nck)
!
! Function to return the N_Choose_K = N!/(K!(N-K)!)
!
! Input:
!   n, k -- Integer: Must be non-negative with n >= k.
!
! Output:
!   nck -- Real(rp): N choose K will return negative number if there is an error.
!-

function n_choose_k(n, k) result (nck)

use physical_constants
use output_mod, dummy => n_choose_k

implicit none

real(rp) nck, lnn
real(rp), parameter :: f(0:20,0:10) = reshape([ &
     1.0_rp,     -1.0_rp,     -1.0_rp,     -1.0_rp,     -1.0_rp,     -1.0_rp,     -1.0_rp,     -1.0_rp,     -1.0_rp,     -1.0_rp,     -1.0_rp,  &
     1.0_rp,      1.0_rp,     -1.0_rp,     -1.0_rp,     -1.0_rp,     -1.0_rp,     -1.0_rp,     -1.0_rp,     -1.0_rp,     -1.0_rp,     -1.0_rp,  &
     1.0_rp,      2.0_rp,      1.0_rp,     -1.0_rp,     -1.0_rp,     -1.0_rp,     -1.0_rp,     -1.0_rp,     -1.0_rp,     -1.0_rp,     -1.0_rp,  &
     1.0_rp,      3.0_rp,      3.0_rp,      1.0_rp,     -1.0_rp,     -1.0_rp,     -1.0_rp,     -1.0_rp,     -1.0_rp,     -1.0_rp,     -1.0_rp,  &
     1.0_rp,      4.0_rp,      6.0_rp,      4.0_rp,      1.0_rp,     -1.0_rp,     -1.0_rp,     -1.0_rp,     -1.0_rp,     -1.0_rp,     -1.0_rp,  &
     1.0_rp,      5.0_rp,     10.0_rp,     10.0_rp,      5.0_rp,      1.0_rp,     -1.0_rp,     -1.0_rp,     -1.0_rp,     -1.0_rp,     -1.0_rp,  &
     1.0_rp,      6.0_rp,     15.0_rp,     20.0_rp,     15.0_rp,      6.0_rp,      1.0_rp,     -1.0_rp,     -1.0_rp,     -1.0_rp,     -1.0_rp,  &
     1.0_rp,      7.0_rp,     21.0_rp,     35.0_rp,     35.0_rp,     21.0_rp,      7.0_rp,      1.0_rp,     -1.0_rp,     -1.0_rp,     -1.0_rp,  &
     1.0_rp,      8.0_rp,     28.0_rp,     56.0_rp,     70.0_rp,     56.0_rp,     28.0_rp,      8.0_rp,      1.0_rp,     -1.0_rp,     -1.0_rp,  &
     1.0_rp,      9.0_rp,     36.0_rp,     84.0_rp,    126.0_rp,    126.0_rp,     84.0_rp,     36.0_rp,      9.0_rp,      1.0_rp,     -1.0_rp,  &
     1.0_rp,     10.0_rp,     45.0_rp,    120.0_rp,    210.0_rp,    252.0_rp,    210.0_rp,    120.0_rp,     45.0_rp,     10.0_rp,      1.0_rp,  &
     1.0_rp,     11.0_rp,     55.0_rp,    165.0_rp,    330.0_rp,    462.0_rp,    462.0_rp,    330.0_rp,    165.0_rp,     55.0_rp,     11.0_rp,  &
     1.0_rp,     12.0_rp,     66.0_rp,    220.0_rp,    495.0_rp,    792.0_rp,    924.0_rp,    792.0_rp,    495.0_rp,    220.0_rp,     66.0_rp,  &
     1.0_rp,     13.0_rp,     78.0_rp,    286.0_rp,    715.0_rp,   1287.0_rp,   1716.0_rp,   1716.0_rp,   1287.0_rp,    715.0_rp,    286.0_rp,  &
     1.0_rp,     14.0_rp,     91.0_rp,    364.0_rp,   1001.0_rp,   2002.0_rp,   3003.0_rp,   3432.0_rp,   3003.0_rp,   2002.0_rp,   1001.0_rp,  &
     1.0_rp,     15.0_rp,    105.0_rp,    455.0_rp,   1365.0_rp,   3003.0_rp,   5005.0_rp,   6435.0_rp,   6435.0_rp,   5005.0_rp,   3003.0_rp,  &
     1.0_rp,     16.0_rp,    120.0_rp,    560.0_rp,   1820.0_rp,   4368.0_rp,   8008.0_rp,  11440.0_rp,  12870.0_rp,  11440.0_rp,   8008.0_rp,  &
     1.0_rp,     17.0_rp,    136.0_rp,    680.0_rp,   2380.0_rp,   6188.0_rp,  12376.0_rp,  19448.0_rp,  24310.0_rp,  24310.0_rp,  19448.0_rp,  &
     1.0_rp,     18.0_rp,    153.0_rp,    816.0_rp,   3060.0_rp,   8568.0_rp,  18564.0_rp,  31824.0_rp,  43758.0_rp,  48620.0_rp,  43758.0_rp,  &
     1.0_rp,     19.0_rp,    171.0_rp,    969.0_rp,   3876.0_rp,  11628.0_rp,  27132.0_rp,  50388.0_rp,  75582.0_rp,  92378.0_rp,  92378.0_rp,  &
     1.0_rp,     20.0_rp,    190.0_rp,   1140.0_rp,   4845.0_rp,  15504.0_rp,  38760.0_rp,  77520.0_rp, 125970.0_rp, 167960.0_rp, 184756.0_rp   &
], [21,11], order = [2,1])

integer n, k

character(*), parameter :: r_name = 'n_choose_k'

! Error check

if (k < 0 .or. n < k) then
  call out_io (s_error$, r_name, 'N_CHOOSE_K CALLED WITH K < 0 or N < K: ' // int_str(n) // ', ' //int_str(k))
  if (global_com%exit_on_error) call err_exit
  nck = -1
endif

! Use Sterling's formula if n is very large

if (n > ubound(f, 1)) then
  if (n > 170) then
    call out_io (s_error$, r_name, 'N_CHOOSE_K(N) CALLED WITH N TOO LARGE (> 170)!')
    if (global_com%exit_on_error) call err_exit
    nck = -1
    return
  endif

  nck = factorial(n) / (factorial(k) * factorial(n-k))
elseif (2*k < n) then
  nck = f(n,k)
else
  nck = f(n,n-k)
endif

end function
