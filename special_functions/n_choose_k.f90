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
     1_rp,     -1_rp,     -1_rp,     -1_rp,     -1_rp,     -1_rp,     -1_rp,     -1_rp,     -1_rp,     -1_rp,     -1_rp,  &
     1_rp,      1_rp,     -1_rp,     -1_rp,     -1_rp,     -1_rp,     -1_rp,     -1_rp,     -1_rp,     -1_rp,     -1_rp,  &
     1_rp,      2_rp,      1_rp,     -1_rp,     -1_rp,     -1_rp,     -1_rp,     -1_rp,     -1_rp,     -1_rp,     -1_rp,  &
     1_rp,      3_rp,      3_rp,      1_rp,     -1_rp,     -1_rp,     -1_rp,     -1_rp,     -1_rp,     -1_rp,     -1_rp,  &
     1_rp,      4_rp,      6_rp,      4_rp,      1_rp,     -1_rp,     -1_rp,     -1_rp,     -1_rp,     -1_rp,     -1_rp,  &
     1_rp,      5_rp,     10_rp,     10_rp,      5_rp,      1_rp,     -1_rp,     -1_rp,     -1_rp,     -1_rp,     -1_rp,  &
     1_rp,      6_rp,     15_rp,     20_rp,     15_rp,      6_rp,      1_rp,     -1_rp,     -1_rp,     -1_rp,     -1_rp,  &
     1_rp,      7_rp,     21_rp,     35_rp,     35_rp,     21_rp,      7_rp,      1_rp,     -1_rp,     -1_rp,     -1_rp,  &
     1_rp,      8_rp,     28_rp,     56_rp,     70_rp,     56_rp,     28_rp,      8_rp,      1_rp,     -1_rp,     -1_rp,  &
     1_rp,      9_rp,     36_rp,     84_rp,    126_rp,    126_rp,     84_rp,     36_rp,      9_rp,      1_rp,     -1_rp,  &
     1_rp,     10_rp,     45_rp,    120_rp,    210_rp,    252_rp,    210_rp,    120_rp,     45_rp,     10_rp,      1_rp,  &
     1_rp,     11_rp,     55_rp,    165_rp,    330_rp,    462_rp,    462_rp,    330_rp,    165_rp,     55_rp,     11_rp,  &
     1_rp,     12_rp,     66_rp,    220_rp,    495_rp,    792_rp,    924_rp,    792_rp,    495_rp,    220_rp,     66_rp,  &
     1_rp,     13_rp,     78_rp,    286_rp,    715_rp,   1287_rp,   1716_rp,   1716_rp,   1287_rp,    715_rp,    286_rp,  &
     1_rp,     14_rp,     91_rp,    364_rp,   1001_rp,   2002_rp,   3003_rp,   3432_rp,   3003_rp,   2002_rp,   1001_rp,  &
     1_rp,     15_rp,    105_rp,    455_rp,   1365_rp,   3003_rp,   5005_rp,   6435_rp,   6435_rp,   5005_rp,   3003_rp,  &
     1_rp,     16_rp,    120_rp,    560_rp,   1820_rp,   4368_rp,   8008_rp,  11440_rp,  12870_rp,  11440_rp,   8008_rp,  &
     1_rp,     17_rp,    136_rp,    680_rp,   2380_rp,   6188_rp,  12376_rp,  19448_rp,  24310_rp,  24310_rp,  19448_rp,  &
     1_rp,     18_rp,    153_rp,    816_rp,   3060_rp,   8568_rp,  18564_rp,  31824_rp,  43758_rp,  48620_rp,  43758_rp,  &
     1_rp,     19_rp,    171_rp,    969_rp,   3876_rp,  11628_rp,  27132_rp,  50388_rp,  75582_rp,  92378_rp,  92378_rp,  &
     1_rp,     20_rp,    190_rp,   1140_rp,   4845_rp,  15504_rp,  38760_rp,  77520_rp, 125970_rp, 167960_rp, 184756_rp   &
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
