!+
! Regression test for ran_gauss_scalar and ran_gauss_vector.
!
! Tests:
!   1. Deterministic reproducibility (fixed seed gives same sequence)
!   2. Mean ~ 0 (z-test)
!   3. Variance ~ 1 (chi-squared test via ratio)
!   4. Skewness ~ 0
!   5. Kurtosis ~ 3 (excess kurtosis ~ 0)
!   6. Anderson-Darling normality test
!   7. Sigma-cut truncation works correctly
!   8. Scalar and vector paths produce identical sequences
!-

program ran_gauss_test

use random_mod

implicit none

integer, parameter :: n_samples = 100000
integer, parameter :: n_vec = 6

real(rp) :: scalar_vals(n_samples)
real(rp) :: vector_vals(n_samples)
real(rp) :: vec_buf(n_vec)
real(rp) :: mean, var, skew, kurt, ad_stat
real(rp) :: sig_cut_val, max_abs
real(rp) :: sorted(n_samples)
real(rp) :: z, phi_z, s_term
real(rp) :: diff_max

integer :: i, j, idx

! Output file
open (1, file = 'output.now')

! ============================================================
! Test 1: Deterministic reproducibility - scalar path
! A fixed seed must produce the exact same sequence.
! ============================================================

call ran_seed_put(42)
call ran_gauss(scalar_vals(1))
call ran_gauss(scalar_vals(2))
call ran_gauss(scalar_vals(3))
call ran_gauss(scalar_vals(4))

write (1, '(a, 4f18.12)') '"Repro-scalar" ABS 1e-12', scalar_vals(1:4)

! ============================================================
! Test 2: Deterministic reproducibility - vector path
! ============================================================

call ran_seed_put(42)
call ran_gauss(vec_buf)
write (1, '(a, 6f18.12)') '"Repro-vector" ABS 1e-12', vec_buf

! ============================================================
! Test 3: Scalar and vector paths produce identical sequences
! Same seed, scalar one-at-a-time vs vector batch must match.
! ============================================================

call ran_seed_put(99)
do i = 1, 12
  call ran_gauss(scalar_vals(i))
enddo

call ran_seed_put(99)
call ran_gauss(vector_vals(1:6))
call ran_gauss(vector_vals(7:12))

diff_max = 0
do i = 1, 12
  diff_max = max(diff_max, abs(scalar_vals(i) - vector_vals(i)))
enddo

write (1, '(a, es16.8)') '"Scalar-vs-vector-diff" ABS 1e-12', diff_max

! ============================================================
! Generate large sample for statistical tests using vector path
! ============================================================

call ran_seed_put(12345)
idx = 1
do while (idx <= n_samples)
  call ran_gauss(vec_buf)
  do j = 1, n_vec
    if (idx <= n_samples) then
      vector_vals(idx) = vec_buf(j)
      idx = idx + 1
    endif
  enddo
enddo

! ============================================================
! Test 4: Mean (z-test, expect |mean| < 4*sigma/sqrt(n))
! For N(0,1), sigma_mean = 1/sqrt(n) ~ 0.00316
! Tolerance: |mean| < 0.02 (~ 6 sigma, extremely conservative)
! ============================================================

mean = sum(vector_vals(1:n_samples)) / n_samples
write (1, '(a, f18.12)') '"Vec-mean" ABS 2e-2', mean

! ============================================================
! Test 5: Variance (expect ~ 1.0)
! Var(sample_var) ~ 2/n for normal, so sigma ~ 0.0045
! Tolerance: |var - 1| < 0.03 (~ 6 sigma)
! ============================================================

var = sum((vector_vals(1:n_samples) - mean)**2) / (n_samples - 1)
write (1, '(a, f18.12)') '"Vec-variance" ABS 3e-2', var - 1.0_rp

! ============================================================
! Test 6: Skewness (expect ~ 0)
! Var(skewness) ~ 6/n, sigma ~ 0.0077
! Tolerance: |skew| < 0.05
! ============================================================

skew = sum(((vector_vals(1:n_samples) - mean) / sqrt(var))**3) / n_samples
write (1, '(a, f18.12)') '"Vec-skewness" ABS 5e-2', skew

! ============================================================
! Test 7: Kurtosis (expect ~ 3, excess kurtosis ~ 0)
! Var(kurtosis) ~ 24/n, sigma ~ 0.0155
! Tolerance: |kurt - 3| < 0.10
! ============================================================

kurt = sum(((vector_vals(1:n_samples) - mean) / sqrt(var))**4) / n_samples
write (1, '(a, f18.12)') '"Vec-kurtosis" ABS 1e-1', kurt - 3.0_rp

! ============================================================
! Test 8: Anderson-Darling test for normality
! AD statistic for N(0,1): critical value at 1% is 1.0348
! We use a generous threshold of 2.0 so deterministic test won't
! flap, while still catching badly broken distributions.
! ============================================================

! Sort the sample
sorted = vector_vals(1:n_samples)
call shell_sort(sorted)

ad_stat = 0.0_rp
do i = 1, n_samples
  z = sorted(i)
  ! Standard normal CDF via erfc: Phi(z) = 0.5 * erfc(-z/sqrt(2))
  phi_z = 0.5_rp * erfc(-z / sqrt(2.0_rp))
  ! Clamp to avoid log(0)
  phi_z = max(phi_z, 1.0e-15_rp)
  phi_z = min(phi_z, 1.0_rp - 1.0e-15_rp)

  s_term = (2.0_rp * i - 1.0_rp) * log(phi_z) + &
           (2.0_rp * (n_samples - i) + 1.0_rp) * log(1.0_rp - phi_z)
  ad_stat = ad_stat + s_term
enddo
ad_stat = -real(n_samples, rp) - ad_stat / real(n_samples, rp)

write (1, '(a, f18.12)') '"Vec-AD-stat" ABS 2.0', ad_stat

! ============================================================
! Generate large sample for statistical tests using scalar path
! ============================================================

call ran_seed_put(54321)
do i = 1, n_samples
  call ran_gauss(scalar_vals(i))
enddo

mean = sum(scalar_vals(1:n_samples)) / n_samples
write (1, '(a, f18.12)') '"Scalar-mean" ABS 2e-2', mean

var = sum((scalar_vals(1:n_samples) - mean)**2) / (n_samples - 1)
write (1, '(a, f18.12)') '"Scalar-variance" ABS 3e-2', var - 1.0_rp

skew = sum(((scalar_vals(1:n_samples) - mean) / sqrt(var))**3) / n_samples
write (1, '(a, f18.12)') '"Scalar-skewness" ABS 5e-2', skew

kurt = sum(((scalar_vals(1:n_samples) - mean) / sqrt(var))**4) / n_samples
write (1, '(a, f18.12)') '"Scalar-kurtosis" ABS 1e-1', kurt - 3.0_rp

! Sort and AD test for scalar path
sorted = scalar_vals(1:n_samples)
call shell_sort(sorted)

ad_stat = 0.0_rp
do i = 1, n_samples
  z = sorted(i)
  phi_z = 0.5_rp * erfc(-z / sqrt(2.0_rp))
  phi_z = max(phi_z, 1.0e-15_rp)
  phi_z = min(phi_z, 1.0_rp - 1.0e-15_rp)

  s_term = (2.0_rp * i - 1.0_rp) * log(phi_z) + &
           (2.0_rp * (n_samples - i) + 1.0_rp) * log(1.0_rp - phi_z)
  ad_stat = ad_stat + s_term
enddo
ad_stat = -real(n_samples, rp) - ad_stat / real(n_samples, rp)

write (1, '(a, f18.12)') '"Scalar-AD-stat" ABS 2.0', ad_stat

! ============================================================
! Test 9: Sigma-cut truncation (vector path)
! All samples must be within [-sigma_cut, sigma_cut]
! ============================================================

sig_cut_val = 3.0_rp

call ran_seed_put(77777)
max_abs = 0.0_rp
do i = 1, n_samples / n_vec
  call ran_gauss(vec_buf, sigma_cut = sig_cut_val)
  do j = 1, n_vec
    max_abs = max(max_abs, abs(vec_buf(j)))
  enddo
enddo

! max_abs must be < sig_cut_val; write the excess (should be 0 or negative)
write (1, '(a, f18.12)') '"Vec-sigcut-ok" ABS 1e-12', max(0.0_rp, max_abs - sig_cut_val)

! ============================================================
! Test 10: Sigma-cut truncation (scalar path)
! ============================================================

call ran_seed_put(88888)
max_abs = 0.0_rp
do i = 1, n_samples
  call ran_gauss(scalar_vals(i), sigma_cut = sig_cut_val)
  max_abs = max(max_abs, abs(scalar_vals(i)))
enddo

write (1, '(a, f18.12)') '"Scalar-sigcut-ok" ABS 1e-12', max(0.0_rp, max_abs - sig_cut_val)

!

close (1)
print *, 'ran_gauss_test: All tests written to output.now'

contains

!+
! Simple shell sort for real(rp) array (ascending order).
!-

subroutine shell_sort(arr)
real(rp), intent(inout) :: arr(:)
real(rp) :: tmp
integer :: n_arr, gap, i_s, j_s

n_arr = size(arr)
gap = n_arr / 2
do while (gap > 0)
  do i_s = gap + 1, n_arr
    tmp = arr(i_s)
    j_s = i_s
    do while (j_s > gap .and. arr(j_s - gap) > tmp)
      arr(j_s) = arr(j_s - gap)
      j_s = j_s - gap
    enddo
    arr(j_s) = tmp
  enddo
  gap = gap / 2
enddo

end subroutine shell_sort

end program ran_gauss_test
