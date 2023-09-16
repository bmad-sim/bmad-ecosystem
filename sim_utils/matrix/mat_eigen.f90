!+
! Subroutine mat_eigen (mat, eigen_val, eigen_vec, error, print_err)
!
! Routine for determining the eigenvectors and eigenvalues of a matrix.
! For complex conjugate pairs, the eigenvectors are normalized so that (v_j^*) * S * (v_j) = I for j odd.
!
! When the eigenvalues are complex conjugate pairs, the eigenvectors and eigenvalues
! are grouped so that the conjugate pairs are in slots (1,2), (3,4), etc.
!
! Also: For complex conjugate pairs (j, j+1), j = 1, 3, or 5, the odd numbered 
! eigenvector/eigenvalues will be such that the product (v{j}^*) (v{j+1}) has 
! positive imaginary part where v{j} is the j^th component of the eigenvector.
! This ensures that the odd numbered eigenvector/eigenvalues are associated
! with the tune and the even numbered ones have the negative of the tune.
! See the Bmad manual for more details.
!
! Also: the eigenvectors will be sorted so that:
!   Eigenvectors (1,2) will have the largest "horizontal" components corresponding to vec(1:2), 
!   Eigenvectors (3,4) will have the largest "vertical" components corresponding to vec(3:4), etc.
!
! Input:
!   mat(n,n)  -- Real(rp): Matrix. n must be even.
!   print_err -- Logical, optional: If present and False then suppress 
!                  "no eigen-system found" messages.
!
! Output:
!   eigen_val(n)   -- Complex(rp): Eigenvalues.
!   eigen_vec(n,n) -- Complex(rp): Eigenvectors are rows: eigen_vec(i,:)
!   error          -- Logical: Set True on error. False otherwise.
!-

subroutine mat_eigen (mat, eigen_val, eigen_vec, error, print_err)

use sim_utils, dummy => mat_eigen
use f95_lapack

implicit none

real(rp) mat(:,:), mat2(size(mat,1), size(mat,2))
real(rp) :: eval_r(size(mat,1)), eval_i(size(mat,1)), vec(size(mat,1), size(mat,1)), tmp(size(mat,1))
real(rp) :: amp(size(mat,1)/2, size(mat,1)/2), amp2(size(mat,1)/2, size(mat,1)/2), dmax(size(mat,1)/2)
real(rp) fnorm, fs

integer :: sort(size(mat,1)/2)

complex(rp) eigen_val(:), eigen_vec(:,:), sgn

integer i, j, k, ii, jj, kk, m, m2, n, nn, ier

logical, optional :: print_err
logical error, err

character(*), parameter :: r_name = 'mat_eigen'

!

n = size(mat,1)
nn = n / 2

mat2 = mat
call la_geev(mat2, eval_r, eval_i, vr = vec, info = ier)
error = (ier /= 0)
if (error) then
  if (logic_option(.true., print_err)) call out_io (s_fatal$, r_name, "la_geev returned error: \i0\ ", ier)
  eigen_val = 0.0d0
  eigen_vec = 0.0d0
  return
endif

! la_geev orders stable eigenvectors to come in pairs but not unstable ones. So reorder.

do ii = 1, nn-1
  j = 2*ii-1
  if (eval_i(j) == 0 .and. eval_i(j+1) /= 0) then
    eval_i(j:) = [eval_i(j+1:), eval_i(j)]
    eval_r(j:) = [eval_r(j+1:), eval_r(j)]
    tmp = vec(:,j)
    vec(:,j:n-1) = vec(:,j+1:n)
    vec(:,n) = tmp
  endif
enddo

! Order eigen with horizontal first etc. 

do ii = 1, nn
  i = 2 * ii
  do kk = 1, nn
    k = 2 * kk
    amp(kk,ii) = sum(vec(k-1:k, i-1:i)**2)
  enddo
  amp(:,ii) = amp(:,ii) / sum(amp(:,ii))
enddo

do ii = 1, nn
  amp2 = amp

  do j = 1, nn
    m = maxloc(amp2(j,:), 1)
    amp2(j, m) = -1
    m2 = maxloc(amp2(j,:), 1)
    dmax(j) = amp(j,m) - amp(j,m2)
  enddo

  k = maxloc(dmax, 1)
  m = maxloc(amp(k,:), 1)
  sort(m) = k
  amp(k,:) = 0
  amp(:,m) = 0
enddo

! With degenerate eigenvalues (happens when eigenvalues = 1), la_geev does not return orthogonal eigenvectors.
! EG: vec5 = [0, 0, 0, 0, 1, 0], vec6 = [0, 0, 0, 0, 1, 1e-18] (yes this is an actual result).
! This can be confusing so orthoganlize.

do ii = 1, nn
  i = 2 * ii
  j = 2 * sort(ii)

  if (eval_i(i-1) == 0) then  ! Unstable mode

    fs = sum(vec(:,i-1)**2)
    if (eval_r(i-1) == eval_r(i) .and. abs(fs) > 1e-100_rp) then ! orthogonalize
      fnorm = dot_product(vec(:,i-1), vec(:,i)) / fs
      vec(:,i) = vec(:,i) - fnorm * vec(:,i-1)
    endif

    fnorm = sqrt(sum(vec(:, i-1)**2))
    if (abs(fnorm) > 1e-100_rp) vec(:, i-1) = vec(:, i-1) / fnorm

    fnorm = sqrt(sum(vec(:, i)**2))
    if (abs(fnorm) > 1e-100_rp) vec(:, i) = vec(:, i) / fnorm

    eigen_val(j-1)    = eval_r(i-1)
    eigen_val(j)      = eval_r(i)
    eigen_vec(j-1, :) = vec(:, i-1)
    eigen_vec(j, :)   = vec(:, i)

  else    ! Stable mode with complex conjugate pairs.
    eigen_val(j-1)    = cmplx(eval_r(i), eval_i(i))
    eigen_val(j)      = cmplx(eval_r(i-1), eval_i(i-1))
    eigen_vec(j-1, :) = cmplx(vec(:, i-1), -vec(:, i))
    eigen_vec(j, :)   = cmplx(vec(:, i-1),  vec(:, i))

    k = j - 1
    fnorm = 0
    do jj = 1, n, 2
      fnorm = fnorm + 2 * aimag(conjg(eigen_vec(k, jj)) * eigen_vec(k,jj+1))
    enddo

    if (fnorm < 0) then  ! flip
      eigen_vec(k:k+1, :) = eigen_vec(k+1:k:-1, :)
      eigen_val(k:k+1)    = eigen_val(k+1:k:-1)
      fnorm = -fnorm
    endif

    if (abs(eigen_vec(k,k)) == 0) then
      sgn = 1.0_rp / sqrt(fnorm)
    else
      sgn = abs(eigen_vec(k,k)) / (sqrt(fnorm) * eigen_vec(k,k))
    endif

    eigen_vec(k,:)   = sgn        * eigen_vec(k,:)
    eigen_vec(k+1,:) = conjg(sgn) * eigen_vec(k+1,:)
  endif
enddo

error = .false.

end subroutine mat_eigen
