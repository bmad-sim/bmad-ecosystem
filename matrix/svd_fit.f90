!+
! Subroutine svd_fit (A, b, cutoff, x_min)
!
! Routine to minimize Ax - b using SVD.
!
! Input:
!   A(:,:)      -- real(rp): Matrix size (n_data,n_var)
!   b(:)        -- real(rp): Vector size (n_data)
!   cutoff      -- real(rp): Sigular value cutoff. A reasonable value is 1e-5_rp.
!
! Output:
!   A(:,:)      -- real(rp): Garbage values.
!   x_min(:)    -- real(rp): X at minimum. Vector size (n_var).
!-

subroutine svd_fit (A, b, cutoff, x_min)

use f95_lapack
use precision_def

implicit none

real(rp) A(:,:), b(:), cutoff, x_min(:)
integer :: n_dat, n_var
real(rp) w_cut, w(size(A,2)), v(size(A,2), size(A,2))

! Lapack95 routine. Note that this returns V^Transpose, not V, so we need an extra step

n_dat = size(b);  n_var = size(A,2)

call dgesdd_f95(A, w(1:min(n_dat,n_var)), VT=v, JOB = 'U') 
v = transpose(v)

! Set any extra components of singluar value list w to zero

if (n_dat < n_var) w(n_dat+1:) = 0

! And compute x_min

w_cut = cutoff * maxval(w)

where (w > w_cut)
  x_min = matmul(b,A)/w
elsewhere
  x_min = 0.0
end where

x_min = matmul(v, x_min)

end subroutine
