!+
! Subroutine svd_fit (A, b, cutoff, x_min, chisq, w_vec, v_mat)
!
! Routine to minimize Ax - b using SVD.
!
! Input:
!   A(:,:)      -- real(rp): Matrix of size (n_data,n_var)
!   b(:)        -- real(rp): Vector of size n_data.
!   cutoff      -- real(rp): Sigular value cutoff. A reasonable value is 1e-5_rp.
!
! Output:
!   A(:,:)      -- real(rp): SVD U matrix.
!   x_min(:)    -- real(rp): X at minimum. Vector of size n_var.
!   chisq       -- real(rp), optional: Chi^2 figure of merit.
!   w_vec(:)    -- real(rp), optional: SVD W vector with cutoff applied. Of size n_var.
!   v_mat(:,:)  -- real(rp), optional: SVD V matrix (not transpose) of size (n_var,n_var).
!-

subroutine svd_fit (A, b, cutoff, x_min, chisq, w_vec, v_mat)

use f95_lapack, only: dgesdd_f95
use precision_def

implicit none

real(rp) A(:,:), b(:), cutoff, x_min(:)
real(rp), optional :: chisq, w_vec(:), v_mat(:,:)
integer :: n_dat, n_var
real(rp) w_cut, w(size(A,2)), v(size(A,2), size(A,2)), A2(size(A,1), size(A,2))

! Lapack95 routine. Note that this returns V^Transpose, not V, so we need an extra step

w = 0
A2 = A
n_dat = size(b);  n_var = size(A,2)

call dgesdd_f95(A, w(1:min(n_dat,n_var)), vt=v, job = 'U') 
v = transpose(v)

! Set any extra components of singluar value list w to zero

if (n_dat < n_var) w(n_dat+1:) = 0

! And compute x_min

w_cut = cutoff * maxval(w)

where (w > w_cut)
  x_min = matmul(b,A)/w
elsewhere
  x_min = 0.0
  w = 0
end where

x_min = matmul(v, x_min)

if (present(chisq)) chisq = norm2(matmul(A2,x_min)-b)**2
if (present(w_vec)) w_vec = w
if (present(v_mat)) v_mat = v

end subroutine
