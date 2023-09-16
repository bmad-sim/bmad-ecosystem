!+
! Subroutine spin_mat_to_eigen (orb_mat, spin_map, orb_eval, orb_evec, n0, spin_evec, error)
!
! Routine to return the spin/orbit eigen values and eigen vectors given the 1-turn spin/orbit map.
!
! Caution: n0 as calculated by this routine may be antiparallel n0 calculated elsewhere.
! To make things consistent, spin_evec needs to be flipped if n0 is flipped.
!
! Input:
!   orb_mat(6,6)      -- real(rp): Orbital matrix.
!   spin_map(0:3,0:6) -- real(rp): Quaternion 0th & 1st order map.
!
! Output:
!   orb_eval(6)       -- complex(rp): Eigenvalues.
!   orb_evec(6,6)     -- complex(rp): Orbital eigenvectors. orb_evec(j,:) is the j^th vector.
!   n0(3)             -- real(rp): n_0 invariant spin
!   spin_evec(6,3)    -- complex(rp): Spin eigenvectors. spin_evec(j,:) is the j^th vector.
!   error             -- logical: Set true if there is an error. False otherwise.
!-

subroutine spin_mat_to_eigen (orb_mat, spin_map, orb_eval, orb_evec, n0, spin_evec, error)

use sim_utils
use f95_lapack

implicit none

real(rp) orb_mat(6,6), spin_map(0:3,0:6), n0(3)
real(rp) w0(3,3), mat1(6,6), x_axis(3), y_axis(3)
complex(rp) orb_eval(6), orb_evec(6,6), spin_evec(6,3)
complex(rp) n0q(0:3), q0(0:3), qv(0:3), qq(0:3), dd(2,2), gv(2,1)

integer k, ipiv2(2), pinfo
logical error

! With RF off the eigen analysis is singular. 
! In this case, modify the matrix to remove the singularity.

mat1 = orb_mat
if (mat1(6,5) == 0) then ! Eigen anal is singular without RF so put in a fix.
  mat1(5,1:4) = 0
  mat1(5,6) = 0
endif

call mat_eigen(mat1, orb_eval, orb_evec, error)
if (error) then
  orb_eval = 0;  orb_evec = 0
  n0 = 0;  spin_evec = 0
  return
endif

w0 = quat_to_w_mat(spin_map(:,0))
n0 = spin_map(1:3,0) / norm2(spin_map(1:3,0))
n0q = [0.0_rp, n0]
q0 = spin_map(0:3,0)

! The equation for spin_evec is degenerate in the n0 direction so round off error is a problem.
! To get around this, reduce the problem to an equation in the plane (x, y) transverse to n0.

k = minloc(abs(n0), 1)
x_axis = n0 * n0(k)
x_axis(k) = x_axis(k) - 1
x_axis = -x_axis / norm2(x_axis)
y_axis = cross_product(n0, x_axis)

do k = 1, 6
  dd(1,:) = -[sum(matmul(w0, x_axis) * x_axis), sum(matmul(w0, y_axis) * x_axis)]
  dd(2,:) = -[sum(matmul(w0, x_axis) * y_axis), sum(matmul(w0, y_axis) * y_axis)]
  dd(1,1) = orb_eval(k) + dd(1,1)
  dd(2,2) = orb_eval(k) + dd(2,2)

  qv = matmul(spin_map(0:3,1:6), orb_evec(k,:))
  qq = quat_mul(quat_mul(qv, n0q), quat_conj(q0)) + quat_mul(quat_mul(q0, n0q), quat_conj(qv))
  gv(:,1) = [sum(qq(1:3) * x_axis), sum(qq(1:3) * y_axis)]

  call zgesv_f95 (dd, gv, ipiv2, pinfo)
  spin_evec(k,:) = gv(1,1) * x_axis + gv(2,1) * y_axis
enddo

end subroutine 
