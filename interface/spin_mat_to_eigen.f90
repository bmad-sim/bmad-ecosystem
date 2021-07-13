!+
! Subroutine spin_mat_to_eigen (orb_mat, spin_map, eigen_val, orb_evec, n0, spin_evec)
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
!   eigen_val(6)      -- complex(rp): Eigenvalues.
!   orb_evec(6,6)     -- complex(rp): Orbital eigenvectors. orb_evec(j,:) is the j^th vector.
!   n0(3)             -- real(rp): n_0 invariant spin
!   spin_evec(6,3)    -- complex(rp): Spin eigenvectors. spin_evec(j,:) is the j^th vector.
!-

subroutine spin_mat_to_eigen (orb_mat, spin_map, eigen_val, orb_evec, n0, spin_evec)

use sim_utils
use eigen_mod
use f95_lapack

implicit none

real(rp) orb_mat(6,6), spin_map(0:3,0:6), n0(3)
real(rp) w0(3,3), mat1(6,6)
complex(rp) eigen_val(6), orb_evec(6,6), spin_evec(6,3)
complex(rp) n0q(0:3), q0(0:3), qv(0:3), qq(0:3), dd(3,3), gv(3,1)

integer k, ipiv2(3), pinfo
logical err

! With RF off the eigen analysis is singular. 
! In this case, modify the matrix to remove the singularity.

mat1 = orb_mat
if (mat1(6,5) == 0) then ! Eigen anal is singular without RF so put in a fix.
  mat1(5,1:4) = 0
  mat1(5,6) = 0
endif

call mat_eigen(mat1, eigen_val, orb_evec, err)

w0 = quat_to_w_mat(spin_map(:,0))
n0 = spin_map(1:3,0) / norm2(spin_map(1:3,0))
n0q = [0.0_rp, n0]
q0 = spin_map(0:3,0)

do k = 1, 6
  call cplx_mat_make_unit(dd)
  dd = eigen_val(k) * dd - w0
  qv = matmul(spin_map(0:3,1:6), orb_evec(k,:))
  qq = quat_mul(quat_mul(qv, n0q), quat_conj(q0)) + quat_mul(quat_mul(q0, n0q), quat_conj(qv))
  gv(:,1) = qq(1:3)
  call zgesv_f95 (dd, gv, ipiv2, pinfo)
  spin_evec(k,:) = gv(:,1)
enddo

end subroutine 
