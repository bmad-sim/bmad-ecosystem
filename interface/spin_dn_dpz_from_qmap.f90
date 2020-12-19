!+
! Function spin_dn_dpz_from_qmap (orb_mat, q_map, dn_dpz_partial) result (dn_dpz)
!
! Routine to calculate dn_dpz from th 1-turn 6x6 orbital matrix and the linear quaternian map
! Also see spin_dn_dpz_from_mat8 for the SLIM formalism equivalent.
!
! Note: The 1-turn spin quaternion q(r) at phase space point r(1:6) is
!   q(r) = q_map(0:3,0) + q_map(0:3,1:6) * r(1:6)
!
! Input:
!   orb_mat(6,6)        -- real(rp): 1-turn orbital matrix.
!   q_map(0:3,0:6)      -- real(rp): 1-turn spin linear quaternion map.
!   dn_dpz_partial(3,3) -- real(rp), optional: dn_dpz_partial(i,:) is dn_dpz with only one osccilation 
!                              mode "excited". So dn_dpz_partial(1,:) represents a-mode excitation, etc.
!
! Output:
!   dn_dpz(3)       -- real(rp): dn_dpz.
!-

function spin_dn_dpz_from_qmap (orb_mat, q_map, dn_dpz_partial) result (dn_dpz)

use sim_utils
use eigen_mod
use f95_lapack

implicit none

real(rp), optional :: dn_dpz_partial(3,3)
real(rp) orb_mat(6,6), q_map(0:3,0:6), dn_dpz(3)
real(rp) w0(3,3), mat1(6,6)
complex(rp) n0(0:3), q0(0:3), n_eigen(6,3), vv(6,6), aa(6,1)
complex(rp) eval(6), evec(6,6), qv(0:3), qq(0:3), dd(3,3), gv(3,1)
integer k, kk, pinfo, ipiv2(3), ipiv6(6)
logical err

! With RF off the eigen analysis is singular. 
! In this case, modify the matrix to remove the singularity.

mat1 = orb_mat
if (mat1(6,5) == 0) then ! Eigen anal is singular without RF so put in a fix.
  mat1(5,1:4) = 0
  mat1(5,6) = 0
endif

call mat_eigen(mat1, eval, evec, err)

w0 = quat_to_w_mat(q_map(:,0))
n0 = [0.0_rp, q_map(1:3,0)] / norm2(q_map(1:3,0))
q0 = q_map(0:3,0)

do k = 1, 6
  call cplx_mat_make_unit(dd)
  dd = eval(k) * dd - w0
  qv = matmul(q_map(0:3,1:6), evec(k,:))
  qq = quat_mul(quat_mul(qv, n0), quat_conj(q0)) + quat_mul(quat_mul(q0, n0), quat_conj(qv))
  gv(:,1) = qq(1:3)
  call zgesv_f95 (dd, gv, ipiv2, pinfo)
  n_eigen(k,:) = gv(:,1)
enddo

vv = transpose(evec)
aa(:,1) = [0, 0, 0, 0, 0, 1]
call zgesv_f95(vv, aa, ipiv6, pinfo)
dn_dpz = real(matmul(aa(:,1), n_eigen), rp)

!

if (present(dn_dpz_partial)) then
  do k = 1, 3
    kk = 2*k - 1
    dn_dpz_partial(k,:) = real(matmul(aa(kk:kk+1,1), n_eigen(kk:kk+1,:)), rp)
  enddo
endif

end function spin_dn_dpz_from_qmap
