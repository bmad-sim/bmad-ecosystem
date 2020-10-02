!+
! Function dn_dpz_from_qmap (orb_mat, q_map) result (dn_dpz)
!
! Rotine to calculate dn_dpz from th 1-turn 6x6 orbital matrix and the linear quaternian map
! This uses the standard SLIM coordinates for the spin motion.
!
! Note: The 1-turn spin quaternion q(r) at phase space point r(1:6) is
!   q(r) = q_map(0:3,0) + q_map(0:3,1:6) * r(1:6)
!
! Input:
!   orb_mat(6,6)    -- real(rp): 1-turn orbital matrix.
!   q_map(0:3,0:6)  -- real(rp): 1-turn spin linear quaternion map.
!
! Output:
!   dn_dpz(3)       -- real(rp): dn_dpz.
!-

function dn_dpz_from_qmap (orb_mat, q_map) result (dn_dpz)

use sim_utils
use eigen_mod
use f95_lapack

implicit none

real(rp) orb_mat(6,6), q_map(0:3,0:6), dn_dpz(3)
real(rp) w0(3,3)
complex(rp) n0(0:3), q0(0:3), n_eigen(6,3), vv(6,6), aa(6,1)
complex(rp) eval(6), evec(6,6), qv(0:3), qq(0:3), dd(3,3), gv(3,1)
integer k, pinfo, ipiv2(3), ipiv6(6)
logical err

!

call mat_eigen(orb_mat, eval, evec, err)
w0 = quat_to_w_mat(q_map(:,0))
n0 = [0.0_rp, q_map(1:3,0)] / norm2(q_map(1:3,0))
q0 = q_map(0:3,0)

do k = 1, 6
  call cmplx_mat_make_unit(dd)
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

end function dn_dpz_from_qmap
