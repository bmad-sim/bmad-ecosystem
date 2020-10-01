!+
! Subroutine dn_dpz_from_mat8 (mat_1turn, dn_dpz)
!
! Rotine to calculate dn_dpz from a 1-turn 8x8 spin-orbital matrix.
!
! This uses the SLIM spin (l,n,m) coordinates as discussed in the Bmad manual.
! Note: Since dn_dpz is perpendicular to the n-axis, it is always the case that dn_dpz = 0.
!
! Input:
!   mat_1turn(8,8) -- real(rp): Spin-orbital matrix.
!
! Output:
!   dn_dpz(3)      -- real(rp): dn_dpz (l,n,m) coordinates.
!
!-

subroutine dn_dpz_from_mat8 (mat_1turn, dn_dpz)

use sim_utils
use eigen_mod
use f95_lapack

implicit none

real(rp) mat_1turn(8,8), dn_dpz(3)
complex(rp) eval(6), evec(6,6), dd(2,2), gv(2,1), w_vec(6,2), vv(6,6), aa(6,1)
integer k, pinfo, ipiv2(2), ipiv6(6)
logical err

!

call mat_eigen(mat_1turn(1:6,1:6), eval, evec, err)

do k = 1, 6
  dd = mat_1turn(7:8,7:8)
  dd(1,1) = dd(1,1) - eval(k)
  dd(2,2) = dd(2,2) - eval(k)
  gv(:,1) = -matmul(mat_1turn(7:8,1:6), evec(k,1:6))
  call zgesv_f95 (dd, gv, ipiv2, pinfo)
  w_vec(k,:) = gv(:,1)
enddo

vv = transpose(evec)
aa(:,1) = [0, 0, 0, 0, 0, 1]
call zgesv_f95(vv, aa, ipiv6, pinfo)
dn_dpz = [real(sum(w_vec(:,1)* aa(:,1)), rp), 0.0_rp, real(sum(w_vec(:,2)* aa(:,1)), rp)]

end subroutine dn_dpz_from_mat8
