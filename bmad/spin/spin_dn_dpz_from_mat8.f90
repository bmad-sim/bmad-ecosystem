!+
! Function spin_dn_dpz_from_mat8 (mat_1turn, dn_dpz_partial, error) result (dn_dpz)
!
! Routine to calculate dn_dpz from a 1-turn 8x8 spin-orbital matrix.
! Also see spin_dn_dpz_from_qmap for the quaternion equivalent.
!
! This uses the SLIM spin (l,n,m) coordinates as discussed in the Bmad manual.
! Note: Since dn_dpz is perpendicular to the n-axis, it is always the case that dn_dpz = 0.
!
! Input:
!   mat_1turn(8,8)      -- real(rp): Spin-orbital matrix.
!   dn_dpz_partial(3,3) -- real(rp), optional: dn_dpz_partial(i,:) is dn_dpz with only one osccilation 
!                              mode "excited". So dn_dpz_partial(1,:) represents a-mode excitation, etc.
!
! Output:
!   dn_dpz(3)         -- real(rp): dn_dpz (l,n,m) coordinates.
!   error             -- logical: Set True if there is an error. False otherwise.
!-

function spin_dn_dpz_from_mat8 (mat_1turn, dn_dpz_partial, error) result (dn_dpz)

use bmad_routine_interface, dummy => spin_dn_dpz_from_mat8
use f95_lapack, only: zgesv_f95

implicit none

real(rp), optional :: dn_dpz_partial(3,3)
real(rp) mat_1turn(8,8), dn_dpz(3), mat1(6,6)
complex(rp) eval(6), evec(6,6), dd(2,2), gv(2,1), w_vec(6,2), vv(6,6), aa(6,1)
integer k, kk, pinfo, ipiv2(2), ipiv6(6)
logical error

! With RF off the eigen analysis is singular. 
! In this case, modify the matrix to remove the singularity.

mat1 = mat_1turn(1:6,1:6)
if (mat1(6,5) == 0) then ! Eigen anal is singular without RF so put in a fix.
  mat1(5,1:4) = 0
  mat1(5,6) = 0
endif

call mat_eigen(mat1(1:6,1:6), eval, evec, error)
if (error) then
  dn_dpz = 0
  return
endif

do k = 1, 6
  call cplx_mat_make_unit(dd)
  dd = eval(k) * dd - mat_1turn(7:8,7:8)
  gv(:,1) = matmul(mat_1turn(7:8,1:6), evec(k,1:6))
  call zgesv_f95 (dd, gv, ipiv2, pinfo)
  w_vec(k,:) = gv(:,1)
enddo

vv = transpose(evec)
aa(:,1) = [0, 0, 0, 0, 0, 1]
call zgesv_f95(vv, aa, ipiv6, pinfo)
dn_dpz = [real(sum(w_vec(:,1)* aa(:,1)), rp), 0.0_rp, real(sum(w_vec(:,2)* aa(:,1)), rp)]

!

if (present(dn_dpz_partial)) then
  do k = 1, 3
    kk = 2*k - 1
    dn_dpz = [real(sum(w_vec(kk:kk+1,1)* aa(kk:kk+1,1)), rp), 0.0_rp, real(sum(w_vec(kk:kk+1,2)* aa(kk:kk+1,1)), rp)]
  enddo
endif

end function spin_dn_dpz_from_mat8
