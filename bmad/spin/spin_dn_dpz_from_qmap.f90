!+
! Function spin_dn_dpz_from_qmap (orb_mat, q_map, dn_dpz_partial, dn_dpz_partial2, error, n0) result (dn_dpz)
!
! Routine to calculate dn_dpz from th 1-turn 6x6 orbital matrix and the linear quaternian map
! Also see spin_dn_dpz_from_mat8 for the SLIM formalism equivalent.
!
! Note: The 1-turn spin quaternion q(r) at phase space point r(1:6) is
!   q(r) = q_map(0:3,0) + q_map(0:3,1:6) * r(1:6)
!
! Input:
!   orb_mat(6,6)          -- real(rp): 1-turn orbital matrix.
!   q_map(0:3,0:6)        -- real(rp): 1-turn spin linear quaternion map.
!   dn_dpz_partial(3,3)   -- real(rp) dn_dpz_partial(i,:) is dn_dpz with only one osccilation 
!                               mode "excited". So dn_dpz_partial(1,:) represents a-mode excitation, etc.
!   dn_dpz_partial2(3,3)  -- real(rp) dn_dpz_partial(i,:) is dn_dpz with only two osccilation 
!                               modes "excited". So dn_dpz_partial(1,:) represents b-mode and c-mode
!                               excitation without the a-mode, etc.
!   n0(3)                 -- real(rp) n0 vector. This needs to be supplied if there is the
!                               possibility that n0 is antiparallel to q_map(0:3,0).
!
! Output:
!   dn_dpz(3)       -- real(rp): dn_dpz.
!   error           -- logical: Set True if there is an error. False otherwise.
!-

function spin_dn_dpz_from_qmap (orb_mat, q_map, dn_dpz_partial, dn_dpz_partial2, error, n0) result (dn_dpz)

use bmad_routine_interface, dummy => spin_dn_dpz_from_qmap
use f95_lapack

implicit none

real(rp), optional :: n0(3)
real(rp) :: dn_dpz_partial(3,3), dn_dpz_partial2(3,3)
real(rp) orb_mat(6,6), q_map(0:3,0:6), dn_dpz(3), n0_calc(3)
complex(rp) n_eigen(6,3), vv(6,6), aa(6,1), aa4(4), n4(4,3)
complex(rp) eval(6), evec(6,6)
integer k, kk, pinfo, ipiv6(6), sgn
logical error

!

call spin_mat_to_eigen (orb_mat, q_map, eval, evec, n0_calc, n_eigen, error)
if (error) then
  dn_dpz = 0
  return
endif

sgn = 1
if (present(n0)) sgn = sign_of(dot_product(n0, n0_calc))

vv = transpose(evec)
aa(:,1) = [0, 0, 0, 0, 0, 1]
call zgesv_f95(vv, aa, ipiv6, pinfo)
dn_dpz = sgn * real(matmul(aa(:,1), n_eigen), rp)

!

do k = 1, 3
  kk = 2*k - 1
  dn_dpz_partial(k,:) = sgn * real(matmul(aa(kk:kk+1,1), n_eigen(kk:kk+1,:)), rp)
enddo

aa4 = aa(3:6,1);  n4 = n_eigen(3:6,:)
dn_dpz_partial2(1,:) = sgn * real(matmul(aa4, n4), rp)

aa4= [aa(1:2,1), aa(5:6,1)];  n4(1:2,:) = n_eigen(1:2,:);  n4(3:4,:) = n_eigen(5:6,:)
dn_dpz_partial2(2,:) = sgn * real(matmul(aa4, n4), rp)

aa4 = aa(1:4,1);  n4 = n_eigen(1:4,:)
dn_dpz_partial2(3,:) = sgn * real(matmul(aa4, n4), rp)

end function spin_dn_dpz_from_qmap
