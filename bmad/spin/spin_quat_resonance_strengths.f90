!+
! Subroutine spin_quat_resonance_strengths (orb_evec, spin_q, xi_sum, xi_diff)
!
! Routine to calculate for linear spin/orbit resonances the resonance strength from the 
! quaternion spin map and a particular eigen mode.
!
! Note: This routine will not be accurate unless the machine is at a sum or difference resonance.
! Also see: spin_mat8_resonance_strengths.
!
! Input:
!   orb_evec(6)       -- complex(rp): Orbital eigenvector.
!   spin_q(0:3,0:6)   -- real(rp): First order spin map.
!
! Output:
!   xi_sum            -- real(rp): Sum resonance strength.
!   xi_diff           -- real(rp): Difference resonance strength.
!-

subroutine spin_quat_resonance_strengths (orb_evec, spin_q, xi_sum, xi_diff)

use sim_utils

implicit none

real(rp) spin_q(0:3,0:6), xi_sum, xi_diff, nn0(3), anorm
complex(rp) orb_evec(6), qv(0:3), qv2(0:3), np(0:3), nm(0:3)

integer k

!

nn0 = spin_q(1:3,0) 
anorm = norm2(nn0)
if (anorm == 0) then
  xi_sum = 0
  xi_diff = 0
  return
endif

nn0 = sign_of(spin_q(0,0)) * nn0 / norm2(nn0)

np(0) = 0.5_rp
np(1:3) = i_imag * nn0 / 2
nm(0) = 0.5_rp
nm(1:3) = -i_imag * nn0 / 2

do k = 0, 3
  qv(k)  = sum(orb_evec(:) * spin_q(k,1:6))
  qv2(k) = sum(conjg(orb_evec(:)) * spin_q(k,1:6))
enddo

xi_sum  = sqrt(2.0) * norm2(abs(quat_mul(np, qv, nm))) / pi
xi_diff = sqrt(2.0) * norm2(abs(quat_mul(np, qv2, nm))) / pi

end subroutine
