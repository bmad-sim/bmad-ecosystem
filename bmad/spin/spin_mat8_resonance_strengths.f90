!+
! Subroutine spin_mat8_resonance_strengths (evec, mat8, xi_sum, xi_diff)
!
! Routine to calculate for linear spin/orbit resonances the resonance strength from the 
! mat8 spin/orbital matrix and a particular eigen mode.
!
! Note: This routine will not be accurate unless the machine is at a sum or difference resonance.
! Also see: spin_quat_resonance_strengths.
!
! Input:
!   orb_evec(6)       -- complex(rp): Orbital eigenvector.
!   mat8(6,6)         -- real(rp): Spin/orbital matrix.
!
! Output:
!   xi_sum            -- real(rp): Sum resonance strength.
!   xi_diff           -- real(rp): Difference resonance strength.
!-

subroutine spin_mat8_resonance_strengths (orb_evec, mat8, xi_sum, xi_diff)

use sim_utils

implicit none

real(rp) mat8(8,8), xi_sum, xi_diff, nn0(3)
complex(rp) orb_evec(6)

!

xi_sum  = abs(sum(orb_evec        * (mat8(7,1:6) + mat8(8,1:6) * i_imag))) / twopi
xi_diff = abs(sum(conjg(orb_evec) * (mat8(7,1:6) + mat8(8,1:6) * i_imag))) / twopi

end subroutine
