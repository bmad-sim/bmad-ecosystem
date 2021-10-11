!+
! Subroutine spin_mat8_resonance_strengths (evec, mat8, xi_mat8)
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
!   xi_mat8(2)        -- real(rp): Resonance strengths for sum and difference resonances.
!-

subroutine spin_mat8_resonance_strengths (orb_evec, mat8, xi_mat8)

use sim_utils

implicit none

real(rp) mat8(8,8), xi_mat8(2), nn0(3)
complex(rp) orb_evec(6)

!

xi_mat8(1) = abs(sum(orb_evec        * (mat8(7,1:6) + mat8(8,1:6) * i_imag))) / twopi
xi_mat8(2) = abs(sum(conjg(orb_evec) * (mat8(7,1:6) + mat8(8,1:6) * i_imag))) / twopi

end subroutine
