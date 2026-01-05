!+
! Function dpc_given_dE(pc_old, mass, dE) result (dpc)
!
! Routine to return the change in pc due to a change in dE calculated in such a way as to avoid round-off error.
!-

function dpc_given_dE(pc_old, mass, dE) result(dpc)

use bmad_struct

implicit none

real(rp) pc_old, mass, dE, dpc, del2

!

del2 = dE**2 + 2 * sqrt(pc_old**2 + mass**2) * dE
!!dpc = del2 / (sqrt(pc_old**2 + del2) + pc_old)
dpc = sqrt(pc_old**2 + del2) - pc_old

end function dpc_given_dE

