!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!+
! Function average_twiss (frac1, twiss1, twiss2) result (ave_twiss)
!
! Routine to average twiss parameters.
!
! Input:
!   frac1  -- real(rp): Fraction of twiss1 to use in the average.
!   twiss1 -- twiss_struct: Twiss parameters to average.
!   twiss1 -- twiss_struct: Twiss parameters to average.
!
! Output:
!   ave_twiss -- twiss_struct: Average twiss.
!-

function average_twiss (frac1, twiss1, twiss2) result (ave_twiss)

use bmad_struct

implicit none

type (twiss_struct) twiss1, twiss2, ave_twiss
real(rp) frac1

!

ave_twiss%beta      = frac1 * twiss1%beta      + (1-frac1) * twiss2%beta
ave_twiss%alpha     = frac1 * twiss1%alpha     + (1-frac1) * twiss2%alpha
ave_twiss%gamma     = frac1 * twiss1%gamma     + (1-frac1) * twiss2%gamma
ave_twiss%phi       = frac1 * twiss1%phi       + (1-frac1) * twiss2%phi
ave_twiss%eta       = frac1 * twiss1%eta       + (1-frac1) * twiss2%eta
ave_twiss%etap      = frac1 * twiss1%etap      + (1-frac1) * twiss2%etap
ave_twiss%deta_ds   = frac1 * twiss1%deta_ds   + (1-frac1) * twiss2%deta_ds
ave_twiss%sigma     = frac1 * twiss1%sigma     + (1-frac1) * twiss2%sigma
ave_twiss%sigma_p   = frac1 * twiss1%sigma_p   + (1-frac1) * twiss2%sigma_p
ave_twiss%emit      = frac1 * twiss1%emit      + (1-frac1) * twiss2%emit
ave_twiss%norm_emit = frac1 * twiss1%norm_emit + (1-frac1) * twiss2%norm_emit

end function average_twiss
