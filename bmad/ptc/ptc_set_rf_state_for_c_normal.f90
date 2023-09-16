!+
! Subroutine ptc_set_rf_state_for_c_normal (nocavity)
!
! Routine to set the appropriate PTC parameters so that routines like c_normal and
! c_full_canonise know if the RF is on or off. That is, whether pz is a constant of
! the motion (RF off) or whether there are longitudinal oscillations (RF on). 
! This will affect the analysis.
!
! Input:
!   nocavity      -- logical: True -> RF is off and vice versa.
!-

subroutine ptc_set_rf_state_for_c_normal (nocavity)

use c_tpsa, only: c_bmad_reinit, ndpt_bmad

implicit none

logical nocavity

!

if (nocavity) then
  call c_bmad_reinit(5+ndpt_bmad)
else
  call c_bmad_reinit(0)
endif

end subroutine ptc_set_rf_state_for_c_normal

