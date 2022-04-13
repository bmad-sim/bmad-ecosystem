!+
! Function strong_beam_strength (ele) result (strength)
!
! Routine to calculate the strength (charge * number of particles) of the strong beam 
! that is part of the beambeam element.
!
! Input:
!   ele   -- ele_struct: Beambeam element.
!
! Output:
!   strength  -- real(rp): Strong beam strength.
!-

function strong_beam_strength (ele) result (strength)

use bmad_struct

implicit none

type (ele_struct) ele
real(rp) strength

! Historically the number of particles was stored in lat%param%n_part.

strength = ele%value(n_particle$) * ele%value(charge$)
if (strength == 0 .and. associated(ele%branch)) strength = ele%branch%lat%param%n_part * ele%value(charge$)

end function
