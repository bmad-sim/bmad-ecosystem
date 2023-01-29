!+
! Function particle_is_moving_backwards (orbit) result (is_moving_backwards)
!
! Routine to determine if a particle is moving in the backward -s direction.
! Note: By convention, reversing orbit%time_dir reverses the motion.
! If not moving backward, the particle is dead or is moving forward.
!
! Remember: +s and +z directions are counteraligned if element being tracked 
! through is reversed.
!
! Input:
!   orbit  -- coord_struct: Particle coordinates
!
! Output:
!   is_moving_backwards -- Logical: True if moving backward. False otherwise.
!-

function particle_is_moving_backwards (orbit) result (is_moving_backwards)

use bmad_struct

implicit none

type (coord_struct) orbit
logical is_moving_backwards

!

is_moving_backwards = (orbit%state == alive$) .and. (orbit%direction*orbit%time_dir == -1)

end function particle_is_moving_backwards
