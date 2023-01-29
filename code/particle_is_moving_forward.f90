!+
! Function particle_is_moving_forward (orbit) result (is_moving_forward)
!
! Routine to determine if a particle is moving in the forward +s direction.
! Note: By convention, reversing orbit%time_dir reverses the motion.
! If not moving forward it is dead or is moving backward.
!
! Remember: +s and +z directions are counteraligned if element being tracked 
! through is reversed.
!
! Input:
!   orbit     -- coord_struct: Particle coordinates
!
! Output:
!   is_moving_forward -- Logical: True if moving forward. False otherwise.
!-

function particle_is_moving_forward (orbit) result (is_moving_forward)

use bmad_struct

implicit none

type (coord_struct) orbit
logical is_moving_forward

!

is_moving_forward = (orbit%state == alive$) .and. (orbit%direction*orbit%time_dir == 1)

end function particle_is_moving_forward
