!+
! Function particle_is_moving_forward (orbit, dir) result (is_moving_forward)
!
! Routine to determine if a particle is tracking in the forward +s direction.
! Note: reversing orbit%time_dir reverses the tracking direction.
! If not tracking forward it is dead or is moving backward.
!
! Remember: +s and +z directions are counteraligned if element being tracked 
! through is reversed.
!
! Input:
!   orbit     -- coord_struct: Particle coordinates
!   dir       -- integer, optional: +1 if tracking forward(default) or -1 to return True if tracking backwards.
!
! Output:
!   is_moving_forward -- Logical: True if moving forward. False otherwise.
!-

function particle_is_moving_forward (orbit, dir) result (is_moving_forward)

use bmad_struct

implicit none

type (coord_struct) orbit
integer, optional :: dir
logical is_moving_forward

!

is_moving_forward = (orbit%state == alive$) .and. (orbit%direction*orbit%time_dir == integer_option(1, dir))

end function particle_is_moving_forward
