!+
! Function patch_flips_propagation_direction (x_pitch, y_pitch) result (is_flip)
!
! Routine to determine if the propagation direction is flipped in a patch.
! This is true if the tranformation matrix element S(3,3) = cos(x_pitch) * cos(y_pitch) 
! is negative.
!
! Input:
!   x_pitch   -- Real(rp): Rotaion around y-axis
!   y_pitch   -- Real(rp): Rotation around x-axis.
!
! Output:
!   is_flip -- Logical: True if patch does a flip
!-

function patch_flips_propagation_direction (x_pitch, y_pitch) result (is_flip)

use bmad_struct

implicit none

real(rp) x_pitch, y_pitch
logical is_flip

!

is_flip = (cos(x_pitch) * cos(y_pitch) < 0)

end function patch_flips_propagation_direction 
