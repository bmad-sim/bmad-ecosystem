!+
! Function rot_2d (vec_in, angle) result (vec_out)
!
! Routine to do a 2D rotation
!
! Input:
!   vec_in(2)   -- real(rp): Init vec
!   angle       -- real(rp): angle in radians.
!
! Output:
!   vec_out(2)  -- real(rp): Rotated vec.
!-

function rot_2d (vec_in, angle) result (vec_out)

use precision_def
implicit none

real(rp) vec_in(2), angle, vec_out(2)

!

if (angle == 0) then
  vec_out = vec_in
  return
endif

vec_out = [cos(angle)*vec_in(1) - sin(angle)*vec_in(2), sin(angle)*vec_in(1) + cos(angle)*vec_in(2)]

end function
