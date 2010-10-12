!+
! Subroutine tilt_coords (tilt_val coord)
!
! Subroutine to effectively tilt (rotate in the x-y plane) an element
! by instead rotating the particle position with negative the angle.
!
! Modules needed:
!   use bmad
!
! Input:
!   tilt_val  -- Real(rp): Tilt value (could be the roll value for a bend)
!
! Output:
!   coord -- Real(:): Coordinates of particles.
!-

subroutine tilt_coords (tilt_val, coord)

use precision_def

implicit none

real(rp), save :: cos_ang, sin_ang, old_ang = 0
real(rp) tilt_val, coord(:), rot_mat(2,2)

logical set

!

if (tilt_val == 0) return

if (tilt_val == -old_ang) then
  sin_ang = -sin_ang
  old_ang = -old_ang
else if (tilt_val /= old_ang) then
  sin_ang = sin(tilt_val)
  cos_ang = cos(tilt_val)
  old_ang = tilt_val
endif

rot_mat(1,1) =  cos_ang
rot_mat(1,2) =  sin_ang
rot_mat(2,1) = -sin_ang
rot_mat(2,2) =  cos_ang

coord(1:3:2) = matmul(rot_mat, coord(1:3:2))
coord(2:4:2) = matmul(rot_mat, coord(2:4:2))
                        
end subroutine
