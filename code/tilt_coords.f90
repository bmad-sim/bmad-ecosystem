!+
! Subroutine tilt_coords (tilt_val coord, set)
!
! Subroutine to effectively tilt (rotate in the x-y plane) an element
! by instead rotating the particle position with negative the angle.
! That is, SET = True rotates the particle from the lab coords to the
! unrotated element frame.
! This subroutine is usually used with the subroutine offset_coords_m.
!
! Modules needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!   tilt_val  -- Real: Tilt value (could be the roll value for a bend)
!   set       -- Logical: If .true. then rotate the element.
!                           If .false. then unrotate the element.
!
! Output:
!   coord -- Real(6): Coordinates of particles.
!
! Note: with SET = .false. It is assumed that the rotation matrix has
!     been computed with a previous call with SET = .true.
!-


subroutine tilt_coords (tilt_val, coord, set)

  implicit none

  real cos_ang, sin_ang, rot_mat(2,2)
  real tilt_val, coord(6)

  logical set

  save cos_ang, sin_ang

!

  if (tilt_val == 0) return

  if (set) then
    sin_ang = sin(tilt_val)
    cos_ang = cos(tilt_val)
    rot_mat(1,1) =  cos_ang
    rot_mat(1,2) =  sin_ang
    rot_mat(2,1) = -sin_ang
    rot_mat(2,2) =  cos_ang
  else
    rot_mat(1,2) = -sin_ang
    rot_mat(2,1) =  sin_ang
  endif

  coord(1:3:2) = matmul(rot_mat, coord(1:3:2))
  coord(2:4:2) = matmul(rot_mat, coord(2:4:2))
                        
end subroutine
