!+
! Subroutine tilt_coords (tilt_val coord, set)
!
! Subroutine to effectively tilt (rotate in the x-y plane) an element
! by instead rotating the particle position with negative the angle.
! That is, SET = True rotates the particle from the lab coords to the
! unrotated element frame.
!
! Modules needed:
!   use bmad
!
! Input:
!   tilt_val  -- Real(rp): Tilt value (could be the roll value for a bend)
!   set       -- Logical: If .true. then rotate the element.
!                           If .false. then unrotate the element.
!
! Output:
!   coord -- Real(:): Coordinates of particles.
!
! Note: with SET = .false. It is assumed that the rotation matrix has
!     been computed with a previous call with SET = .true.
!-

!$Id$
!$Log$
!Revision 1.6  2003/07/09 01:38:20  dcs
!new bmad with allocatable ring%ele_(:)
!
!Revision 1.5  2002/02/23 20:32:25  dcs
!Double/Single Real toggle added
!
!Revision 1.4  2002/01/08 21:44:43  dcs
!Aligned with VMS version  -- DCS
!
!Revision 1.3  2001/11/29 19:39:54  helms
!Updates from DCS including (*) -> (:)
!
!Revision 1.2  2001/09/27 18:31:58  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"

subroutine tilt_coords (tilt_val, coord, set)

  use precision_def

  implicit none

  real(rp) cos_ang, sin_ang, rot_mat(2,2)
  real(rp) tilt_val, coord(:)

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
    rot_mat(1,1) =  cos_ang
    rot_mat(1,2) = -sin_ang
    rot_mat(2,1) =  sin_ang
    rot_mat(2,2) =  cos_ang
  endif

  coord(1:3:2) = matmul(rot_mat, coord(1:3:2))
  coord(2:4:2) = matmul(rot_mat, coord(2:4:2))
                        
end subroutine
