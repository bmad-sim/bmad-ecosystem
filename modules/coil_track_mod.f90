module coil_track_mod

contains

subroutine b_field_mult (ring, coord, first, last, s_pos, b_vector)

  use bmad
  implicit none

  type (ring_struct)  ring
  type (coord_struct)  coord

  real(rdef) s_pos(:), b_vector(3)
  real(rdef) b_loop(3)

  integer first, last, loop$, i

!

  loop$ = 999  ! temp so compiler will not complain

  b_vector = 0.0

  do i = first, last
    if (ring%ele_(i)%key == loop$) then
      call b_field_loop (coord, ring%ele_(i), s_pos(i - first + 1), b_loop)
      b_vector = b_vector + b_loop
    elseif (ring%ele_(i)%key /= drift$) then
      type *, 'ERROR IN B_FIELD_MULT.F:'
      type *, 'COIL CONTAINER DOES NOT KNOW HOW TO USE A ',  &
        key_name(ring%ele_(i)%key),'.'
      type *, 'EXITING.'
      call err_exit
    endif
  enddo

end subroutine

!+
! Subroutine B_FIELD_LOOP (COORD, ELE, S_POS, B_LOOP)
!
!   Subroutine to calculate the magnetic field vector due to a circular current
! loop.
! -- Created by Daniel Fromowitz, September 1999.
!
! Input:
!     COORD -- Coord_struct: TRUE Cartesian coordinates of particle, i.e. z is
!                            relative to the (linac) origin; it is not a
!                            displacement!
!     ELE   -- Ele_struct: Element
!     S_POS -- Real(rdef): Longitudinal position of coil component
!
! Output:
!     B_LOOP(3) -- Real(rdef): (Cartesian) Magnetic field vector x, y, and z
!                        components (in units of mu_0 / 2)
!-

!$Id$
!$Log$
!Revision 1.1  2002/06/13 14:54:59  dcs
!Interfaced with FPP/PTC
!
!Revision 1.3  2002/02/23 20:32:10  dcs
!Double/Single Real toggle added
!
!Revision 1.2  2001/09/27 18:31:48  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine b_field_loop (coord, ele, s_pos, b_loop)
  use bmad
  implicit none

  type (coord_struct)  coord
  type (ele_struct)  ele
  real(rdef) s_pos, b_loop(3)

  real(rdef) coef, e, e2_inv, inv_rot_pyr(3,3), r, rel_coords(3)
  real(rdef) rot_py(3,3), tworx, v1, v2, x_rel, y_rel, z_rel
  real(rdef) pitch, cos_p, sin_p, yaw, cos_y, sin_y, roll, cos_r, sin_r
  real(rdef) diameter, r2, ri, r2i

  integer x$, y$, z$
  parameter (x$ = 1, y$ = 2, z$ = 3)



!------------------------------------------------------------
! This subroutine needs to be rewritten.  -- DCS

  r = ele%value(radius$)
  diameter = 2 * r
  r2 = r**2
  ri = r * ele%value(current$)
  r2i = r2 * ele%value(current$)

  rel_coords(1) = coord%x%pos - ele%value(x_offset$)
  rel_coords(2) = coord%y%pos - ele%value(y_offset$)
  rel_coords(3) = coord%z%pos - s_pos

! Pitch and subsequent yaw rotation matrix:
! Note: it is assumed that the original (unrotated) loop was _first_ yawed and
! _then_ pitched (about space-fixed axes) to its actual rotated orientation.

  pitch = ele%value(y_pitch$)
  yaw = ele%value(x_pitch$)
  sin_p = sin(pitch)
  cos_p = cos(pitch)
  sin_y = sin(yaw)
  cos_y = cos(yaw)

  rot_py(1,1) = cos_y
  rot_py(1,2) = -sin_y * sin_p
  rot_py(1,3) = sin_y * cos_p
  rot_py(2,1) = 0
  rot_py(2,2) = cos_p
  rot_py(2,3) = sin_p
  rot_py(3,1) = -sin_y
  rot_py(3,2) = -cos_y * sin_p
  rot_py(3,3) = cos_y * cos_p

  rel_coords = matmul (rot_py, rel_coords)

! Roll the temporary particle coordinates about the loop axis so that the
! y-coordinate relative to the loop axis is zero.

  x_rel = rel_coords(1)
  y_rel = rel_coords(2)

  if ((y_rel == 0.0) .and. (x_rel == 0.0)) then
    roll = 0.0
  else
    roll = atan2(-y_rel, x_rel)
  endif
  sin_r = sin(roll)
  cos_r = cos(roll)
  x_rel = x_rel * cos_r - y_rel * sin_r
  if (x_rel/r > 0.001) then
    z_rel = rel_coords(3)
    tworx = diameter * x_rel
    e = (x_rel**2 + z_rel**2 + r2) / tworx
    e2_inv = 1.0/e**2
    coef = ri / tworx**1.5
    v1 = hypergeom(1, e2_inv) / e**1.5
    v2 = hypergeom(2, e2_inv) / e**2.5
    b_loop(x$) = coef * z_rel * v2
    b_loop(z$) = coef * (r * v1 - x_rel * v2)
  else
    b_loop(x$) = 0.0
    b_loop(z$) = r2i / (rel_coords(3)**2 + r2)**1.5
  endif
  b_loop(y$) = 0.0

! Roll, yaw, and pitch rotation matrix for the B-field:

  inv_rot_pyr(1,1) = cos_y * cos_r
  inv_rot_pyr(1,2) = cos_y * sin_r
  inv_rot_pyr(1,3) = sin_y
  inv_rot_pyr(2,1) = -sin_p * sin_y * cos_r - cos_p * sin_r
  inv_rot_pyr(2,2) = -sin_p * sin_y * sin_r + cos_p * cos_r
  inv_rot_pyr(2,3) = sin_p * cos_y
  inv_rot_pyr(3,1) = -cos_p * sin_y * cos_r + sin_p * sin_r
  inv_rot_pyr(3,2) = -cos_p * sin_y * sin_r - sin_p * cos_r
  inv_rot_pyr(3,3) = cos_p * cos_y

! Rotate the magnetic field vector to cancel the roll, yaw, and pitch rotations.

  b_loop = matmul (inv_rot_pyr, b_loop)

end subroutine

end module


