!+
! Subroutine track_to_surface (ele, orbit, param, w_surface)
!
! Routine to track a photon to the surface of the element.
!
! After calling this routine, if the surface is curved, the routine rotate_for_curved_surface should
! be called to rotate the photon's velocity coordinates to the local surface coordinate system.
!
! Input:
!   ele            -- ele_struct: Element
!   orbit          -- coord_struct: Coordinates in the element coordinate frame
!   param          -- lat_param_struct: Branch parameters.
!
! Output:
!   orbit          -- coord_struct: At surface in local surface coordinate frame
!     %state         -- set to lost$ if the photon does not intersect the surface (can happen when surface is curved).
!   w_surface(3,3) -- real(rp), rotation matrix to transform to surface coords.
!-

subroutine track_to_surface (ele, orbit, param, w_surface)

use super_recipes_mod
use photon_utils_mod, dummy => track_to_surface

implicit none

type (ele_struct) ele
type (coord_struct) orbit
type (lat_param_struct) param

real(rp) :: w_surface(3,3), s_len, s1, s2, s_center, x0, y0, z
integer status 

character(*), parameter :: r_name = 'track_to_surface'

! If there is curvature, compute the reflection point which is where 
! the photon intersects the surface.

if (has_curvature(ele%photon)) then

  ! Assume flat crystal, compute s required to hit the intersection

  s_center = orbit%vec(5) / orbit%vec(6)

  s1 = s_center
  s2 = s_center
  z = photon_depth_in_element(s_center, status); if (status /= 0) return
  if (z > 0) then
    do
      s1 = s1 - 0.1
      z = photon_depth_in_element(s1, status); if (status /= 0) return
      if (z < 0) exit
      if (s1 < -10) then
        !! call out_io (s_warn$, r_name, &
        !!      'PHOTON INTERSECTION WITH SURFACE NOT FOUND FOR ELEMENT: ' // ele%name)
        orbit%state = lost$
        return
      endif
    enddo
  else
    do
      s2 = s2 + 0.1
      z = photon_depth_in_element(s2, status); if (status /= 0) return
      if (z > 0) exit
      if (s1 > 10) then
        !! call out_io (s_warn$, r_name, &
        !!      'PHOTON INTERSECTION WITH SURFACE NOT FOUND FOR ELEMENT: ' // ele%name)
        orbit%state = lost$
        return
      endif
    enddo
  endif

  s_len = super_zbrent (photon_depth_in_element, s1, s2, 0.0_rp, 1d-10, status)
  if (orbit%state == lost$) return

  ! Compute the intersection point

  orbit%vec(1:5:2) = s_len * orbit%vec(2:6:2) + orbit%vec(1:5:2)
  orbit%t = orbit%t + s_len / c_light

else
  s_len = -orbit%vec(5) / orbit%vec(6)
  orbit%vec(1:5:2) = orbit%vec(1:5:2) + s_len * orbit%vec(2:6:2) ! Surface is at z = 0
  orbit%t = orbit%t + s_len / c_light
endif

! Check aperture and rotate to surface coords.

if (ele%aperture_at == surface$) call check_aperture_limit (orbit, ele, surface$, param)
if (orbit%state /= alive$) return

if (has_curvature(ele%photon)) call rotate_for_curved_surface (ele, orbit, set$, w_surface)

contains
!-----------------------------------------------------------------------------------------------
!+
! Function photon_depth_in_element (s_len, status) result (delta_h)
! 
! Private routine to be used as an argument in zbrent. Propagates
! photon forward by a distance s_len. Returns delta_h = z-z0
! where z0 is the height of the element surface. 
! Since positive z points inward, positive delta_h => inside element.
!
! Input:
!   s_len   -- real(rp): Place to position the photon.
!
! Output:
!   delta_h -- real(rp): Depth of photon below surface in crystal coordinates.
!   status  -- integer: 0 -> Calculation OK.
!                       1 -> Calculation not OK.
!-

function photon_depth_in_element (s_len, status) result (delta_h)

real(rp), intent(in) :: s_len
real(rp) :: delta_h
real(rp) :: point(3)

integer status  ! Need to use status arg due to super_zbrent.
logical err_flag

! Extend_grid = True is needed since test points may well be outside of the grid.

point = s_len * orbit%vec(2:6:2) + orbit%vec(1:5:2)
delta_h = point(3) - z_at_surface(ele, point(1), point(2), err_flag, .true.)
if (err_flag) then
  orbit%state = lost$
  status = 1
else
  status = 0
endif

end function photon_depth_in_element

end subroutine track_to_surface

