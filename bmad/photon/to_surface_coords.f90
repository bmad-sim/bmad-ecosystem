!+
! Subroutine to_surface_coords (lab_orbit, ele, surface_orbit)
!
! Routine to convert lab_orbit laboratory coordinates to surface body coordinates.
!
! Input:
!   lab_orbit     -- coord_struct: Photon position in laboratory coords.
!   ele           -- ele_struct: Detector element.
!
! Output:
!   surface_orbit -- coord_struct: Photon position in element body coordinates.
!     %state      -- Set to lost$ if orbit outside of surface (can happen with sperical surface).
!-

subroutine to_surface_coords (lab_orbit, ele, surface_orbit)

use photon_utils_mod, dummy => to_surface_coords
use super_recipes_mod, only: super_qromb

implicit none

type (coord_struct) lab_orbit, surface_orbit
type (ele_struct) ele

real(rp) w_surf(3,3), z
integer idim
logical err_flag

!

surface_orbit = lab_orbit
call offset_photon (ele, surface_orbit, set$)  ! Go to coordinates of the detector
call track_to_surface (ele, surface_orbit, ele%branch%param, w_surf)

if (.not. has_curvature(ele%photon)) return

z = z_at_surface(ele, surface_orbit%vec(1), surface_orbit%vec(3), err_flag)
if (err_flag) then
  surface_orbit%state = lost$
  return
endif

!

idim = 1
surface_orbit%vec(1) = super_qromb(qfunc, 0.0_rp, surface_orbit%vec(1), 1e-6_rp, 1e-8_rp, 4, err_flag)

idim = 2
surface_orbit%vec(3) = super_qromb(qfunc, 0.0_rp, surface_orbit%vec(3), 1e-6_rp, 1e-8_rp, 4, err_flag)

surface_orbit%vec(5) = z_at_surface(ele, surface_orbit%vec(1), surface_orbit%vec(3), err_flag, .true.)

!------------------------------------------------------
contains

function qfunc (x) result (value)

real(rp), intent(in) :: x(:)
real(rp) :: value(size(x))
real(rp) z, xy(2), dz_dxy(2)
logical err_flag

integer i

!

xy = 0

do i = 1, size(x)
  xy(idim) = x(i)
  z = z_at_surface(ele, xy(1), xy(2), err_flag, .true., dz_dxy)
  value(i) = sqrt(1.0_rp + dz_dxy(idim)**2)
enddo

end function qfunc

end subroutine to_surface_coords

