!+
! Subroutine check_aperture_limit_custom (orb, ele, at, param)
!
! Dummy routine. Will generate an error if called.
! A valid check_aperture_limit_custom is needed only if ele%aperture_type is set to
! custom$.
!
! Modules needed:
!   use bmad
!
! Input:
!   orb   -- Coord_struct: coordinates of a particle.
!   ele   -- Ele_struct: Element holding the aperture
!     %value(x1_limit$) -- Horizontal negative side aperture.
!     %value(x2_limit$) -- Horizontal positive side aparture.
!     %value(y1_limit$) -- Vertical negative side aperture.
!     %value(y2_limit$) -- Vertical positive side aparture.
!     %offset_moves_aperture -- If True then aperture moves with the element.
!   at    -- Integer: entrance_end$ or exit_end$
!   param -- lat_param_struct: Parameter structure
!     %aperture_limit_on -- The aperture limit is only checked if this is true.
!               The exception is when the orbit is larger than 
!               bmad_com%max_aperture_limit. In this case param%lost will
!               be set to True.
!
! Output:
!   param -- lat_param_struct: Parameter structure:
!     %lost -- Set True if the orbit is outside the aperture. 
!              Note: %lost is NOT set False if the orbit is inside 
!                the aperture.
!     %plane_lost_at   -- Integer: Plane where particle is lost:
!                           x_plane$ or y_plane$
!     %unstable_factor -- Real(rp): |orbit_amp/limit|
!-

subroutine check_aperture_limit_custom (orb, ele, at, param)

use bmad_struct

implicit none

type (coord_struct) :: orb
type (ele_struct) :: ele
type (lat_param_struct) :: param

integer at

!

print *, 'ERROR IN CHECK_APERTURE_LIMIT_CUSTOM: THIS DUMMY ROUTINE SHOULD NOT HAVE'
print *, '      BEEN CALLED IN THE FIRST PLACE.'
call err_exit

ele%value(1) = 0   ! so compiler will not complain

end subroutine
