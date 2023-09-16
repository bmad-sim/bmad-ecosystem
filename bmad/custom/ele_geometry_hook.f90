!+
! Subroutine ele_geometry_hook (floor0, ele, floor, finished, len_scale)
!
! Routine that can be customized to calculate the floor position of an element.
! This routine is always called by ele_geometry.
!
! This routine is useful, for example, to provide a way to calculate the orientation of a support
! structure that is, say, on a kinematic mount.
! 
! Input:
!   floor0        -- Starting floor coordinates at upstream end.
!                      Not used for fiducial and girder elements.
!   ele           -- Ele_struct: Element to propagate the geometry through.
!   len_scale     -- Real(rp), optional: factor to scale the length of the element.
!                       1.0_rp => Output is geometry at end of element (default).
!                       0.5_rp => Output is geometry at center of element. [Cannot be used for crystals.]
!                      -1.0_rp => Used to propagate geometry in reverse.
!
! Output:
!   floor       -- floor_position_struct: Floor position at downstream end.
!     %r(3)              -- X, Y, Z Floor position at end of element
!     %theta, phi, %psi  -- Orientation angles 
!   finished    -- logical: Set True to prevent ele_geometry from doing the geometry calculation.
!-

subroutine ele_geometry_hook (floor0, ele, floor, finished, len_scale)

use bmad_interface, dummy => ele_geometry_hook

implicit none

type (ele_struct) ele
type (floor_position_struct) floor0, floor
real(rp) len_scale
logical finished

!

finished = .false.   ! Must set this

end subroutine ele_geometry_hook
