module photon_utils

use photon_struct

contains

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine wall_at_s (wall, s, wall_pt)
!
! Routine to calculate the wall dimensions at a given longitudinal point s.
!
! Modules needed:
!   use photon_utils
!
! Input:
!   wall -- wall_2d_struct: Wall
!   s    -- Real(rp): Longitudinal position.
!
! Output:
!   wall_pt -- wall_2d_pt_struct: Wall dimensions at s.
!-

subroutine wall_at_s (wall, s, wall_pt)

implicit none

type (wall_2d_struct) wall
type (wall_2d_pt_struct) wall_pt

real(rp) s, r
integer ix

!

call bracket_index (wall%pt%s, 0, wall%n_pt_max, s, ix)

if (ix == wall%n_pt_max) then
  wall_pt = wall%pt(ix)
else
  wall_pt%s = s
  r = (s - wall%pt(ix)%s) / (wall%pt(ix+1)%s - wall%pt(ix)%s)
  wall_pt%width2 = (1 - r) * wall%pt(ix)%width2 + r * wall%pt(ix)%width2
endif

end subroutine

end module
