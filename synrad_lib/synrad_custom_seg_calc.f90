!+
! Subroutine synrad_custom_seg_calc (wall, ray, seg, frac_illum)
!
! Routine for doing a custom synrad calculation.
! To implement: Code a routine of this name and compile with your program.
!
! Input:
!   wall        -- wall_struct: +x or -x side of the wall.
!   ray         -- ray_struct: Photon hitting segment.
!   seg         -- wall_seg_struct: Segment being hit.
!   frac_illum  -- real(rp): Fraction of the segment illuminated.

subroutine synrad_custom_seg_calc (wall, ray, seg, frac_illum)

use synrad_interface, dummy => synrad_custom_seg_calc

implicit none

type (wall_struct) wall
type (ray_struct) ray
type (wall_seg_struct) seg

real(rp) frac_illum

! The default is to do nothing.

end subroutine
