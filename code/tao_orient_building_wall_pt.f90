!+
! Function tao_oreint_building_wall_pt(pt_in) result (pt_out)
!
! Routine to apply the building wall orientation parameters from s%building_wall%orientation 
! to a wall point.
!
! Input:
!   pt_in   -- tao_building_wall_point_struct: Building wall point.
!
! Output:
!   pt_out  -- tao_building_wall_point_struct: Building wall point with orientation params applied.
!-

function tao_oreint_building_wall_pt(pt_in) result (pt_out)

use tao_struct

implicit none

type (tao_building_wall_point_struct) pt_in, pt_out
real(rp) ct, st

!

if (s%building_wall%orientation%theta == 0 .and. s%building_wall%orientation%x_offset == 0 .and. &
                                                       s%building_wall%orientation%z_offset == 0) then
  pt_out = pt_in
  return
endif

ct = cos(s%building_wall%orientation%theta)
st = sin(s%building_wall%orientation%theta)

pt_out%z = s%building_wall%orientation%z_offset + pt_in%z * ct - pt_in%x * st
pt_out%x = s%building_wall%orientation%x_offset + pt_in%z * st + pt_in%x * ct 
pt_out%radius = pt_in%radius
pt_out%z_center = s%building_wall%orientation%z_offset + pt_in%z_center * ct - pt_in%x_center * st
pt_out%x_center = s%building_wall%orientation%x_offset + pt_in%z_center * st + pt_in%x_center * ct 

end function tao_oreint_building_wall_pt
