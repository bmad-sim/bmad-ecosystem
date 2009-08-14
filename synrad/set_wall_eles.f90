!+
! Subroutine set_wall_eles (wall, lat)
!
! Routine to find the lattice element and twiss parameters at each wall segment.
!-

subroutine set_wall_eles (wall, lat)

use synrad_struct
use synrad_interface

implicit none

type (wall_struct) :: wall
type (lat_struct) :: lat
type (ele_struct) :: ele

integer i, is, ix_ele

!

do is = 1, wall%n_seg_tot

  ! find the element at the s midpoint of the wall segment
  call ele_at_s (lat, wall%seg(is)%s_mid, ix_ele)
  wall%seg(is)%ix_ele = ix_ele

  ! Get the betas at the s midpoint
  call twiss_and_track_at_s (lat, wall%seg(is)%s_mid, ele)
  wall%seg(is)%a = ele%a
  wall%seg(is)%b = ele%b
  
end do

end subroutine set_wall_eles
