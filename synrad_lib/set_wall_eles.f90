!+
! Subroutine set_wall_eles (wall, branch)
!
! Routine to find the lattice element and twiss parameters at each wall segment.
!-

subroutine set_wall_eles (wall, branch)

use synrad_struct
use twiss_and_track_mod

implicit none

type (wall_struct) :: wall
type (branch_struct) :: branch
type (ele_struct) :: ele

integer i, is, ix_ele

!

do is = 1, wall%n_seg_max

  ! find the element at the s midpoint of the wall segment
  ix_ele = element_at_s (branch%lat, wall%seg(is)%s_mid, .true., branch%ix_branch)
  wall%seg(is)%ix_ele = ix_ele

  ! Get the betas at the s midpoint
  call twiss_and_track_at_s (branch%lat, wall%seg(is)%s_mid, ele, ix_branch = branch%ix_branch)
  wall%seg(is)%a = ele%a
  wall%seg(is)%b = ele%b
  
end do

end subroutine set_wall_eles
