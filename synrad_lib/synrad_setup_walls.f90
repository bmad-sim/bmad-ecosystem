!+
! Subroutine synrad_setup_walls (walls, branch, seg_len_max)
!
! Routine initialize the vacuum chamber walls structure after the wall points have been setup.
!
! Input:
!   walls   -- Walls_struct: wall structure with wall points already setup.
!   branch  -- branch_struct: lattice branch to use
!
! Output:
!   walls   -- Walls_struct: wall structure.
!-

subroutine synrad_setup_walls (walls, branch, seg_len_max)

use synrad_mod, except => synrad_setup_walls

implicit none

type (walls_struct), target :: walls
type (branch_struct) branch
type (wall_struct), pointer :: outside, inside

real(rp) end_n_x, end_p_x, seg_len_max

!

outside => walls%positive_x_wall
inside  => walls%negative_x_wall

outside%side = positive_x$
inside%side = negative_x$

walls%s_max = branch%param%total_length
walls%lat_geometry = branch%param%geometry

! check that endpoints are correct

if (abs(outside%pt(outside%n_pt_max)%s - walls%s_max) > 0.01) then
  print *, 'Note: Outside wall ends at:', outside%pt(outside%n_pt_max)%s
  print *, '      And not at lattice end of:', walls%s_max
  print *, '      [But last point is always adjusted to have s = s_lat]'
endif

if (abs(inside%pt(inside%n_pt_max)%s - walls%s_max) > 0.01) then
  print *, 'Note: Inside wall ends at:', inside%pt(inside%n_pt_max)%s
  print *, '      And not at lattice end of:', walls%s_max
  print *, '      [But last point is always adjusted to have s = s_lat]'
endif

outside%pt(outside%n_pt_max)%s = walls%s_max
inside%pt(inside%n_pt_max)%s = walls%s_max

! delete overlapping wall points and add wall points as needed so 
! that the triangle points can be added.

call synrad_adjust_wall_points (outside, branch)
call synrad_adjust_wall_points (inside, branch)

! segment wall

call break_wall_into_segments (inside, seg_len_max, branch)
call break_wall_into_segments (outside, seg_len_max, branch)

! Wall ends

walls%start_end%s = 0
walls%start_end%x = walls%positive_x_wall%pt(0)%x
walls%start_end%s_mid = 0
walls%start_end%x_mid = walls%negative_x_wall%pt(0)%x + walls%start_end%len/2.
walls%start_end%len = abs(walls%positive_x_wall%pt(0)%x - walls%negative_x_wall%pt(0)%x)

end_n_x = walls%negative_x_wall%seg(walls%negative_x_wall%n_seg_max)%x
end_p_x = walls%positive_x_wall%seg(walls%positive_x_wall%n_seg_max)%x

walls%exit_end%s = walls%s_max
walls%exit_end%x = end_n_x
walls%exit_end%s_mid = walls%s_max
walls%exit_end%x_mid = (end_n_x + end_p_x)/2.
walls%exit_end%len = end_p_x - end_n_x 

! do some checking

call check_wall (inside, walls%s_max, walls%lat_geometry)
call check_wall (outside, walls%s_max, walls%lat_geometry)

end subroutine
