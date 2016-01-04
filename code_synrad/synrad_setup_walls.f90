!+
! Subroutine synrad_setup_walls (walls, branch, seg_len_max, seg_len_phantom_max)
!
! Routine initialize the vacuum chamber walls structure after the wall points have been setup.
!
! Input:
!   walls               -- Walls_struct: wall structure with wall points already setup.
!   branch              -- branch_struct: lattice branch to use
!   seg_len_phantom_max -- real(rp), optional: If present then use this number for phantom segments
!                             instead of seg_len_max.
!
! Output:
!   walls   -- Walls_struct: wall structure.
!-

subroutine synrad_setup_walls (walls, branch, seg_len_max, seg_len_phantom_max)

use synrad_mod, except => synrad_setup_walls

implicit none

type (walls_struct), target :: walls
type (branch_struct) branch
type (wall_struct), pointer :: plus_side, minus_side

real(rp) end_n_x, end_p_x, seg_len_max
real(rp), optional :: seg_len_phantom_max

!

plus_side  => walls%positive_x_wall
minus_side => walls%negative_x_wall

plus_side%side = positive_x$
minus_side%side = negative_x$

walls%s_max = branch%param%total_length
walls%lat_geometry = branch%param%geometry

! check that endpoints are correct

if (abs(plus_side%pt(plus_side%n_pt_max)%s - walls%s_max) > 0.01) then
  print *, 'Note: Plus_side wall ends at:', plus_side%pt(plus_side%n_pt_max)%s
  print *, '      And not at lattice end of:', walls%s_max
  print *, '      [But last point is always adjusted to have s = s_lat]'
endif

if (abs(minus_side%pt(minus_side%n_pt_max)%s - walls%s_max) > 0.01) then
  print *, 'Note: Minus_side wall ends at:', minus_side%pt(minus_side%n_pt_max)%s
  print *, '      And not at lattice end of:', walls%s_max
  print *, '      [But last point is always adjusted to have s = s_lat]'
endif

plus_side%pt(plus_side%n_pt_max)%s = walls%s_max
minus_side%pt(minus_side%n_pt_max)%s = walls%s_max

! delete overlapping wall points and add wall points as needed so 
! that the triangle points can be added.

call synrad_adjust_wall_points (plus_side, branch)
call synrad_adjust_wall_points (minus_side, branch)

! segment wall

call break_wall_into_segments (minus_side, seg_len_max, branch, seg_len_phantom_max)
call break_wall_into_segments (plus_side, seg_len_max, branch, seg_len_phantom_max)

! Wall ends

walls%start_end%s = 0
walls%start_end%x = plus_side%pt(0)%x
walls%start_end%len = abs(plus_side%pt(0)%x - minus_side%pt(0)%x)

end_n_x = minus_side%seg(minus_side%n_seg_max)%x
end_p_x = plus_side%seg(plus_side%n_seg_max)%x

walls%exit_end%s = walls%s_max
walls%exit_end%x = end_n_x
walls%exit_end%len = end_p_x - end_n_x 

! do some checking

call check_wall (minus_side, walls%s_max, walls%lat_geometry)
call check_wall (plus_side, walls%s_max, walls%lat_geometry)

end subroutine
