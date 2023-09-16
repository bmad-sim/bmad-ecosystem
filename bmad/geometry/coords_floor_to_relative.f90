!+
! Function coords_floor_to_relative (floor0, global_position, calculate_angles, is_delta_position) result (local_position)
!
! Returns local floor position relative to floor0 given a global floor position.
! This is an essentially an inverse of routine coords_relative_to_floor.
!
! Input:
!   floor0            -- floor_position_struct: reference position
!   global_position   -- floor_position_struct: global position 
!   calculate_angles  -- logical, optional: calculate angles for local_position 
!                          Default: True.
!                          False returns local_position angles (%theta, %phi, %psi) = 0.
!   is_delta_position -- logical, optional: If True then treat global_position%r as a difference
!                           position in global space and only rotate the position but not shift it.
!                           Default: False.
!
! Output:
!  local_position -- floor_position_struct: position relative to floor0
!-

function coords_floor_to_relative (floor0, global_position, calculate_angles, is_delta_position) result (local_position)

use bmad_interface, dummy => coords_floor_to_relative

implicit none

type (floor_position_struct) floor0, global_position, local_position
real(rp) :: w0_mat_T(3,3), w_mat(3,3)
logical, optional :: calculate_angles, is_delta_position

! transpose
w0_mat_T = transpose(floor0%W)

! Solve for r_local = [x, y, z]_local
   
if (logic_option(.false., is_delta_position)) then
  local_position%r = matmul(w0_mat_T, global_position%r)
else
  local_position%r = matmul(w0_mat_T, global_position%r - floor0%r)
endif

local_position%w =  matmul(w0_mat_T, global_position%w)

! If angles are not needed, just return zeros; 
if (logic_option(.true., calculate_angles)) then
  call update_floor_angles(local_position, floor0)
else
  local_position%theta = 0
  local_position%phi = 0
  local_position%psi = 0
endif 


end function coords_floor_to_relative
