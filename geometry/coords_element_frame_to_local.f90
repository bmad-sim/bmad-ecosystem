
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function coords_element_frame_to_local (body_position, ele, w_mat, calculate_angles) result (local_position)
!
! Returns Cartesian coordinates relative to the exit end lab frame floor coordinates (which
! is ele%floor if ele%orientation == 1).
! This routine takes into account element misalignments.
!
! Input:
!   body_position     -- floor_position_struct: Element body frame coordinates.
!     %r                  [x, y, s] position with s = Position from entrance end of element .
!   ele               -- ele_struct: element that local_position coordinates are relative to.
!   calculate_angles  -- logical, optional: calculate angles for local_position 
!                          Default: True.
!                          False returns local_position angles (%theta, %phi, %psi) = 0.
!  
! Output         
!   local_position  -- floor_position_struct: Cartesian coordinates relative to exit of the element.
!   w_mat(3,3)      -- real(rp), optional: W matrix at to transform vectors. 
!                                  v_local     = w_mat . v_ele_frame
!                                  v_ele_frame = transpose(w_mat) . v_local
!-

function coords_element_frame_to_local (body_position, ele, w_mat, calculate_angles) result(local_position)

use bmad_interface, dummy => coords_element_frame_to_local

implicit none

type (floor_position_struct) :: body_position, local_position
type (ele_struct) :: ele
real(rp) :: s, g, theta, Sb(3,3), Lb(3)
real(rp) ::  L_mis(3), S_mis(3,3), S_mat0(3,3)
real(rp), optional :: w_mat(3,3)
logical, optional :: calculate_angles
!

local_position = body_position
s = body_position%r(3)

if (ele%key == sbend$) then
  ! Get coords relative to center
  theta = ele%value(g$)*s

  ! In ele frame at s relative to beginning. 
  ! Center coordinates at s, and move to center frame by ds = L/2 -s
  local_position%r(3) = 0
  local_position = bend_shift(local_position, ele%value(g$), ele%value(L$)/2 - s)
  
  ! Put into tilted frame
  call rotate_vec(local_position%r, z_axis$, ele%value(ref_tilt_tot$))
  call rotate_mat(local_position%w, z_axis$, ele%value(ref_tilt_tot$))

  ! Misalign
  call ele_misalignment_L_S_calc(ele, L_mis, S_mis)
  local_position%r = matmul(s_mis, local_position%r) + L_mis
  local_position%W = matmul(s_mis, local_position%W)

  ! Transform from center frame to element end frame
  local_position = bend_shift(local_position, ele%value(g$), ele%value(L$)/2, &
                                                                 w_mat = Sb, ref_tilt = ele%value(ref_tilt_tot$))
 
  if (present(w_mat)) then
    ! Initial rotation to ele's center frame and the tilt
    S_mat0 = w_mat_for_x_pitch (-(theta-ele%value(angle$)/2))
    call rotate_mat(S_mat0, z_axis$, ele%value(ref_tilt_tot$))
    w_mat = matmul(S_mis, S_mat0)
    w_mat = matmul(Sb, w_mat)
  endif

else
  if (ele%value(x_pitch_tot$) /= 0 .or. ele%value(y_pitch_tot$) /= 0 .or. ele%value(tilt_tot$) /= 0) then
    ! need to rotate element frame to align with local frame. 
    call floor_angles_to_w_mat (ele%value(x_pitch_tot$), ele%value(y_pitch_tot$), ele%value(tilt_tot$), S_mis)
    local_position%r(3) = local_position%r(3) - ele%value(L$)/2 ! Angle offsets are relative to the center of the element
    local_position%r = matmul(S_mis, local_position%r)
    local_position%w = matmul(S_mis, local_position%w)
    local_position%r(3)  = local_position%r(3) - ele%value(L$)/2 ! Shift relative to end of element for output
  else
    ! Just shift relative to end
    call mat_make_unit(S_mis)
    local_position%r(3) = local_position%r(3)  - ele%value(L$)
  endif
  
  ! Add offsets
  local_position%r = local_position%r + [ele%value(x_offset_tot$), ele%value(y_offset_tot$), ele%value(z_offset_tot$)]
  
  if (present(w_mat)) w_mat = S_mis
endif

! If angles are not needed, just return zeros; 
if (logic_option(.true., calculate_angles)) then
  call update_floor_angles(local_position)
else
  local_position%theta = 0
  local_position%phi = 0
  local_position%psi = 0
endif 

end function coords_element_frame_to_local
