!+
! Function coords_body_to_local (body_position, ele, w_mat, calculate_angles) result (local_position)
!
! converts from body coordinates to local (laboratory) coordinates.
!
! Input:
!   body_position     -- floor_position_struct: Element body frame coordinates.
!     %r(3)               [x, y, s] position with s = Position from entrance end of element.
!   ele               -- ele_struct: element that local_position coordinates are relative to.
!   calculate_angles  -- logical, optional: calculate angles for local_position 
!                          Default: True.
!                          False returns local_position angles (%theta, %phi, %psi) = 0.
!  
! Output         
!   local_position  -- floor_position_struct: Local laboratory coordinates.
!     %r(3)               [x, y, s] position with s = Position from entrance end of element.
!   w_mat(3,3)      -- real(rp), optional: W matrix at to transform vectors. 
!                                  v_local  = w_mat . v_body
!                                  v_body   = transpose(w_mat) . v_local
!-

function coords_body_to_local (body_position, ele, w_mat, calculate_angles) result(local_position)

use bmad_interface, dummy => coords_body_to_local

implicit none

type (floor_position_struct) :: body_position, local_position
type (ele_struct) :: ele
real(rp) :: s0, g, dd, ds, L2, ref_tilt, Sb(3,3), Sb2(3,3), Lb(3)
real(rp) ::  L_mis(3), S_mis(3,3), S_mat0(3,3)
real(rp), optional :: w_mat(3,3)
logical, optional :: calculate_angles
!

local_position = body_position
L2 = ele%value(L$) / 2

select case (ele%key)
case (patch$)
  if (present(w_mat)) call mat_make_unit(w_mat)

case (sbend$, rf_bend$)
  g = ele%value(g$)
  s0 = body_position%r(3)
  ref_tilt = ele%value(ref_tilt_tot$)

  ! Rotate to entrance frame
  local_position%r(3) = 0
  local_position = bend_shift(local_position, g, -s0, Sb)

  ! Put into tilted frame
  call rotate_vec(local_position%r, z_axis$, ref_tilt)
  call rotate_mat(local_position%w, z_axis$, ref_tilt)
  call rotate_mat(Sb, z_axis$, ref_tilt)

  ! Transfer to center coords
  local_position = bend_shift(local_position, g, L2, Sb2, ref_tilt)
  Sb = matmul(Sb2, Sb)

  ! Misalign
  call ele_misalignment_L_S_calc(ele, L_mis, S_mis)
  local_position%r = matmul(s_mis, local_position%r) + L_mis
  local_position%W = matmul(s_mis, local_position%W)
  Sb = matmul(S_mis, Sb)

  ! Transform from center frame to frame where particle is with %r(3) = 0 and then set 
  ! %r(3) to distance from entrance end
  if (g /= 0) then
    dd = (local_position%r(1) * cos(ref_tilt) + local_position%r(2) * sin(ref_tilt)) * g
    ds = atan(local_position%r(3) * g / (1 + dd)) / g
    local_position = bend_shift(local_position, g, ds, Sb2, ref_tilt)
    local_position%r(3) = ds + L2
    Sb = matmul(Sb2, Sb)
  else
    local_position%r(3) = local_position%r(3) + L2
  endif

  if (present(w_mat)) w_mat = Sb

case default
  if (ele%value(x_pitch_tot$) /= 0 .or. ele%value(y_pitch_tot$) /= 0 .or. ele%value(tilt_tot$) /= 0) then
    ! need to rotate element frame to align with local frame. 
    call floor_angles_to_w_mat (ele%value(x_pitch_tot$), ele%value(y_pitch_tot$), ele%value(tilt_tot$), S_mis)
    local_position%r(3) = local_position%r(3) - L2 ! Angle offsets are relative to the center of the element
    local_position%r = matmul(S_mis, local_position%r)
    local_position%w = matmul(S_mis, local_position%w)
    local_position%r(3) = local_position%r(3) + L2 ! Shift relative to end of element for output
  else
    ! Just shift relative to end
    S_mis = mat3_unit$
  endif
  
  ! Add offsets
  local_position%r = local_position%r + [ele%value(x_offset_tot$), ele%value(y_offset_tot$), ele%value(z_offset_tot$)]
  
  if (present(w_mat)) w_mat = S_mis
end select

! If angles are not needed, just return zeros; 
if (logic_option(.true., calculate_angles)) then
  call update_floor_angles(local_position)
else
  local_position%theta = 0
  local_position%phi = 0
  local_position%psi = 0
endif 

end function coords_body_to_local
