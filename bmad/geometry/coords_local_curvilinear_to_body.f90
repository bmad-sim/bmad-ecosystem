!+
! Function coords_local_curvilinear_to_body (local_position, ele, w_mat, calculate_angles) result (body_position)
!
! Returns body frame coords relative to the exit end.
!
! Input:
!   local_position    -- floor_position_struct: local coordinates.
!     %r(3)               [x, y, s] position with s = Position from entrance end of element.
!   ele               -- ele_struct: element that coordinates are relative to.
!   calculate_angles  -- logical, optional: calculate angles for body_position 
!                          Default: True.
!                          False returns body_position angles (%theta, %phi, %psi) = 0.
!  
! Output         
!   body_position   -- floor_position_struct: Element coordinates relative to exit of the element.
!     %r(3)               [x, y, s] position with s = Position from entrance end of element.
!   w_mat(3,3)      -- real(rp), optional: W matrix at to transform vectors. 
!                                  v_local  = w_mat . v_body
!                                  v_body   = transpose(w_mat) . v_local
!-

function coords_local_curvilinear_to_body (local_position, ele, w_mat, calculate_angles) result(body_position)

use bmad_interface, dummy => coords_local_curvilinear_to_body

implicit none

type (floor_position_struct) :: body_position, local_position
type (ele_struct) :: ele
real(rp) :: s, g, ds, ref_tilt, theta, Sb(3,3), Sb2(3,3), Lb(3)
real(rp) ::  L_mis(3), S_mis(3,3), S_mat0(3,3)
real(rp), optional :: w_mat(3,3)
logical, optional :: calculate_angles
!

body_position = local_position
s = body_position%r(3)

select case (ele%key)
case (patch$)
  if (present(w_mat)) call mat_make_unit(w_mat)

case (sbend$, rf_bend$)
  ! Get coords relative to center
  g = ele%value(g$)
  ref_tilt = ele%value(ref_tilt_tot$)
  theta = g*s

  ! Transform from element end frame to center frame
  body_position%r(3) = 0
  body_position = bend_shift(body_position, g, ele%value(L$)/2-s, w_mat = Sb, ref_tilt = ref_tilt)

  ! Misalign
  call ele_misalignment_L_S_calc(ele, L_mis, S_mis)
  call mat_inverse(S_mis, S_mis)
  body_position%r = matmul(s_mis, body_position%r - L_mis)
  body_position%W = matmul(s_mis, body_position%W)
  Sb = matmul(S_mis, Sb)


  ! Transform to beginning
  body_position = bend_shift(body_position, g, -ele%value(L$)/2, w_mat = Sb2, ref_tilt = ref_tilt)
  Sb = matmul(Sb2, Sb)

  ! Put into tilted frame
  call rotate_vec(body_position%r, z_axis$, -ref_tilt)
  call rotate_mat(body_position%w, z_axis$, -ref_tilt)
  call rotate_mat(Sb, z_axis$, -ref_tilt)

  ! Transform to position with %r(3) = 0 and then shift %r(3) to be s-position from entrance.
  if (g /= 0) then
    ds = atan(body_position%r(3) * g / (1 + body_position%r(1) * g)) / g
    body_position = bend_shift(body_position, g, ds, Sb2)
    body_position%r(3) = ds
    Sb = matmul(Sb2, Sb)
  endif

  if (present(w_mat)) w_mat = transpose(Sb)

case default
  ! Add offsets
  body_position%r = body_position%r - [ele%value(x_offset_tot$), ele%value(y_offset_tot$), ele%value(z_offset_tot$)]

  if (ele%value(x_pitch_tot$) /= 0 .or. ele%value(y_pitch_tot$) /= 0 .or. ele%value(tilt_tot$) /= 0) then
    ! need to rotate element frame to align with local frame. 
    call floor_angles_to_w_mat (ele%value(x_pitch_tot$), ele%value(y_pitch_tot$), ele%value(tilt_tot$), w_mat_inv = S_mis)
    body_position%r(3) = body_position%r(3) - ele%value(L$)/2 ! Angle offsets are relative to the center of the element
    body_position%r = matmul(S_mis, body_position%r)
    body_position%w = matmul(S_mis, body_position%w)
    body_position%r(3) = body_position%r(3) + ele%value(L$)/2 ! Shift relative to end of element for output
  else
    ! Just shift relative to end
    call mat_make_unit(S_mis)
    body_position%r(3) = body_position%r(3)
  endif
  
  if (present(w_mat)) w_mat = transpose(S_mis)
end select

! If angles are not needed, just return zeros; 
if (logic_option(.true., calculate_angles)) then
  call update_floor_angles(body_position)
else
  body_position%theta = 0
  body_position%phi = 0
  body_position%psi = 0
endif 

end function coords_local_curvilinear_to_body
