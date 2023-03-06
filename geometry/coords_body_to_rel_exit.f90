!+
! Function coords_body_to_rel_exit (body_position, ele, w_mat, calculate_angles) result (rel_exit)
!
! Returns global coordinates relative to the exit end lab frame floor coordinates (which
! is ele%floor if ele%orientation == 1).
! This routine takes into account element misalignments.
!
! Input:
!   body_position     -- floor_position_struct: Element body frame coordinates.
!     %r                  [x, y, s] position with s = Position from entrance end of element .
!   ele               -- ele_struct: element that rel_exit coordinates are relative to.
!   calculate_angles  -- logical, optional: calculate angles for rel_exit 
!                          Default: True.
!                          False returns rel_exit angles (%theta, %phi, %psi) = 0.
!  
! Output         
!   rel_exit        -- floor_position_struct: Cartesian coordinates relative to exit of the element.
!   w_mat(3,3)      -- real(rp), optional: W matrix at to transform vectors. 
!                                  v_rel_exit = w_mat . v_body
!                                  v_body     = transpose(w_mat) . v_rel_exit
!-

function coords_body_to_rel_exit (body_position, ele, w_mat, calculate_angles) result(rel_exit)

use bmad_interface, dummy => coords_body_to_rel_exit

implicit none

type (floor_position_struct) :: body_position, rel_exit
type (ele_struct) :: ele
real(rp) :: s, g, theta, Sb(3,3), Lb(3)
real(rp) ::  L_mis(3), S_mis(3,3), S_mat0(3,3)
real(rp), optional :: w_mat(3,3)
logical, optional :: calculate_angles
!

rel_exit = body_position
s = body_position%r(3)

if (ele%key == sbend$ .or. ele%key == rf_bend$) then
  ! Get coords relative to center
  theta = ele%value(g$)*s

  ! In ele frame at s relative to beginning. 
  ! Center coordinates at s, and move to center frame by ds = L/2 -s
  rel_exit%r(3) = 0
  rel_exit = bend_shift(rel_exit, ele%value(g$), ele%value(L$)/2 - s)
  
  ! Put into tilted frame
  call rotate_vec(rel_exit%r, z_axis$, ele%value(ref_tilt_tot$))
  call rotate_mat(rel_exit%w, z_axis$, ele%value(ref_tilt_tot$))

  ! Misalign
  call ele_misalignment_L_S_calc(ele, L_mis, S_mis)
  rel_exit%r = matmul(s_mis, rel_exit%r) + L_mis
  rel_exit%W = matmul(s_mis, rel_exit%W)

  ! Transform from center frame to element end frame
  rel_exit = bend_shift(rel_exit, ele%value(g$), ele%value(L$)/2, &
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
    rel_exit%r(3) = rel_exit%r(3) - ele%value(L$)/2 ! Angle offsets are relative to the center of the element
    rel_exit%r = matmul(S_mis, rel_exit%r)
    rel_exit%w = matmul(S_mis, rel_exit%w)
    rel_exit%r(3)  = rel_exit%r(3) - ele%value(L$)/2 ! Shift relative to end of element for output
  else
    ! Just shift relative to end
    call mat_make_unit(S_mis)
    rel_exit%r(3) = rel_exit%r(3)  - ele%value(L$)
  endif
  
  ! Add offsets
  rel_exit%r = rel_exit%r + [ele%value(x_offset_tot$), ele%value(y_offset_tot$), ele%value(z_offset_tot$)]
  
  if (present(w_mat)) w_mat = S_mis
endif

! If angles are not needed, just return zeros; 
if (logic_option(.true., calculate_angles)) then
  call update_floor_angles(rel_exit)
else
  rel_exit%theta = 0
  rel_exit%phi = 0
  rel_exit%psi = 0
endif 

end function coords_body_to_rel_exit
