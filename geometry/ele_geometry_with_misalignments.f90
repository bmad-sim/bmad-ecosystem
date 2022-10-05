!+
! Function ele_geometry_with_misalignments (ele, len_scale) result (floor)
!
! Routine to calculate the element body coordinates (that is, coordinates with misalignments) 
! for an element at some distance s_offset from the upstream end.
!
! Note: For crystal and mirror elements, len_scale ~ 0.5 means compute the floor coordinates of the
! element surface.
!
! Input:
!   ele       -- ele_struct: Lattice element under consideration.
!   len_scale        -- Real(rp), optional: factor to scale the length of the element.
!                          1.0_rp => Output is geometry at end of element (default).
!                          0.5_rp => Output is geometry at center of element. 
!                         -1.0_rp => Used to propagate geometry in reverse.
!
! Output:
!   floor     -- floor_position_struct: Floor position with misalignments
!-

function ele_geometry_with_misalignments (ele, len_scale) result (floor)

use bmad_interface, dummy => ele_geometry_with_misalignments

implicit none

type (ele_struct), target :: ele
type (ele_struct), pointer :: ele0
type (floor_position_struct) floor, f0

real(rp), optional :: len_scale
real(rp) s_rel, s_mis(3,3)

!

select case (ele%key)

case (floor_shift$, group$, overlay$, hybrid$, beginning_ele$, match$, null_ele$, patch$)
  ! These elements do not have misalignments
  floor = ele%floor

case (crystal$, mirror$, multilayer_mirror$)
  ele0 => pointer_to_next_ele(ele, -1)
  f0 = ele0%floor
  call floor_angles_to_w_mat(ele%value(x_pitch_tot$), ele%value(y_pitch_tot$), ele%value(tilt_tot$), s_mis)
  f0%r = f0%r + matmul(f0%W, [ele%value(x_offset_tot$), ele%value(y_offset_tot$), ele%value(z_offset_tot$)])
  f0%W = matmul(f0%W, s_mis)
  call floor_w_mat_to_angles (f0%W, f0%theta, f0%phi, f0%psi, f0)
  call ele_geometry (f0, ele, floor, len_scale)

  ! Misalignments are referenced to beginning of element
  floor = coords_relative_to_floor (ele0%floor, [ele%value(x_offset_tot$), ele%value(y_offset_tot$), ele%value(z_offset_tot$)], &
                                      ele%value(x_pitch_tot$), ele%value(y_pitch_tot$), ele%value(tilt_tot$))
  call ele_geometry (floor, ele, floor, len_scale)

case (girder$)
  floor = coords_relative_to_floor (ele%floor, [ele%value(x_offset_tot$), ele%value(y_offset_tot$), ele%value(z_offset_tot$)], &
                                      ele%value(x_pitch_tot$), ele%value(y_pitch_tot$), ele%value(tilt_tot$)) 
case default
  ! Misalignments referenced to center of element
  call ele_geometry (ele%floor, ele, floor, -0.5_rp)
  floor = coords_relative_to_floor (floor, [ele%value(x_offset_tot$), ele%value(y_offset_tot$), ele%value(z_offset_tot$)], &
                                      ele%value(x_pitch_tot$), ele%value(y_pitch_tot$), ele%value(tilt_tot$)) 
  if (ele%value(l$) == 0) then
    call ele_geometry (floor, ele, floor, 0.5_rp)
  else
    s_rel = real_option(1.0_rp, len_scale) - 0.5_rp
    call ele_geometry (floor, ele, floor, s_rel)
  endif

end select

end function ele_geometry_with_misalignments
