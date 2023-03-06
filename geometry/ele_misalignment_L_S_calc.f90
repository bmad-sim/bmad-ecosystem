!+
! Subroutine ele_misalignment_L_S_calc (ele, L_mis, S_mis)
! 
! Calculates transformation vector L_mis and matrix S_mis due to misalignments for an ele
! Used to transform coordinates and vectors relative to the center of the element
!
! Input:
!   ele         -- real(rp): Element
!
! Output:
!   L_mis(3)    -- real(rp): Misalignment vector relative to center of element
!   S_mis(3,3)  -- real(rp): Misalignment matrix relative to center of element
!
!-

subroutine ele_misalignment_L_S_calc (ele, L_mis, S_mis)

use bmad_interface, dummy => ele_misalignment_L_S_calc

implicit none

type(ele_struct) :: ele 
real(rp) :: Lc(3), Sb(3,3), s0
real(rp) :: L_mis(3), S_mis(3,3)

!

L_mis = [ele%value(x_offset_tot$), ele%value(y_offset_tot$), ele%value(z_offset_tot$)]

select case(ele%key)
case(sbend$, rf_bend$)
  ! L_mis at ele center:
  ! L_mis = L_offsets + [Rz(roll) - 1] . Rz(tilt) . Ry(bend_angle/2) . rho . [cos(bend_angle/2) -1, 0, sin(bend_angle/2)]
  if (ele%value(roll_tot$) /= 0) then
    Lc = ele%value(rho$) * [cos_one(ele%value(angle$)/2), 0.0_rp, sin(ele%value(angle$)/2)]
    call rotate_vec(Lc, y_axis$, ele%value(angle$)/2)  ! rotate to entrance about y axis by half angle 
    call rotate_vec(Lc, z_axis$, ele%value(ref_tilt_tot$)) ! rotate about z axis by tilt
    L_mis = L_mis - Lc
    call rotate_vec(Lc, z_axis$, ele%value(roll_tot$))        ! rotate about z axis for roll
    L_mis = L_mis + Lc
  endif

  ! S_mis at ele center
  call floor_angles_to_w_mat (ele%value(x_pitch$), ele%value(y_pitch$), ele%value(roll_tot$), s_mis)
  
case default
    call floor_angles_to_w_mat (ele%value(x_pitch_tot$), ele%value(y_pitch_tot$), ele%value(tilt_tot$), S_mis)
end select

end subroutine ele_misalignment_L_S_calc
