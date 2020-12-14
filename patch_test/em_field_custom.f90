!+
! Subroutine em_field_custom (ele, param, s_rel, t_rel, here, local_ref_frame, field, calc_dfield, err_flag)
!
! Routine for handling custom (user supplied) EM fields.
! This routine is called when ele%field_calc = custom$ or when ele is a custom element (ele%key = custom$)
! In order to be used, this stub file must be modified appropriately. See the Bmad manual for more details.
!
! Note: Unlike all other elements, "s_rel" and "here" areguments for a patch element are with respect to 
! the exit reference frame of the element. See the Bmad manual for more details.
!
! Note: Fields should not have any unphysical discontinuities. 
! Discontinuities may cause Runge-Kutta integration to fail resulting in particles getting marked as "lost".
! The mode of failure here is that RK will try smaller and smaller steps to integrate through the 
! discontinuity until the step size gets lower than bmad_com%min_ds_adaptive_tracking. At this
! point the particle gets marked as lost.
!
! General rule: Your code may NOT modify any argument that is not listed as
! an output agument below.
!
! Input:
!   ele         -- Ele_struct: Custom element.
!   param       -- lat_param_struct: Lattice parameters.
!   s_rel       -- Real(rp): Longitudinal position relative to the start of the element.
!   here        -- Coord_struct: Coords with respect to the reference particle.
!   local_ref_frame 
!               -- Logical, If True then take the 
!                     input coordinates and output fields as being with 
!                     respect to the frame of referene of the element. 
!   calc_dfield -- Logical, optional: If present and True then the field 
!                     derivative matrix is wanted by the calling program.
!
! Output:
!   field    -- Em_field_struct: Structure hoding the field values.
!   err_flag -- Logical, optional: Set true if there is an error. False otherwise.
!-

recursive subroutine em_field_custom (ele, param, s_rel, orb, local_ref_frame, field, calc_dfield, err_flag)

use em_field_mod, except_dummy2 => em_field_custom

implicit none

type (ele_struct), target :: ele
type (ele_struct), pointer :: ele2
type (lat_param_struct) param
type (coord_struct), intent(in) :: orb
type (coord_struct) :: orb2

real(rp), intent(in) :: s_rel
logical local_ref_frame
type (em_field_struct) :: field
logical, optional :: calc_dfield, err_flag
character(*), parameter :: r_name = 'em_field_custom'

real(rp) w_mat(3,3), w_mat_inv(3,3), r_vec(3), r0vec(3)
real(rp), pointer :: v(:)

!

if (present(err_flag)) err_flag = .false.

! Convert particle coordinates from exit to entrance frame.
v => ele%value   ! v helps makes code compact
call floor_angles_to_w_mat (v(x_pitch$), v(y_pitch$), v(tilt$), w_mat, w_mat_inv)
r0vec = [v(x_offset$), v(y_offset$), v(z_offset$)]
r_vec = [orb%vec(1), orb%vec(3), s_rel]  ! coords in exit frame
r_vec = matmul(w_mat, r_vec) +  r0vec     ! coords in entrance frame

!

ele2 => ele%branch%ele(ele%ix_ele-1)
orb2 = orb
orb2%vec(1:3:2) = r_vec(1:2)
call em_field_calc (ele2, param, ele2%value(l$)/2, orb2, .false., field, calc_dfield, err_flag)

! Convert field from entrance to exit frame
field%E = matmul(w_mat_inv, field%E)
field%B = matmul(w_mat_inv, field%B)
if (logic_option(.false., calc_dfield)) then
  field%dE = matmul(w_mat_inv, matmul(field%dE, w_mat))
  field%dB = matmul(w_mat_inv, matmul(field%dB, w_mat))
endif

end subroutine
