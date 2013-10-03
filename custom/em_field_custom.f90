!+
! Subroutine em_field_custom (ele, param, s_rel, time, here, local_ref_frame, field, calc_dfield, err_flag)
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
!   time   -- Real(rp): Particle time.
!                 For absolute time tracking this is the absolute time.
!                 For relative time tracking this is relative to the reference particle entering the element.
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

subroutine em_field_custom (ele, param, s_rel, time, orb, local_ref_frame, field, calc_dfield, err_flag)

use bmad_struct
use bmad_interface, except_dummy => em_field_custom

implicit none

type (ele_struct) :: ele
type (lat_param_struct) param
type (coord_struct), intent(in) :: orb
real(rp), intent(in) :: s_rel, time
logical local_ref_frame
type (em_field_struct) :: field
logical, optional :: calc_dfield, err_flag
character(32) :: r_name = 'em_field_custom'

!

call out_io (s_fatal$, r_name, 'THIS DUMMY ROUTINE SHOULD NOT HAVE BEEN CALLED IN THE FIRST PLACE.')
if (present(err_flag)) err_flag = .true.


end subroutine
