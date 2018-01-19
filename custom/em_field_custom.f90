!+
! Subroutine em_field_custom (ele, param, s_rel, orbit, local_ref_frame, field, calc_dfield, err_flag, calc_potential, grid_allow_s_out_of_bounds, rf_time)
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
!   ele              -- Ele_struct: Custom element.
!   param            -- lat_param_struct: Lattice parameters.
!   s_rel            -- Real(rp): Longitudinal position relative to the start of the element.
!   orbit            -- Coord_struct: Coords with respect to the reference particle.
!   local_ref_frame  -- Logical, If True then take the 
!                          input coordinates and output fields as being with 
!                          respect to the frame of referene of the element. 
!   calc_dfield      -- Logical, optional: If present and True then the field 
!                          derivative matrix is wanted by the calling program.
!   calc_potential   -- logical, optional: Calc electric and magnetic potentials? Default is false. 
!   grid_allow_s_out_of_bounds 
!                    -- logical, optional: For grids, allow s-coordinate to be grossly out of bounds 
!                          and return zero instead of an error? Default: False. Used internally for overlapping fields.
!   rf_time          -- real(rp), optional: Set the time relative to the RF clock. Normally this time is calculated using
!                          orbit%t or orbit%vec(5) but sometimes it is convenient to be able to override this.
!                          For example, time_runge_kutta uses this.
!
! Output:
!   field    -- Em_field_struct: Structure hoding the field values.
!   err_flag -- Logical, optional: Set true if there is an error. False otherwise.
!-

subroutine em_field_custom (ele, param, s_rel, orbit, local_ref_frame, field, calc_dfield, err_flag, calc_potential, grid_allow_s_out_of_bounds, rf_time)

use bmad_struct
use bmad_interface, except_dummy => em_field_custom

implicit none

type (ele_struct) :: ele
type (lat_param_struct) param
type (coord_struct), intent(in) :: orbit
type (em_field_struct) :: field

real(rp), intent(in) :: s_rel
real(rp), optional :: rf_time

logical local_ref_frame
logical, optional :: calc_dfield, err_flag, calc_potential, grid_allow_s_out_of_bounds

character(*), parameter :: r_name = 'em_field_custom'

!

call out_io (s_fatal$, r_name, 'THIS DUMMY ROUTINE SHOULD NOT HAVE BEEN CALLED IN THE FIRST PLACE.')
if (present(err_flag)) err_flag = .true.


end subroutine
