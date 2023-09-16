!+
! Function gradient_shift_sr_wake (ele, param) result (grad_shift)
! 
! Function to return the shift in the accelerating gradient due to:
!   1) Short range longitudinal wake forces.
!   2) Adjustment in the external applied RF power attempting to compensate for the SR wake loss.
!
! Note: This routine will only return a non-zero value when bmad_com%sr_wakes_on = True
!
! Input:
!   ele           -- ele_struct: Lcavity element.
!   param         -- lat_param_struct: Lattice parameters
!     %n_part        -- Number of particles in a bunch
!     %particle      -- Type of particle
!
! Output:
!   grad_shift -- Real(rp): Shift in gradient
!-

function gradient_shift_sr_wake (ele, param) result (grad_shift)

use bmad_struct
implicit none

type (ele_struct) ele
type (lat_param_struct) param
real(rp) grad_shift

! 

if (bmad_com%sr_wakes_on .and. ele%value(l$) /= 0) then
  grad_shift = ele%value(e_loss$) * param%n_part * abs(charge_of(ele%ref_species)) * e_charge / ele%value(l$) 
else
  grad_shift = 0
endif

end function gradient_shift_sr_wake

