!+
! Function e_accel_field (ele, voltage_or_gradient) result (field)
!
! Routine to return the gradient or voltage through an e_gun, lcavity or rfcavity element.
!
! Input:
!   ele                 -- ele_struct: Lcavity or rfcavity element.
!   voltage_or_gradient -- integer: voltage$ or gradient$
!
! Output:
!   field    -- real(rp): Cavity field or gradient.
!-

function e_accel_field (ele, voltage_or_gradient) result (field)

use bmad_struct

implicit none

type (ele_struct) ele
real(rp) field
integer voltage_or_gradient 

!

if (.not. ele%is_on) then
  field = 0
  return
endif

select case (ele%key)
case (lcavity$)
  select case (voltage_or_gradient)
  case (voltage$)
    field = (ele%value(voltage$) + ele%value(voltage_err$)) * ele%value(field_autoscale$)
  case (gradient$)
    field = (ele%value(gradient$) + ele%value(gradient_err$)) * ele%value(field_autoscale$)
  end select

case (rfcavity$, e_gun$, em_field$, crab_cavity$)
  select case (voltage_or_gradient)
  case (voltage$)
    field = ele%value(voltage$) * ele%value(field_autoscale$)
  case (gradient$)
    field = ele%value(gradient$) * ele%value(field_autoscale$)
  end select
end select

end function e_accel_field

