!+
! Function e_accel_field (ele, voltage_or_gradient, bmad_standard_tracking) result (field)
!
! Routine to return the gradient or voltage through an e_gun, lcavity or rfcavity element.
!
! Input:
!   ele                     -- ele_struct: Lcavity or rfcavity element.
!   voltage_or_gradient     -- integer: voltage$ or gradient$
!   bmad_standard_tracking  -- logical, optional: Using bmad_standard tracking? Default is False.
!
! Output:
!   field                   -- real(rp): Cavity field or gradient.
!-

function e_accel_field (ele, voltage_or_gradient, bmad_standard_tracking) result (field)

use bmad_struct

implicit none

type (ele_struct) ele
real(rp) field, f_auto
integer voltage_or_gradient 
logical, optional :: bmad_standard_tracking

!

if (.not. ele%is_on) then
  field = 0
  return
endif

! With bmad_standard tracking field_autoscale is not used since bmad_standard tracking by design gives 
! the correct energy change. In fact, using phi0_autoscale would be a mistake 
! if, say, tracking_method = runge_kutta, mat6_calc_method = bmad_standard.

f_auto = ele%value(field_autoscale$)
if (logic_option(.false., bmad_standard_tracking)) f_auto = 1

select case (ele%key)
case (lcavity$)
  select case (voltage_or_gradient)
  case (voltage$)
    field = (ele%value(voltage$) + ele%value(voltage_err$)) * f_auto
  case (gradient$)
    field = (ele%value(gradient$) + ele%value(gradient_err$)) * f_auto
  end select

case (rfcavity$, e_gun$, em_field$, crab_cavity$)
  select case (voltage_or_gradient)
  case (voltage$)
    field = ele%value(voltage$) * f_auto
  case (gradient$)
    field = ele%value(gradient$) * f_auto
  end select
end select

end function e_accel_field

