!+
! Function master_parameter_value(master_parameter, ele) result (value)
!
! Routine to calculate the "master parameter" value for a fieldmap.
!
! Input:
!   master_parameter  -- integer: Index of the master parameter.
!   ele               -- ele_struct: Element containing the fieldmap.
!
! Output:
!   value             -- real(rp): Value of the master parameter.
!-

function master_parameter_value (master_parameter, ele) result (value)

use bmad, dummy => master_parameter_value

implicit none

type (ele_struct) ele
real(rp) value
integer master_parameter

! master_parameter < 1 means it does not exist so the scale value is 1.

if (master_parameter < 1) then
  value = 1
  return
endif

select case (attribute_name(ele, master_parameter))
case ('VOLTAGE');   value = ele%value(voltage$) + ele%value(voltage_err$)
case ('GRADIENT');  value = ele%value(gradient$) + ele%value(gradient_err$)
case ('PHI0');      value = ele%value(phi0$) + ele%value(phi0_err$)
case ('G');         value = ele%value(g$) + ele%value(dg$)
case ('B_FIELD');   value = ele%value(b_field$) + ele%value(db_field$)
case default;       value = ele%value(master_parameter)
end select

end function
