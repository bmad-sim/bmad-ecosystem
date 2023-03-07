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

elseif (master_parameter > num_ele_attrib$) then  ! Custom master parameter
  value = ele%custom(master_parameter-custom_attribute0$)

else  ! Normal
  value = ele%value(master_parameter)
endif

end function
