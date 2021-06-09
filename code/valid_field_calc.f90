!+
! Function valid_field_calc (ele, field_calc) result (is_valid)
!
! Routine to return whether a given field_calc method is valid for a given element.
!
! Input:
!   ele         -- ele_struct: Lattice element.
!   field_calc  -- integer: bmad_standard$, etc.
!
! Output:
!   is_valid    -- logical: True if a valid method. False otherwise.
!-

function valid_field_calc (ele, field_calc) result (is_valid)

use bmad_struct

implicit none

type (ele_struct) ele
integer field_calc
logical is_valid

! 

is_valid = .false.

select case (ele%key)

case (group$, overlay$, girder$)
  select case (field_calc)
  case (no_field$)
    is_valid = .true.
  end select

case (null_ele$)
  ! No valid calc

case (wiggler$, undulator$)
  select case (field_calc)
  case (bmad_standard$, fieldmap$, planar_model$, refer_to_lords$, helical_model$, custom$)
    is_valid = .true.
  end select

case (solenoid$)
  select case (field_calc)
  case (bmad_standard$, fieldmap$, refer_to_lords$, custom$, soft_edge$)
    is_valid = .true.
  end select

case default
  select case (field_calc)
  case (bmad_standard$, fieldmap$, refer_to_lords$, custom$)
    is_valid = .true.
  end select
end select

end function valid_field_calc 

