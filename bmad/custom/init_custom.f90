!+
! Subroutine init_custom (ele, err_flag)
!
! Prototype routine for initializing custom elements or elements that do custom
! calculations. Custom calculations are done if any one of the following
! ele_struct components is set to custom$:
!   ele%tracking_method
!   ele%mat6_calc_method
!   ele%field_calc
!   ele%aperture_type
!
! Input:
!   ele    -- Ele_struct: Element to init.
!
! Output:
!   ele      -- Ele_struct: Initalized element.
!   err_flag -- Logical: Set true if there is an error. False otherwise.
!+

subroutine init_custom (ele, err_flag)

use bmad_struct
use bmad_interface, except_dummy => init_custom

implicit none

type (ele_struct), target :: ele
logical err_flag

!

err_flag = .false.

end subroutine
