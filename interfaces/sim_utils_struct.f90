!+
! sim_utils_struct
!-

#include "CESR_platform.inc"

module sim_utils_struct

use precision_def
use physical_constants

type spline_struct
  real(rp) x, y       ! data points
  real(rp) coef(0:3)  ! coefficients for cubic spline
end type

! A real_pointer_struct is just a pointer to a real number.
! This is used to construct arrays of real pointers.

type real_pointer_struct
  real(rp), pointer :: r
end type 

! This is to suppress the ranlib "has no symbols" message
integer, private :: private_dummy

end module
