!+
! sim_utils_struct
!-

module sim_utils_struct

use precision_def
use physical_constants

! A real_pointer_struct is just a pointer to a real number.
! This is used to construct arrays of real pointers.

type real_pointer_struct
  real(rp), pointer :: r => null()
end type 

type all_pointer_struct
  real(rp), pointer :: r => null()
  integer, pointer :: i => null()
  logical, pointer :: l => null()
end type 

! This is to suppress the ranlib "has no symbols" message
integer, private :: private_dummy

end module
