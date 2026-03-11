!+
! Module private_pointer_mod
!
! This module is used to pass variables to nested functions to avoid the dreaded GCC executable stack issue.
!-

module private_pointer_mod

use bmad_struct

!

type (ele_struct), pointer :: ele_ptr
type (lat_param_struct), pointer :: param_ptr
type (coord_struct), pointer :: start_ptr


end module
