module bmad_base_mod

use precision_def
use output_mod

! common flags
! status structure

type bmad_status_struct
  logical :: ok             = .true.   ! Error flag
  logical :: type_out       = .true.   ! Print error messages?
  logical :: exit_on_error  = .true.   ! Exit program on error?
end type

type (bmad_status_struct), save :: bmad_status

end module
