module precision_def

#include "CESR_precision.inc"

integer, parameter :: sp = kind(1e0)
integer, parameter :: dp = selected_real_kind(2*precision(1e0_sp))

! This is to suppress the ranlib "has no symbols" message
integer, private :: private_dummy

type global_common_struct
  logical :: be_thread_safe = .false.                   ! Avoid thread unsafe practices?
end type

type (global_common_struct), save :: global_com

end module
