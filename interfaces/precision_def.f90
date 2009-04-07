module precision_def

#include "CESR_precision.inc"

integer, parameter :: sp = kind(1e0)
integer, parameter :: dp = selected_real_kind(2*precision(1e0_sp))

! This is to suppress the ranlib "has no symbols" message
integer, private :: private_dummy


end module
