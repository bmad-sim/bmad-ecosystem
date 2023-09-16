module moga_struct_mod

use precision_def

implicit none

type mag_struct
  character(1) type  ! 'c' for chromatic, 'h' for harmonic
  character(18) name
  character(5) property  ! element property to set
  real(rp) lb    !lower bound for constraint calculation
  real(rp) ub    !upper bound for constraint calculation
  real(rp) lir   !lower bound for initial population range
  real(rp) uir   !upper bound for initial population range
  real(rp) mutate_delta
end type


end module
