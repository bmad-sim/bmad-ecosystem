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

type crm_struct
  logical stale
  type(mag_struct), pointer :: l_mags(:)
  type(mag_struct), pointer :: c_mags(:)
  real(rp) set_chrom_x
  real(rp) set_chrom_y
  real(rp), allocatable :: ApC(:)
  real(rp), allocatable :: Q1(:,:)
  real(rp), allocatable :: Q1t(:,:)
end type


end module
