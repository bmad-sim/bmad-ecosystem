!+
! Subroutine tao_data_coupling_init (u)
!
! Routine to initialize the coupling structure for a lattice branch.
! This routine is called by tao_lattic_calc and is not meant for general use.
!
! Input:
!   branch -- branch_struct: New lattice branch.
!-

subroutine tao_data_coupling_init (branch)

use tao_struct

implicit none

type (branch_struct) branch
integer m

! 

m = branch%n_ele_max
if (.not. allocated(scratch%cc)) allocate (scratch%cc(0:m))
if (ubound(scratch%cc, 1) < m) then
  deallocate(scratch%cc)
  allocate(scratch%cc(0:m))
endif

scratch%cc%coupling_calc_done = .false.
scratch%cc%amp_calc_done = .false.

end subroutine tao_data_coupling_init

