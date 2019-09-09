!+
! Subroutine shift_eles (branch, ix0, ix1, n_shift)
!
! Routine to shift the lattice elements in a branch:
!   branch%ele(ix0+n_shift:ix1+n_shift) = branch%ele(ix0:ix1). 
! This is a plain transfer of information not using the overloaded equal operator.
!
! NOTE: Do not use this routine unless you know what you are doing!
!
! Input:
!   branch    -- branch_struct: Branch containing elements.
!   ix0       -- integer: Lower bound index for the shift region.
!   ix1       -- integer: Upper bound index for the shift region.
!   n_shift   -- integer: Number of slots to shift.
!
! Output:
!   branch    -- branch_struct: Branch with shifted elements.
!-

subroutine shift_eles (branch, ix0, ix1, n_shift)

use bmad_routine_interface, dummy => shift_eles

implicit none

type (branch_struct) :: branch
integer ix0, ix1, n_shift, ie

!

branch%ele(ix0+n_shift:ix1+n_shift) = branch%ele(ix0:ix1)

do ie = ix0+n_shift, ix1+n_shift
  branch%ele(ie)%ix_ele = ie
enddo

end subroutine shift_eles

