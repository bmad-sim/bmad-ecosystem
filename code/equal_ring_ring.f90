!+
! Subroutine equal_ring_ring (ring1, ring2)
!
! Subroutine that is used to set one ring equal to another. 
! This routine takes care of the pointers in ring1. 
!
! Note: This subroutine is called by the overloaded equal sign:
!		ring1 = ring2
!
! Input:
!   ring2 -- ring_struct: Input ring.
!
! Output:
!   ring1 -- ring_struct: Output ring.
!-


#include "CESR_platform.inc"

subroutine equal_ring_ring (ring1, ring2)

use bmad_struct  ! do not use bmad_interface since "=" is overloaded with this routine

  implicit none
	
  type (ring_struct), intent(out) :: ring1
  type (ring_struct), intent(in) :: ring2

  integer i
  
! deallocate ring1 pointers

  call deallocate_ele_pointers (ring1%ele_init)

  do i = 0, ubound(ring1%ele_, 1)
    call deallocate_ele_pointers (ring1%ele_(i))
  enddo

! set ring1 = ring2.
! if ring2 has allocated pointers then create new storate in ring1

  ring1 = ring2

  call transfer_ele_pointers (ring1%ele_init, ring2%ele_init)

  do i = 0, ring2%n_ele_max
    call transfer_ele_pointers (ring1%ele_(i), ring2%ele_(i))
  enddo


 
end subroutine
