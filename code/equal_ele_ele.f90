!+
! Subroutine equal_ele_ele (ele1, ele2)
!
! Subroutine that is used to set one element equal to another. 
! This routine takes care of the pointers in ele1. 
!
! Note: This subroutine is called by the overloaded equal sign:
!		ele1 = ele2
!
! Input:
!   ele2 -- Ele_struct: Input element.
!
! Output:
!   ele1 -- Ele_struct: Output element.
!-

#include "CESR_platform.inc"

subroutine equal_ele_ele (ele1, ele2)

! do not use bmad_interface since "=" is overloaded with this routine

  use bmad_struct  
	
  implicit none
	
  type (ele_struct), intent(out) :: ele1
  type (ele_struct), intent(in) :: ele2

! deallocate ele1 pointers

  call deallocate_ele_pointers (ele1)

! set ele1 = ele2.
! if ele2 has allocated pointers then create new storage in ele1.
! if not allocated then init the pointers.

  ele1 = ele2
  call transfer_ele_pointers (ele1, ele2)

end subroutine
