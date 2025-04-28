!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine transfer_ele (ele1, ele2, nullify_pointers)
!
! Subroutine to set ele2 = ele1. 
! This is a plain transfer of information not using the overloaded equal operator.
! The result is that ele2's pointers will point to the same memory as ele1's.
!
! NOTE: Do not use this routine unless you know what you are doing!
!
! Input:
!   ele1             -- Ele_struct:
!   nullify_pointers -- Logical, optional: If present and True then nullify the 
!                         pointers in ele2 except for the ele2%lat and ele2%lord pointers. 
!                         This gives a "bare bones" copy where one does not have to 
!                         worry about deallocating allocated structure components later.
!
! Output:
!   ele2 -- Ele_struct:
!-

subroutine transfer_ele (ele1, ele2, nullify_pointers)

! Important! The use statement here is constructed to  avoid the use of 
! the overloaded equal sign for ele_structs in bmad_routine_interface.

use bmad_routine_interface, only: ele_struct, deallocate_ele_pointers, logic_option

implicit none

type (ele_struct), target :: ele1
type (ele_struct) :: ele2
logical, optional :: nullify_pointers

!

ele2 = ele1

if (logic_option (.false., nullify_pointers)) then
  call deallocate_ele_pointers (ele2, .true.)
  ele2%branch => ele1%branch  ! Reinstate
  ele2%lord   => ele1%lord    ! Reinstate
endif

end subroutine transfer_ele

