!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!+
! Subroutine deallocate_ele_array_pointers (eles)
!
! Routine to deallocate the pointers of all the elements in an 
! element array and the array itself.
!
! Input:
!   eles(:) -- Ele_struct, pointer: Array of elements.
!
! Output:
!   eles(:) -- Ele_struct, pointer: Deallocated array.
!-

subroutine deallocate_ele_array_pointers (eles)

use bmad_routine_interface, dummy => deallocate_ele_array_pointers

implicit none

type (ele_struct), pointer :: eles(:)
integer i

!

if (.not. associated(eles)) return

do i = lbound(eles, 1), ubound(eles, 1)
  call deallocate_ele_pointers (eles(i))
enddo

deallocate (eles)

end subroutine deallocate_ele_array_pointers

