!+
! Subroutine elements_locator (key, ring, indx)
!
! Subroutine to locate all the elements of a certain kind in a ring. 
! Note: super_slave elements are not included in the list since super_slaves
! cannot have their attributes changed directly.
!
! Modules Needed:
!   use bmad
!
! Input:
!   key  -- Integer: Element key. For example quadrupole$, etc.
!   ring -- Ring_struct: Ring to search through.
!
! Output:
!   indx(:) -- Integer, pointer: Array of indexes of the elements found.
!              Note: This routine does not try to deallocate indx.
!               It is up to you to deallocate indx if needed.
!-

#include "CESR_platform.inc"

subroutine elements_locator (key, ring, indx)

  use bmad_struct
  use bmad_interface, except => elements_locator
  
  implicit none

  type (ring_struct) ring

  integer key
  integer, pointer :: indx(:)
  integer i, ix

!

  ix = count(ring%ele_(:)%key == key .and. &
                                  ring%ele_(:)%control_type /= super_slave$)

  allocate(indx(ix))

  ix = 0
  do i = 0, ring%n_ele_max
    if (ring%ele_(i)%key == key .and. &
                              ring%ele_(i)%control_type /= super_slave$) then
      ix = ix + 1
      indx(ix) = i
    endif
  enddo

end subroutine
