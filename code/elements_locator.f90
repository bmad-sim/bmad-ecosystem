!+
! Subroutine elements_locator (key, ring, indx)
!
! Subroutine to locate all the elements of a certain kind in a ring. 
! Note: super_slave elements are not included in the list.
!
! Modules Needed:
!   use bmad
!
! Input:
!   key  -- Integer: Element key. For example quadrupole$, etc.
!   ring -- Ring_struct: Ring to search through.
!
! Output:
!   indx(:) -- Integer, allocatable: Array of indexes of the elements found..
!-

#include "CESR_platform.inc"

subroutine elements_locator (key, ring, indx)

  use bmad
  
  implicit none

  type (ring_struct) ring

  integer key
  integer, allocatable :: indx(:)
  integer i, ix

!

  ix = count(ring%ele_(:)%key == key .and. &
                                  ring%ele_(:)%control_type /= super_slave$)

  if (allocated(indx)) deallocate(indx)
  allocate(indx(ix))

  ix = 0
  do i = 1, ring%n_ele_max
    if (ring%ele_(i)%key == key .and. &
                              ring%ele_(i)%control_type /= super_slave$) then
      ix = ix + 1
      indx(ix) = i
    endif
  enddo

end subroutine
