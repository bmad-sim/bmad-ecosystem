#include "CESR_platform.inc"

module equal_mod

  interface assignment (=)

    subroutine ele_equal_ele (ele1, ele2)
      use bmad_struct
      implicit none
      type (ele_struct), intent(inout) :: ele1
      type (ele_struct), intent(in) :: ele2
    end subroutine

    subroutine ele_vec_equal_ele_vec (ele1, ele2)
      use bmad_struct
      implicit none
      type (ele_struct), intent(inout) :: ele1(:)
      type (ele_struct), intent(in) :: ele2(:)
    end subroutine

    subroutine ring_equal_ring (ring1, ring2)
      use bmad_struct
      implicit none
      type (ring_struct), intent(inout) :: ring1
      type (ring_struct), intent(in) :: ring2
    end subroutine

    subroutine ring_vec_equal_ring_vec (ring1, ring2)
      use bmad_struct
      implicit none
      type (ring_struct), intent(inout) :: ring1(:)
      type (ring_struct), intent(in) :: ring2(:)
    end subroutine

!    elemental subroutine coord_equal_coord (coord1, coord2)
!      use bmad_struct
!      implicit none
!      type (coord_struct), intent(out) :: coord1
!      type (coord_struct), intent(in) :: coord2
!    end subroutine

  end interface

end module

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine ele_equal_ele (ele1, ele2)
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

subroutine ele_equal_ele (ele1, ele2)

! do not use bmad_interface since "=" is overloaded with this routine

  use bmad_struct

  implicit none
	
  type (ele_struct), intent(inout) :: ele1
  type (ele_struct), intent(in) :: ele2

! deallocate ele1 pointers

  call deallocate_ele_pointers (ele1)

! set ele1 = ele2.
! if ele2 has allocated pointers then create new storage in ele1.
! if not allocated then init the pointers.

  ele1 = ele2
  call transfer_ele_pointers (ele1, ele2)

end subroutine

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine ele_vec_equal_ele_vec (ele1, ele2)
!
! Subroutine that is used to set one element vector equal to another. 
! This routine takes care of the pointers in ele1. 
!
! Note: This subroutine is called by the overloaded equal sign:
!		ele1(:) = ele2(:)
!
! Input:
!   ele2(:) -- Ele_struct: Input element vector.
!
! Output:
!   ele1(:) -- Ele_struct: Output element vector.
!-

subroutine ele_vec_equal_ele_vec (ele1, ele2)

! do not use bmad_interface since "=" is overloaded with this routine

  use bmad_struct

  implicit none
	
  type (ele_struct), intent(inout) :: ele1(:)
  type (ele_struct), intent(in) :: ele2(:)

  integer i

! error check

  if (size(ele1) /= size(ele2)) then
    print *, 'ERROR IN ELE_VEC_EQUAL_ELE_VEC: ARRAY SIZES ARE NOT THE SAME!'
    call err_exit
  endif

! transfer

  do i = 1, size(ele1)
    ele1(i) = ele2(i)
  enddo

end subroutine

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine ring_equal_ring (ring1, ring2)
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

subroutine ring_equal_ring (ring1, ring2)

! do not use bmad_interface since "=" is overloaded with this routine

  use bmad_struct

  implicit none
	
  type (ring_struct), intent(inout) :: ring1
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

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine ring_vec_equal_ring_vec (ring1, ring2)
!
! Subroutine that is used to set one ring vector equal to another. 
! This routine takes care of the pointers in ring1. 
!
! Note: This subroutine is called by the overloaded equal sign:
!		ring1(:) = ring2(:)
!
! Input:
!   ring2(:) -- Ring_struct: Input ring vector.
!
! Output:
!   ring1(:) -- Ring_struct: Output ring vector.
!-

subroutine ring_vec_equal_ring_vec (ring1, ring2)

! do not use bmad_interface since "=" is overloaded with this routine

  use bmad_struct

  implicit none
	
  type (ring_struct), intent(inout) :: ring1(:)
  type (ring_struct), intent(in) :: ring2(:)

  integer i

! error check

  if (size(ring1) /= size(ring2)) then
    print *, 'ERROR IN RING_VEC_EQUAL_RING_VEC: ARRAY SIZES ARE NOT THE SAME!'
    call err_exit
  endif

! transfer

  do i = 1, size(ring1)
    ring1(i) = ring2(i)
  enddo

end subroutine

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine coord_equal_coord (coord1, coord2)
!
! Subroutine that is used to set one coord equal to another. 
!
! Note: This subroutine is called by the overloaded equal sign:
!		coord1 = coord2
!
! Input:
!   coord2 -- coord_struct: Input coord.
!
! Output:
!   coord1 -- coord_struct: Output coord.
!-

elemental subroutine coord_equal_coord (coord1, coord2)

  use bmad_struct

  implicit none
	
  type (coord_struct), intent(out) :: coord1
  type (coord_struct), intent(in) :: coord2

  coord1%vec = coord2%vec
 
end subroutine
