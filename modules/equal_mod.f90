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
  type (ele_struct) ele_save

  integer i

! If ele2 has not been properly initialized assume there is a problem somewhere

  if (ele2%pointer_init /= has_been_inited$) then
    print *, 'ERROR IN ELE_EQUAL_ELE: NO POINTER_INIT DONE FOR ELE2!'
    call err_exit
  endif

  if (ele1%pointer_init /= has_been_inited$) call deallocate_ele_pointers (ele1)

! save ele1 pointers and set ele1 = ele2.

  ele_save = ele1
  ele1 = ele2

! Transfer pointer info.
! When finished ele1's pointers will be pointing to a different memory
! location from ele2's so that the elements are truely separate.
! The exception is the %gen_field which is not transfered because it is
! part of PTC.

  if (associated(ele2%wig_term)) then
    if (associated (ele_save%wig_term)) then
      if (size(ele_save%wig_term) /= size(ele2%wig_term)) then
        deallocate (ele_save%wig_term)
        allocate (ele1%wig_term(size(ele2%wig_term)))
      endif
    else
      allocate (ele1%wig_term(size(ele2%wig_term)))
    endif
    ele1%wig_term = ele2%wig_term
  else
    if (associated (ele_save%wig_term)) deallocate (ele_save%wig_term)
  endif

  if (associated(ele2%const)) then
    if (associated (ele_save%const)) then
      if (size(ele_save%const) /= size(ele2%const)) then
        deallocate (ele_save%const)
        allocate (ele1%const(size(ele2%const)))
      endif
    else
      allocate (ele1%const(size(ele2%const)))
    endif
    ele1%const = ele2%const
  else
    if (associated (ele_save%const)) deallocate (ele_save%const)
  endif

  do i = 1, 6
    if (associated(ele2%taylor(i)%term)) then
      if (associated (ele_save%taylor(i)%term)) then
        if (size(ele_save%taylor(i)%term) /= size(ele2%taylor(i)%term)) then
          deallocate (ele_save%taylor(i)%term)
          allocate (ele1%taylor(i)%term(size(ele2%taylor(i)%term)))
        endif
      else
        allocate (ele1%taylor(i)%term(size(ele2%taylor(i)%term)))
      endif
      ele1%taylor(i)%term = ele2%taylor(i)%term
    else
      if (associated (ele_save%taylor(i)%term)) deallocate (ele_save%taylor(i)%term)
    endif
  enddo

  if (associated(ele2%a)) then
    if (.not. associated (ele_save%a)) then
      allocate (ele1%a(0:n_pole_maxx), ele1%b(0:n_pole_maxx))
    endif
    ele1%a = ele2%a
    ele1%b = ele2%b
  else
    if (associated (ele_save%a)) deallocate (ele_save%a, ele_save%b)
  endif

  if (associated(ele2%descrip)) then
    if (.not. associated (ele_save%descrip)) allocate (ele1%descrip)
    ele1%descrip = ele2%descrip
  else
    if (associated (ele_save%descrip)) deallocate (ele_save%descrip)
  endif

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
!               ele1(:) = ele2(:)
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

  interface
    subroutine ele_equal_ele (ele1, ele2)
      use bmad_struct
      implicit none
      type (ele_struct), intent(inout) :: ele1
      type (ele_struct), intent(in) :: ele2
    end subroutine
  end interface

! error check

  if (size(ele1) /= size(ele2)) then
    print *, 'ERROR IN ELE_VEC_EQUAL_ELE_VEC: ARRAY SIZES ARE NOT THE SAME!'
    call err_exit
  endif

! transfer

  do i = 1, size(ele1)
    call ele_equal_ele (ele1(i), ele2(i))
  enddo

end subroutine

!!----------------------------------------------------------------------
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
