#include "CESR_platform.inc"

module equal_mod

  use bmad_utils_mod

  interface assignment (=)
    module procedure ele_equal_ele
    module procedure ele_vec_equal_ele_vec
    module procedure ring_equal_ring 
    module procedure ring_vec_equal_ring_vec 
!    module procedure coord_equal_coord 
  end interface

contains

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

  implicit none
	
  type (ele_struct), intent(inout) :: ele1
  type (ele_struct), intent(in) :: ele2
  type (ele_struct) ele_save

  integer i

! save ele1 pointers and set ele1 = ele2.

  call transfer_ele (ele1, ele_save)
  call transfer_ele (ele2, ele1)

! Transfer pointer info.
! When finished ele1's pointers will be pointing to a different memory
! location from ele2's so that the elements are truely separate.
! The exception is the %gen_field which is not transfered because it is
! part of PTC.

  if (associated(ele2%wig_term)) then
    if (associated (ele_save%wig_term)) then
      if (size(ele_save%wig_term) == size(ele2%wig_term)) then
        ele1%wig_term => ele_save%wig_term
      else
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
      if (size(ele_save%const) == size(ele2%const)) then
        ele1%const => ele_save%const
      else
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

  if (associated(ele2%r)) then
    if (associated (ele_save%r)) then
      if (size(ele_save%r) == size(ele2%r)) then
        ele1%r => ele_save%r
      else
        deallocate (ele_save%r)
        allocate (ele1%r(size(ele2%r)))
      endif
    else
      allocate (ele1%r(size(ele2%r)))
    endif
    ele1%r = ele2%r
  else
    if (associated (ele_save%r)) deallocate (ele_save%r)
  endif

  do i = 1, 6
    if (associated(ele2%taylor(i)%term)) then
      if (associated (ele_save%taylor(i)%term)) then
        if (size(ele_save%taylor(i)%term) == size(ele2%taylor(i)%term)) then
          ele1%taylor(i)%term => ele_save%taylor(i)%term
        else
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
    if (associated (ele_save%a)) then
      ele1%a => ele_save%a 
      ele1%b => ele_save%b 
    else
      allocate (ele1%a(0:n_pole_maxx), ele1%b(0:n_pole_maxx))
    endif
    ele1%a = ele2%a
    ele1%b = ele2%b
  else
    if (associated (ele_save%a)) deallocate (ele_save%a, ele_save%b)
  endif

  if (associated(ele2%descrip)) then
    if (associated (ele_save%descrip)) then
      ele1%descrip => ele_save%descrip
    else
      allocate (ele1%descrip)
    endif
    ele1%descrip = ele2%descrip
  else
    if (associated (ele_save%descrip)) deallocate (ele_save%descrip)
  endif

  if (associated(ele2%wake%sr_file)) then
    if (associated (ele_save%wake%sr_file)) then
      ele1%wake%sr_file => ele_save%wake%sr_file
    else
      allocate (ele1%wake%sr_file)
    endif
    ele1%wake%sr_file = ele2%wake%sr_file
  else
    if (associated (ele_save%wake%sr_file)) deallocate (ele_save%wake%sr_file)
  endif

  if (associated(ele2%wake%lr_file)) then
    if (associated (ele_save%wake%lr_file)) then
      ele1%wake%lr_file => ele_save%wake%lr_file
    else
      allocate (ele1%wake%lr_file)
    endif
    ele1%wake%lr_file = ele2%wake%lr_file
  else
    if (associated (ele_save%wake%lr_file)) deallocate (ele_save%wake%lr_file)
  endif

  if (associated(ele2%wake%sr)) then
    if (associated (ele_save%wake%sr)) then
      if (size(ele_save%wake%sr) == size(ele2%wake%sr)) then
        ele1%wake%sr => ele_save%wake%sr
      else
        deallocate (ele_save%wake%sr)
        allocate (ele1%wake%sr(0:ubound(ele2%wake%sr,1)))
      endif
    else
      allocate (ele1%wake%sr(0:ubound(ele2%wake%sr,1)))
    endif
    ele1%wake%sr = ele2%wake%sr
  else
    if (associated (ele_save%wake%sr)) deallocate (ele_save%wake%sr)
  endif

  if (associated(ele2%wake%lr)) then
    if (associated (ele_save%wake%lr)) then
      if (size(ele_save%wake%lr) == size(ele2%wake%lr)) then
        ele1%wake%lr => ele_save%wake%lr
      else
        deallocate (ele_save%wake%lr)
        allocate (ele1%wake%lr(0:ubound(ele2%wake%lr,1)))
      endif
    else
      allocate (ele1%wake%lr(0:ubound(ele2%wake%lr,1)))
    endif
    ele1%wake%lr = ele2%wake%lr
  else
    if (associated (ele_save%wake%lr)) deallocate (ele_save%wake%lr)
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
    call ele_equal_ele (ele1(i), ele2(i))
  enddo

end subroutine

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine ring_equal_ring (ring_out, ring_in)
!
! Subroutine that is used to set one ring equal to another. 
! This routine takes care of the pointers in ring_in. 
!
! Note: This subroutine is called by the overloaded equal sign:
!		ring_out = ring_in
!
! Input:
!   ring_in -- ring_struct: Input ring.
!
! Output:
!   ring_out -- ring_struct: Output ring.
!-

subroutine ring_equal_ring (ring_out, ring_in)

  implicit none

  type (ring_struct), intent(inout) :: ring_out
  type (ring_struct), intent(in) :: ring_in
  type (ele_struct), allocatable :: temp_ele(:)

  integer i, n, n_old, n_wanted

! If ring_in has not been properly initialized then assume there is a problem somewhere

  if (.not. associated (ring_in%ele_)) then
    print *, 'ERROR IN RING_EQUAL_RING: NO RING%ELE_ ON RHS NOT ASSOCIATED!'
    call err_exit
  endif

! resize %ele_ array if needed

  if (ring_out%n_ele_maxx < ring_in%n_ele_maxx) &
                                  call allocate_ring_ele_(ring_out, ring_in%n_ele_maxx)
  
  do i = 0, ring_in%n_ele_max
    ring_out%ele_(i) = ring_in%ele_(i)
  enddo
  ring_out%ele_init = ring_in%ele_init

! handle ring%control_ array

  n = size(ring_in%control_)
  if (associated (ring_in%control_)) then
    if (associated(ring_out%control_)) then
      if (size (ring_in%control_) < size (ring_out%control_)) then
        ring_out%control_(1:n) = ring_in%control_(1:n)
      else
        deallocate (ring_out%control_)
        allocate (ring_out%control_(n))
        ring_out%control_ = ring_in%control_
      endif
    else
      allocate (ring_out%control_(n))
      ring_out%control_ = ring_in%control_
    endif
  else
    if (associated(ring_out%control_)) deallocate (ring_out%control_)
  endif

! handle ring%ic_ array

  n = size(ring_in%ic_)
  if (associated (ring_in%ic_)) then
    if (associated(ring_out%ic_)) then
      if (size (ring_in%ic_) < size (ring_out%ic_)) then
        ring_out%ic_(1:n) = ring_in%ic_(1:n)
      else
        deallocate (ring_out%ic_)
        allocate (ring_out%ic_(n))
        ring_out%ic_ = ring_in%ic_
      endif
    else
      allocate (ring_out%ic_(n))
      ring_out%ic_ = ring_in%ic_
    endif
  else
    if (associated(ring_out%ic_)) deallocate (ring_out%ic_)
  endif

! non-pointer transfer

  call transfer_ring_parameters (ring_in, ring_out)

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
    call ring_equal_ring (ring1(i), ring2(i))
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

  implicit none
	
  type (coord_struct), intent(out) :: coord1
  type (coord_struct), intent(in) :: coord2

  coord1%vec = coord2%vec
 
end subroutine

end module

