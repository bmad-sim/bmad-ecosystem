#include "CESR_platform.inc"

module equal_mod

  use bmad_utils_mod

  interface assignment (=)
    module procedure ele_equal_ele
    module procedure ele_vec_equal_ele_vec
    module procedure ring_equal_ring 
    module procedure ring_vec_equal_ring_vec 
!    module procedure coord_equal_coord
    module procedure mp_slice_equal_mp_slice
    module procedure mp_bunch_equal_mp_bunch
    module procedure mp_beam_equal_mp_beam
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

  use tpsalie_analysis 
  use multipole_mod

  implicit none
	
  type (ele_struct), intent(inout) :: ele1
  type (ele_struct), intent(in) :: ele2
  type (ele_struct) ele_save

  integer i, n_sr, n_lr

! 1) Save ele1 pointers in ele_save
! 2) Set ele1 = ele2.

  call transfer_ele (ele1, ele_save)
  call transfer_ele (ele2, ele1)
  ele1%ix_ele = ele_save%ix_ele     ! this should not change.

! Transfer pointer info.
! When finished ele1's pointers will be pointing to a different memory
! location from ele2's so that the elements are truely separate.

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
      if (all(lbound(ele_save%r) == lbound(ele2%r)) .and. &
          all(ubound(ele_save%r) == ubound(ele2%r)) ) then
        ele1%r => ele_save%r
      else
        deallocate (ele_save%r)
        allocate (ele1%r(lbound(ele2%r,1):ubound(ele2%r,1), &
                         lbound(ele2%r,2):ubound(ele2%r,2)))
      endif
    else
      allocate (ele1%r(lbound(ele2%r,1):ubound(ele2%r,1), &
                       lbound(ele2%r,2):ubound(ele2%r,2)))
    endif
    ele1%r = ele2%r
  else
    if (associated (ele_save%r)) deallocate (ele_save%r)
  endif

  do i = 1, 6
    ele1%taylor(i)%term => ele_save%taylor(i)%term ! reinstate
    ele1%taylor(i) = ele2%taylor(i)      ! use overloaded taylor_equal_taylor
  enddo

  if (associated(ele2%a)) then
    ele1%a => ele_save%a   ! reinstate
    ele1%b => ele_save%b   ! reinstate
    call multipole_init (ele1)
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

  n_sr = sr_wake_array_ubound(ele2); n_lr = lr_wake_array_size(ele2)
  ele1%wake => ele_save%wake  ! reinstate
  call init_wake (ele1%wake, n_sr, n_lr)
  if (n_sr /= -1)  ele1%wake%sr = ele2%wake%sr
  if (n_lr /= 0)   ele1%wake%lr = ele2%wake%lr
  if (associated(ele1%wake)) then
    ele1%wake%sr_file = ele2%wake%sr_file
    ele1%wake%lr_file = ele2%wake%lr_file
  endif

! gen_fields are hard because it involves pointers in PTC.
! just kill the gen_field in ele1 for now.

  if (associated(ele_save%gen_field)) call kill_gen_field (ele_save%gen_field)
  if (associated(ele1%gen_field)) nullify (ele1%gen_field)

!  if (associated(ele2%gen_field)) then
!    if (associated (ele_save%gen_field)) then
!      ele1%gen_field => ele_save%gen_field
!      call kill (ele1%gen_field)
!    else
!      allocate (ele1%gen_field)
!    endif
!    ele1%gen_field = ele2%gen_field
!  else
!    if (associated (ele_save%gen_field)) &
!                              call kill_gen_field (ele_save%gen_field)
!  endif

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

  integer i, n

! If ring_in has not been properly initialized then assume there is 
! a problem somewhere

  if (.not. associated (ring_in%ele_)) then
    print *, 'ERROR IN RING_EQUAL_RING: RING%ELE_(:) ON RHS NOT ASSOCIATED!'
    call err_exit
  endif

! resize %ele_ array if needed

  if (ubound(ring_out%ele_, 1) < ubound(ring_in%ele_, 1)) &
               call allocate_ring_ele_(ring_out, ubound(ring_in%ele_, 1))
  
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

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine mp_slice_equal_mp_slice (slice1, slice2)
!
! Subroutine to set one macroparticle slice equal to another taking care of
! pointers so that they don't all point to the same place.
!
! Note: This subroutine is called by the overloaded equal sign:
!		slice1 = slice2
!
! Input: 
!  slice2 -- macro_slice_struct: Input slice
!
! Output
!  slice1 -- macro_slice_struct: Output slice
!
!-

subroutine mp_slice_equal_mp_slice (slice1, slice2)

  implicit none

  type (macro_slice_struct), intent(inout) :: slice1
  type (macro_slice_struct), intent(in)    :: slice2


  if (associated(slice1%macro)) deallocate(slice1%macro)
  allocate(slice1%macro(size(slice2%macro)))

  slice1%macro(:)  = slice2%macro(:)
  slice1%charge    = slice2%charge

end subroutine mp_slice_equal_mp_slice

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine mp_bunch_equal_mp_bunch (bunch1, bunch2)
!
! Subroutine to set one macroparticle bunch equal to another taking care of
! pointers so that they don't all point to the same place.
!
! Note: This subroutine is called by the overloaded equal sign:
!		bunch1 = bunch2
!
! Input: 
!  bunch2 -- macro_bunch_struct: Input bunch
!
! Output
!  bunch1 -- macro_bunch_struct: Output bunch
!
!-

subroutine mp_bunch_equal_mp_bunch (bunch1, bunch2)

  implicit none

  type (macro_bunch_struct), intent(inout) :: bunch1
  type (macro_bunch_struct), intent(in)    :: bunch2

  integer i

  if (associated(bunch1%slice)) then
    do i = 1, size(bunch1%slice)
      if (associated(bunch1%slice(i)%macro)) deallocate(bunch1%slice(i)%macro)
    enddo  
    deallocate(bunch1%slice)
  endif
  allocate(bunch1%slice(size(bunch2%slice)))

  do i = 1, size(bunch2%slice)
    call mp_slice_equal_mp_slice(bunch1%slice(i), bunch2%slice(i))
  enddo
  bunch1%charge    = bunch2%charge
  bunch1%s_center  = bunch2%s_center

end subroutine mp_bunch_equal_mp_bunch

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine mp_beam_equal_mp_beam (beam1, beam2)
!
! Subroutine to set one macroparticle beam equal to another taking care of
! pointers so that they don't all point to the same place.
!
! Note: This subroutine is called by the overloaded equal sign:
!		beam1 = beam2
!
! Input: 
!  beam2 -- macro_beam_struct: Input beam
!
! Output
!  beam1 -- macro_beam_struct: Output beam
!
!-

subroutine mp_beam_equal_mp_beam (beam1, beam2)

  implicit none

  type (macro_beam_struct), intent(inout) :: beam1
  type (macro_beam_struct), intent(in)    :: beam2

  integer i, j

  if (associated(beam1%bunch)) then
    do i = 1, size(beam1%bunch)
      if (associated(beam1%bunch(i)%slice)) then
	do j = 1, size(beam1%bunch(i)%slice)
       	  if (associated(beam1%bunch(i)%slice(j)%macro)) &
	                   deallocate(beam1%bunch(i)%slice(j)%macro)
        enddo
	deallocate(beam1%bunch(i)%slice)
      endif
    enddo
    deallocate(beam1%bunch)
  endif

  allocate(beam1%bunch(size(beam2%bunch)))
  
  do i = 1, size(beam2%bunch)
    call mp_bunch_equal_mp_bunch (beam1%bunch(i), beam2%bunch(i))
  enddo

end subroutine mp_beam_equal_mp_beam

end module

