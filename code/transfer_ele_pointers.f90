!+
! Subroutine transfer_ele_pointers (ele1, ele2, ele_save)
!
! This 
! Subroutine to transfer the information in the pointers from ele2 to ele1.
! ele_save is holds
! When finished ele1's pointers will be pointing to a different memory 
! location from ele2's so that the elements are truely separate.
!
! The exception is the %gen_field which is not transfered because it is
! part of PTC.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele2 -- Ele_struct: Input element holding the information.
!
! Output:
!   ele1 -- Ele_struct: Output element.
!-

subroutine transfer_ele_pointers (ele1, ele2)

  use bmad

  implicit none

  type (ele_struct), intent(in)  :: ele2
  type (ele_struct), intent(inout) :: ele1

  integer i

!

  if (ele2%pointer_init /= has_been_inited$) then
    print *, 'ERROR IN TRANSFER_ELE_POINTERS: NO POINTER_INIT DONE FOR ELE2!'
    call err_exit
  endif

  if (ele1%pointer_init /= has_been_inited$) call deallocate_ele_pointers (ele1)
  
  if (associated(ele2%wig_term)) then
    allocate (ele1%wig_term(size(ele2%wig_term)))
    ele1%wig_term = ele2%wig_term
  endif

  if (associated(ele2%taylor(1)%term)) then
    do i = 1, 6
      allocate (ele1%taylor(i)%term(size(ele2%taylor(i)%term)))
      ele1%taylor(i)%term = ele2%taylor(i)%term
    enddo
  endif

  if (associated(ele2%a)) then
    allocate (ele1%a(0:n_pole_maxx), ele1%b(0:n_pole_maxx))
    ele1%a = ele2%a
    ele1%b = ele2%b
  endif

  if (associated(ele2%descrip)) then
    allocate (ele1%descrip)
    ele1%descrip = ele2%descrip
  endif

!  if (associated(ele2%gen_field)) then
!    allocate (ele1%gen_field)
!    ele1%gen_field = ele2%gen_field
!  endif


end subroutine   
