!+
! Subroutine deallocate_ele_pointers (ele)
!
! Subroutine to deallocate the pointers in an element.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele -- ele_struct: Element with pointers.
!
! Output:
!   ele -- Ele_struct: Element with deallocated pointers.
!-

subroutine deallocate_ele_pointers (ele)

  use bmad

  implicit none

  type (ele_struct) ele
  integer ix

!

  if (ele%pointer_init == has_been_inited$) then
    if (associated (ele%wig_term)) deallocate (ele%wig_term)
    if (associated (ele%const)) deallocate (ele%const)
    if (associated (ele%taylor(1)%term)) deallocate &
             (ele%taylor(1)%term, ele%taylor(2)%term, ele%taylor(3)%term, &
              ele%taylor(4)%term, ele%taylor(5)%term, ele%taylor(6)%term)
    if (associated (ele%descrip)) deallocate (ele%descrip)
    if (associated (ele%a)) deallocate (ele%a, ele%b)
    call kill_gen_field (ele%gen_field)
  else
    nullify (ele%wig_term)
    nullify (ele%const)
    nullify (ele%taylor(1)%term, ele%taylor(2)%term, ele%taylor(3)%term, &
              ele%taylor(4)%term, ele%taylor(5)%term, ele%taylor(6)%term)
    nullify (ele%descrip)
    nullify (ele%a, ele%b)
    nullify (ele%gen_field)
    ele%pointer_init = has_been_inited$
  endif

end subroutine

