!+
! Subroutine type_map (y)
!
! Subroutine to type the transfer maps of a real_8 array.
!
! Modules needed:
!   use accelerator
!
! Input:
!   y(:)  -- Real_8: 
!-

subroutine type_map (y)

  use accelerator

  implicit none

  type (real_8), intent(in) :: y(:)
  type (universal_taylor) ut

  integer :: i, j, k

!

  do i = 1, size(y)
    ut = 0
    ut = y(i)%t
    type *, '!-----------------'
    type *, '! Term:', i
    type *, 'Order            Coef    Exponents'
    do j = 1, ut%n
      type '(i6, f18.14, 20i3)', sum(ut%j(j,:)), ut%c(j), &
                                          (ut%j(j,k), k = 1, ut%nv)
    enddo
    ut = -1
  enddo

end subroutine
