!+
! subroutine sort_taylor_terms (taylor_in, taylor_sorted)
!
! Subroutine to sort the taylor terms from "lowest" to "highest".
! This subroutine is needed since what comes out of PTC is not sorted.
!
! The number associated with a taylor_term that is used for the sort is:
!     number = sum(exp(i))*10^6 + exp(6)*10^5 + ... + exp(1)*10^0
! Where exp(1) is the exponent for x, exp(2) is the exponent for P_x, etc.
!
! Note: taylor_sorted needs to have been initialized.
! Note: taylor_sorted cannot be taylor_in. That is it is not legal to write:
!           call sort_taylor_terms (this_taylor, this_taylor)
!
! Modules needed:
!   use bmad
!
! Input:
!   taylor_in     -- Taylor_struct: Unsorted taylor series.
!
! Output:
!   taylor_sorted -- Taylor_struct: Sorted taylor series.
!-

subroutine sort_taylor_terms (taylor_in, taylor_sorted)

  use bmad
  use nr

  implicit none

  type (taylor_struct), intent(in)  :: taylor_in
  type (taylor_struct) :: taylor_sorted
  type (taylor_term_struct), allocatable :: tt(:)
  
  integer, allocatable :: ord_(:), ix_(:)

  integer i, j, n, expn(6)

! init

  n = size(taylor_in%term)
  if (associated(taylor_sorted%term)) deallocate(taylor_sorted%term)
  allocate(taylor_sorted%term(n), ix_(n), ord_(n), tt(n))

!

  tt = taylor_in%term

  do i = 1, n
    expn = tt(i)%exp
    ord_(i) = sum(expn)*10**6 + expn(6)*10**5 + expn(5)*10**4 + &
                expn(4)*10**3 + expn(3)*10**2 + expn(2)*10**1 + expn(1)
  enddo

  call indexx (ord_, ix_)

  do i = 1, n
    taylor_sorted%term(i)= tt(ix_(i))
  enddo

  deallocate(ord_, ix_, tt)

end subroutine
