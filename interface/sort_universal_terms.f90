!+
! Subroutine sort_universal_terms (ut_in, ut_sorted)
!
! Subroutine to sort the taylor terms from "lowest" to "highest".
! This subroutine is needed since what comes out of PTC is not sorted.
!
! The number associated with a taylor_term that is used for the sort is:
!     number = sum(exp(i))*10^6 + exp(6)*10^5 + ... + exp(1)*10^0
! Where exp(1) is the exponent for x, exp(2) is the exponent for P_x, etc.
!
! Note: ut_sorted needs to have been initialized.
! Note: ut_sorted cannot be ut_in. That is it is not legal to write:
!           call sort_universal_terms (this_ut, this_ut)
!
! Modules needed:
!   use s_tracking
!
! Input:
!   ut_in     -- Universal_taylor: Unsorted taylor series.
!
! Output:
!   ut_sorted -- Universal_taylor: Sorted taylor series.
!-

subroutine sort_universal_terms (ut_in, ut_sorted)

  use s_tracking
  use nr

  implicit none

  type (universal_taylor), intent(in)  :: ut_in
  type (universal_taylor) :: ut_sorted

  integer, allocatable :: ix_(:), ord_(:)
  integer i, j, n, nv, expn(6)

! init

  n = ut_in%n
  nv = ut_in%nv

  if (nv /= 6) then
    print *, 'ERROR IN SORT_UNIVERSAL_TERMS: I AM NOT SET UP FOR NV /= 6'
    call err_exit
  endif

  if (associated(ut_sorted%n)) &
              deallocate(ut_sorted%n, ut_sorted%nv, ut_sorted%c, ut_sorted%j)
  allocate(ut_sorted%n, ut_sorted%nv, ut_sorted%c(n), ut_sorted%j(n,nv), &
                                                              ix_(n), ord_(n))

  ut_sorted%n = n
  ut_sorted%nv = nv

!

  do i = 1, n
    expn = ut_in%j(i,:)
    ord_(i) = sum(expn)*10**6 + expn(6)*10**5 + expn(5)*10**4 + &
                expn(4)*10**3 + expn(3)*10**2 + expn(2)*10**1 + expn(1)
  enddo

  call indexx (ord_, ix_)

  do i = 1, n
    ut_sorted%c(i)= ut_in%c(ix_(i))
    ut_sorted%j(i,:)= ut_in%j(ix_(i),:)
  enddo

  deallocate(ord_, ix_)

end subroutine
