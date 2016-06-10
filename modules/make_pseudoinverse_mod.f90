module make_pseudoinverse_mod

contains 
!+
! Subroutine make pseudoinverse(A,Ap)
!-
subroutine make_pseudoinverse(A,Ap,svd_condition)
  use f95_lapack

  implicit none

  integer, parameter :: rp = 8

  real(rp) A(:,:)
  real(rp) Ap(:,:)
  real(rp), optional :: svd_condition

  real(rp) Atemp(size(A(:,1)),size(A(1,:)))
  integer n_bpms, n_correctors
  integer lwork, n_svals
  integer n_good_svals
  integer info
  integer i

  real(rp), allocatable :: S(:)
  real(rp), allocatable :: Ut(:,:)
  real(rp), allocatable :: V(:,:)
  real(rp), allocatable :: work(:)
  integer, allocatable :: iwork(:)

  n_bpms = SIZE(A(:,1))
  n_correctors = SIZE(A(1,:))
  lwork = n_correctors*(6+4*n_correctors)+n_bpms
  n_svals = MIN(n_bpms,n_correctors)

  allocate(S(n_svals))
  allocate(Ut(n_bpms,n_bpms))
  allocate(V(n_correctors,n_correctors))
  allocate(work(lwork))
  allocate(iwork(8*n_correctors))

  Atemp = A
  CALL dgesdd('A',n_bpms,n_correctors,Atemp,n_bpms,S,Ut,n_bpms,V,n_correctors,work,lwork,iwork,info) 
  V=transpose(V)
  Ut=transpose(Ut)

  ! Count good singular values
  n_good_svals = 0
  do i=1, n_svals
    if( S(i) .gt. 1.0E-8 ) THEN
      if(present(svd_condition)) then
        if( S(i)/S(1) .gt. svd_condition) then
          n_good_svals = n_good_svals + 1
        endif
      else
        n_good_svals = n_good_svals + 1
      endif
    endif
  enddo

  !Make pseudoinverse using only good singular values
  Ap(:,:) = 0
  do i=1,n_good_svals
    Ap(i,:) = Ut(i,:)/S(i)
  enddo
  Ap = matmul(V,Ap)

  deallocate(S)
  deallocate(Ut)
  deallocate(V)
  deallocate(work)
  deallocate(iwork)

end subroutine

end module
