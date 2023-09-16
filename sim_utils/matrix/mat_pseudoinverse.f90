!+
! Subroutine mat_pseudoinverse(A,Ap,svd_condition,print_err,ok)
!
! Makes Monroe-Penrose pseudoinverse of a matrix A.
!
! Input:
!   A(:,:)        -- Real(rp): Input matrix array
!   svd_condition -- Real(rp), optional: use only those singular values where
!                                        s(i)/s(1) > svd_condition for the pseudoinverse.
!   print_err     -- Logical, optional: If true, print info return from dgesdd if nonzero.
!   ok            -- Logical, optional: True if pseudoinverse calculated OK.  Else false.
!
! Output:
!   Ap(:,:)      -- Real(rp): Pseudoinverse of A.  dim(transpose(A))
!-

subroutine mat_pseudoinverse(A,Ap,svd_condition,print_err,ok)

use output_mod, except => mat_pseudoinverse
use f95_lapack

implicit none

real(rp) A(:,:)
real(rp) Ap(:,:)
real(rp), optional :: svd_condition
logical, optional :: print_err
logical, optional :: ok

real(rp) Atemp(size(A(:,1)),size(A(1,:)))
integer ndata, nvars
integer lwork, n_svals
integer n_good_svals
integer info
integer i

real(rp), allocatable :: S(:)
real(rp), allocatable :: Ut(:,:)
real(rp), allocatable :: V(:,:)
real(rp), allocatable :: work(:)
integer, allocatable :: iwork(:)

character(4) info_str
character(*), parameter :: r_name = 'mat_pseudoinverse'

if(present(ok)) ok = .true.

ndata = size(A(:,1))
nvars = size(A(1,:))
lwork = nvars*(6+4*nvars)+ndata
n_svals = min(ndata,nvars)

allocate(S(n_svals))
allocate(Ut(ndata,ndata))
allocate(V(nvars,nvars))
allocate(work(lwork))
allocate(iwork(8*nvars))

Atemp = A
CALL dgesdd('A',ndata,nvars,Atemp,ndata,S,Ut,ndata,V,nvars,work,lwork,iwork,info) 
if(info .ne. 0) then
  if (logic_option(.false., print_err)) then
    write(info_str,*) info
    call out_io (s_error$, r_name, 'dgesdd failed. info = '//info_str)
  endif
  if (present(ok)) ok = .false.
endif
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
