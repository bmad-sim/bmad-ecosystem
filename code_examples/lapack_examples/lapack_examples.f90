program lapack_examples
!program lapack_examples
!
! FIXME
! usage:
!	svd_test M N ntests
!            (integers)
!   
!   Tests the timing and error of ntests * (MxN) random matrices
!   for
!   LA1: DGESVD_F95( A, S, U, VT) from LAPACK95 
!   LA2: DGESDD_F95( A, S, U, VT) from LAPACK95
!   NR : svdcmp(A, S, V)          from Numerical Recipes
!
!output:
!   svd_report.txt
!


!use la_precision
USE LA_PRECISION, ONLY: WP => DP
use f95_lapack

use bmad

implicit none


real(rp), allocatable :: A(:,:)      
real(rp), allocatable :: B(:,:)      
real(rp), allocatable :: C1(:,:)      
real(rp), allocatable :: C2(:,:)     
real(rp), allocatable :: X(:), Y(:), X1(:)
real(rp), allocatable :: error(:)      
real(rp) time, mm_time, gemm_time


character(10) arg1, arg2
integer i, j, N
integer test, ntests

!

ntests = 100
N = 200


! Get matrix size and ntests from command line
call getarg(1, arg1)
call getarg(2, arg2)
read(arg1, '(i5)') N
read(arg2, '(i5)') ntests


allocate(A(1:N, 1:N))
allocate(B(1:N, 1:N))
allocate(C1(1:N, 1:N))
allocate(C2(1:N, 1:N))
allocate(error(1:ntests))
allocate(X(N), X1(N), Y(N))


mm_time = 0
gemm_time = 0
do, test = 1, ntests

! Make random NxN matrices A and B
do i=1, N
	do j=1, N
		call ran_gauss( A(i,j) )
	enddo
enddo
do i=1, N
	do j=1, N
		call ran_gauss( B(i,j) )
	enddo
enddo

! Get matmul time
call run_timer ('START')
C1 = matmul(A,B)
call run_timer ('STOP', time)
mm_time = mm_time + time

! Get dgemm time
call run_timer ('START')
call dgemm('N', 'N', N, N, N, 1.0_rp, A, N, B, N, 0, C2, N) 
call run_timer ('STOP', time)
gemm_time = gemm_time + time


!Error calc
error(test) = maxval(C1-C2)

enddo

open(1, file = 'matmul_report.txt')
write(1, '(a, i8)' )          'Matrix size (N X N)       N = ', N
write(1, '(a, i8    )' )      'Number of random matrices   = ', ntests
write(1, '(a, 2es13.5, a )' ) 'total time for (matmul, dgemm) = ', mm_time, gemm_time, ' seconds'
close(1)

print '(i8, i8, 3es13.3)', N, ntests, mm_time, gemm_time, maxval(error)
                 
!Linear Solve

call run_timer ('START')
do test=1, ntests
  do i=1, N
	  do j=1, N
	   call ran_gauss( A(i,j) )
	  enddo
  enddo

  do i=1, N
   call ran_gauss( X(i) )
  enddo

  Y = matmul(A, X)

  call LA_GESVX(A, Y, X1)

  error(test) =  maxval(abs(X1-X))

enddo

call run_timer ('STOP', time)

open(1, file = 'linearsolve_report.txt')
write(1, '(a, i8)' )          'Matrix size (N X N)       N = ', N
write(1, '(a, i8    )' )      'Number of random matrices   = ', ntests
write(1, '(a, 1es13.5, a )' ) 'total time for LA_GESVX = ', time, ' seconds'
write(1, '(a, 1es13.5)')      'Error = ', maxval(error)
close(1)



end program
