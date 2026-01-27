!+
!  Module qr_mod
!
!  This module contains routines for performing the QR factorization on m*n, m>=n, matrices.
!  One of the uses of QR is in least squares fitting of a polynomial to data.
!  For a detailed description of how to use QR for least squares fitting, see section 5.3.3
!  of Matrix Computations by Golub and Van Loan.  In short, QR decomposes an n*m matrix
!  into a square orthogonal matrix (Q) times a upper triangular matrix (R).  R is a 
!  rectangular matrix in the sense that only the upper triangular portion of its top
!  m*m elements are nonzero.
!  With QR in hand, the least squares solution to Ax-b is obtained by solving Rx-Qb for x.
!
!  This module contains the following subroutines and functions
!   SUBROUTINE qr(A,Q,R) - Returns Q and R for a given A
!   SUBROUTINE thin_qr(A,Q1,R1) - Returns only those parts of Q and R necessary for LS fitting
!   FUNCTION twonorm2(v) - Returns the 2-norm of a vector.
!   FUNCTION house(v) - Returns the householder vector of a vector.  (Used by QR decomposition)
!   SUBROUTINE RbackS(R,b,x) - Solved Rx=b for x assuming R is upper triangular.
!   SUBROUTINE print_mat(A) - Prints A in row,column format.
!-
MODULE qr_mod

USE sim_utils_interface

IMPLICIT NONE

PUBLIC qr
PUBLIC thin_qr
PUBLIC RbackS
PUBLIC print_mat
PRIVATE twonorm2
PRIVATE house

CONTAINS

!+
! Subroutne qr(Ain,Q,R)
!
! Produces QR factorization of an m*n, m>=n, matrix using Algorithm 5.2.1 from Golub and Van Loan.
! This subroutine is flop efficient, but has extra overhead since it builds the full Q matrix, rather
! than returning just the householder vectors.
!
! Input:
!   Ain(:,:)   -- REAL(rp), INTENT(IN):  m*n, m>=n, array of reals.
! 
! Output:
!   Q(:,:)     -- Double precision m*m array of reals.
!   R(:,:)     -- Double precision n*n array of reals.
!-
SUBROUTINE qr(Ain,Q,R)
  REAL(rp), INTENT(IN) :: Ain(:,:)
  REAL(rp) A(SIZE(Ain,1),SIZE(Ain,2))
  REAL(rp), INTENT(OUT) :: Q(:,:)
  REAL(rp), INTENT(OUT) :: R(:,:)
  REAL(rp) V(SIZE(Ain,1),SIZE(Ain,2))
  INTEGER j,m,n
  INTEGER k

  A = Ain
  m = SIZE(A,1)  !number of rows
  n = SIZE(A,2)  !number of columns

  V = 0.0
  DO j=1,n
    V(j:m,j) = house(A(j:m,j))
    A(j:m,j:n) = A(j:m,j:n) - outer_product(V(j:m,j),MATMUL(V(j:m,j),A(j:m,j:n)))
  ENDDO

  R=0.0
  DO j=1,n
    R(j,j:n) = A(j,j:n)
  ENDDO

  Q=0.0
  DO j=1,m
    Q(j,j) = 1.0
  ENDDO
  DO j=n,1,-1
    Q(j:m,j:m) = Q(j:m,j:m) - outer_product(V(j:m,j),MATMUL(V(j:m,j),Q(j:m,j:m)))
  ENDDO
END SUBROUTINE qr

!+
! Subroutine thin_qr(Ain,Q1,R1)
! 
! Computes thin QR factorization of a m*n, m >= n,  matrix.
! Implements modified Gram-Schmidt algorithm 5.2.5 from Golub and Van Loan, which
! produces only the first n columns of Q.
! Only the first n columns of Q are necessary for LS fitting via QR.  This algorithm
! is more flop and memory efficient than the full QR factorization.
!
! Input:
!   Ain(:,:)   -- REAL(rp), INTENT(IN): m*n, m>=n, array of reals.
! 
! Output:
!   Q1(:,:)    -- REAL(rp), INTENT(OUT): m*n array of reals.
!   R1(:,:)    -- REAL(rp), INTENT(OUT): n*n array of reals.
!-
SUBROUTINE thin_qr(Ain,Q1,R1)
  REAL(rp), INTENT(IN) :: Ain(:,:)
  REAL(rp) A(SIZE(Ain,1),SIZE(Ain,2))
  REAL(rp), INTENT(OUT) :: Q1(:,:)
  REAL(rp), INTENT(OUT) :: R1(:,:)
  INTEGER m,n
  INTEGER k,j

  A = Ain
  m = SIZE(A,1)  !number of rows
  n = SIZE(A,2)  !number of columns

  Q1 = 0.0
  R1 = 0.0
  DO k = 1,n
    R1(k,k) = SQRT(twonorm2(A(1:m,k)))
    Q1(1:m,k) = A(1:m,k)/R1(k,k)
    DO j=k+1,n
      R1(k,j) = DOT_PRODUCT(Q1(1:m,k),A(1:m,j))
      A(1:m,j) = A(1:m,j)-Q1(1:m,k)*R1(k,j)
    ENDDO
  ENDDO
END SUBROUTINE thin_qr

!+
! Function twonorm2(v)
!
! Returns the 2-norm of a vector v.  The 2-norm of a vector is the sum of the squares of
! its elements.
! 
! Input:
!   v(:)   -- REAL(rp), INTENT(IN): 1-dim array of reals.
!
! Output:
!   z      -- REAL(rp): 2-norm of v.
!-
FUNCTION twonorm2(v) RESULT(z)
  REAL(rp), INTENT(IN) :: v(:) 
  INTEGER i,m
  REAL(rp) z

  m = SIZE(v)

  z = 0.0
  DO i = 1,m
    z = z + v(i)*v(i)
  ENDDO
END FUNCTION twonorm2

!+
! Function house(x)
!
! Returns the householder vector v for the input vector x.
! Used in calculation of full QR factorization.
! Based on algorithm 5.1.1 from Golub and Van Loan.
!
! Input:
!   x(:)  -- REAL(rp), INTENT(IN): length n vector of reals.
!
! Output:
!   v(:)  -- REAL(rp), INTENT(IN): length n vector of reals.
!-
FUNCTION house(x) RESULT(v)
  REAL(rp), INTENT(IN) :: x(:)
  REAL(rp) s
  REAL(rp) v(SIZE(x))
  REAL(rp) beta, mu
  INTEGER n

  n = SIZE(x)
  s = dot_product(x(2:n),x(2:n))
  v(1) = 1.0
  v(2:n) = x(2:n)

  IF(s .eq. 0.0) THEN
    beta = 0
  ELSE
    mu = sqrt( x(1)*x(1) + s)
    IF( x(1) .le. 0 ) THEN
      v(1) = x(1) - mu
    ELSE
      v(1) = -s/(x(1)+mu)
    ENDIF
  ENDIF
  beta = 2.0*v(1)*v(1)/(s + v(1)*v(1))
  v = v/v(1)
  v = sqrt(beta)*v
END FUNCTION house

!+
! Subroutine RbackS(R,b,x)
!
! Solves Rx=b for x assuming R is upper triangular.
!
! Input:
!   R(:,:)    -- REAL(rp), INTENT(IN): n*n matrix.  Assumed upper triangular.
!   b(:)      -- REAL(rp), INTENT(IN): length n vector.
!
! Output:
!   x(:)      -- REAL(rp), INTENT(IN): length n vector. Solution to Rx=b.
!-
SUBROUTINE RbackS(R,b,x)
  IMPLICIT NONE

  REAL(rp), INTENT(IN) :: R(:,:)
  REAL(rp), INTENT(IN) :: b(:)
  REAL(rp), INTENT(OUT) :: x(:)
  REAL(rp) term1
  INTEGER m
  INTEGER i,j

  m = SIZE(R,1)

  x(m) = b(m)/R(m,m)
  DO i=m-1,1,-1
    term1 = 0.0
    DO j=i+1,m
      term1 = term1 + R(i,j)*x(j) 
    ENDDO 
    x(i) = ( b(i)-term1 ) / R(i,i)
  ENDDO
END SUBROUTINE RbackS

!+
! Subroutine print_mat(A)
!
! Formatted display of a n,m matrix in A(row,column) format.
!
! Input:
!   A(:,:)    -- REAL(rp), INTENT(IN): n*m matrix.
!
! Output:
!   none -- (writes to stdout)
!-
SUBROUTINE print_mat(A)
  REAL(rp), INTENT(IN) :: A(:,:)
  INTEGER m,n
  INTEGER i
  CHARACTER*7 format
  CHARACTER*1 n_str

  m = SIZE(A,1)
  n = SIZE(A,2)

  WRITE(n_str,'(I1)') n
  format = '('//n_str//'F8.4)'

  DO i=1,m
    WRITE(*,format) A(i,:)
  ENDDO
END SUBROUTINE print_mat

END MODULE qr_mod
