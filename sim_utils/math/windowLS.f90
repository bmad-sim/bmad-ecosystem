!+
! Module windowLS_mod
!
! This module uses the QR factorization to perform a least squares fit
! of a polynomial to a number of data points.  It is assumed that the data points
! are in (t1,y1), (t2,y2), ... ,(tn,yn) format
! where t1=t0, t2=t0+Dt, t3=t0+2Dt,...,tn=t0+(n-1)Dt.  Where the data points are
! separeted by fixed time intervals, the same QR factorization that can be used to fit points
! 1 through 10, can be used to fit points 2 through 11, and 3 through 12, and so on.
! This makes for an efficient method to fit polynomials to incoming data.
!
! This module takes as parameters:
!   1) The number of points to fit over.
!   2) The order of the polynomial.
!   3) The order of the derivative to return.
!
! This module is coded as an object and supports multiple instances.  An outline of how to use
! this module is:
! id = initFixedWindowLS(number_of_pts, delta_t_between_each_point, order_of_fit, order_of_derivative_returned)
! DO
!   <obtain new data point>
!   z = fixedWindowLS(new_y_data_point, id)
! ENDDO
! CALL destFixedWindowLS(id)
!
! initFixedWindowLS is the constructor and returns an id which must be used when calling fixedWindowLS.
! fixedWindowLS updates the fit with a new data point and returns the derivative of the polynomial at the new data point.
! destFixedWindowLS is the destructor.
!-
MODULE windowLS_mod

USE precision_def

IMPLICIT NONE

INTEGER, PARAMETER :: max_wls = 3 !maximum number of FixedWindowLS instances

! FixedWindowLS private data
TYPE wls_struct
  REAL(rp), ALLOCATABLE :: R1(:,:)
  REAL(rp), ALLOCATABLE :: Q1(:,:)
  INTEGER :: N     !number of data points to fit to
  INTEGER :: order !order of fit polynomial
  INTEGER :: der   !order of derivative to be returned by fixedWindowLS
  REAL(rp) :: xend
  REAL(rp), ALLOCATABLE :: y(:)  ! holds data
END TYPE wls_struct

TYPE(wls_struct), PRIVATE, SAVE :: wls(max_wls)
INTEGER, SAVE :: wls_ids = 0  !number of window LS instances

PUBLIC initFixedWindowLS
PUBLIC destFixedWindowLS
PUBLIC fixedWindowLS
PRIVATE shortFactorial

CONTAINS

!+
! Function initFixedWindowLS
!
! Initializes an instance of the fixed window least squares module.
! See module documentation (getf windowLS_mod) for use details.
! Any instance of windowLS created with this module should be destroyed with destFixedWindowLS.
!
! Input:
!   N     -- INTEGER, INTENT(IN): Number of data points to fit over. aka window size.
!   dt    -- REAL(rp), INTENT(IN): Time interval between data points. It is assumed that the data is 
!                                  separated by fixed time intervals.
!   order -- INTEGER, INTENT(IN): Order of fit polynomial.  Must be greater than or equal to der.
!   der   -- INTEGER, INTENT(IN): Order of derivative to be returned. Set der=0 to obtain the fit.
!
! Output:
!   <return value>  -- INTEGER: id of windowLS instance created.
!-
FUNCTION initFixedWindowLS(N,dt,order,der) RESULT(id)
  USE qr_mod

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: N      ! number of data points aka window size
  REAL(rp), INTENT(IN) :: dt    ! time between data points
  INTEGER, INTENT(IN) :: order  ! order of fit polynomial
  INTEGER, INTENT(IN) :: der    ! order of derivative to be returned by fixedWindowLS
                                ! set der = 0 to return the fit (zeroth derivative)
  INTEGER id
  INTEGER i,j ! loop counters
  REAL(rp), ALLOCATABLE :: A(:,:)
  REAL(rp) x

  wls_ids = wls_ids + 1
  IF(wls_ids .gt. max_wls) THEN
    WRITE(*,*) "Maximum number of fixed window LS instanced hard-coded to max_wls = ", max_wls
    STOP
  ENDIF
  id = wls_ids

  ! Copy constructor data to module private data
  wls(id)%N = N
  wls(id)%order = order
  wls(id)%der = der

  IF( der > order ) THEN
    WRITE(*,*) "initFixedWindowLS ERROR: derivative order cannot be greater than fit order"
    WRITE(*,*) "Halting ..."
    STOP
  ENDIF

  ALLOCATE(A(N,order+1))
  ALLOCATE(wls(id)%Q1(N,order+1))
  ALLOCATE(wls(id)%R1(order+1,order+1))
  ALLOCATE(wls(id)%y(N))
  wls(id)%y = 0.0_rp

  wls(id)%xend = N*dt

  A(1,:) = 0.0_rp
  A(1,1) = 1.0_rp
  DO i=2,N
    x = (i-1)*dt
    A(i,1) = 1.0
    DO j=2,order+1
      A(i,j) = x**(j-1)
    ENDDO
  ENDDO

  CALL thin_qr(A,wls(id)%Q1,wls(id)%R1)

  DEALLOCATE(A)

END FUNCTION initFixedWindowLS

!+
! Subroutine destFixedWindowLS
!
! Destroys an instance of fixedWindowLS.  Deallocates module memory.  This subroutine
! should be called before the program closes.
!
! Input:
!   id    -- INTEGER, INTENT(IN): id of instance to destroy
!
! Output:
!   none
!-
SUBROUTINE destFixedWindowLS(id)
  INTEGER, INTENT(IN) :: id

  DEALLOCATE(wls(id)%Q1)
  DEALLOCATE(wls(id)%R1)
  DEALLOCATE(wls(id)%y)
END SUBROUTINE destFixedWindowLS

!+
! Function fixedWindowLS
!
! Main function of the windowLS modult.  Each call to this function adds a data point to the fit
! and returns the derivative evaluated at the end of the window.  It is assumed that all data points
! are separeted by the same interval.  
! This module is initialized with zeros for all data points, and so the results are unreliable until
! a number of data points equal to N has been entered.
!
! initFixedWindowLS must be called prior to calling this function.  destFixedWindowLS should be
! called when the instance is no longer needed.
!
! Input:
!   ynew:     -- REAL(rp), INTENT(IN): New data point.
!
! Output:
!   <return value) -- REAL(rp): Derivative of fit polynomial evaluated at end of window.
!-
FUNCTION fixedWindowLS(ynew,id) RESULT(z)
  USE QR_mod

  IMPLICIT NONE

  REAL(rp), INTENT(IN) :: ynew
  INTEGER, INTENT(IN) :: id
  REAL(rp) P(wls(id)%order+1)
  REAL(rp) z  

  REAL(rp) b(wls(id)%N)  ! hold Q1'y
  INTEGER i

  !Drop oldest data point from bottom, add new data point to the top
  DO i=1,wls(id)%N-1
    wls(id)%y(i) = wls(id)%y(i+1) 
  ENDDO
  wls(id)%y(wls(id)%N) = ynew

  b = MATMUL(TRANSPOSE(wls(id)%Q1),wls(id)%y)
  CALL RbackS(wls(id)%R1,b,P)

  z = 0.0
  DO i=wls(id)%der, wls(id)%order
    z = z + shortFactorial(i,wls(id)%der)*P(i+1)*(wls(id)%xend**(i-wls(id)%der))
  ENDDO
END FUNCTION fixedWindowLS

!+
! Function shortFactorial
!
! Computes a truncated factorial.  Eg: shortFactorial(10,4) would return 10*9*8*7
!
! Input:
!   z    -- INTEGER, INTENT(IN): number to obtain truncated factorial of.
!   N    -- INTEGER, INTENT(IN): number of products to evaluate.
!
! Output:
!   <return value> -- INTEGER: Result of truncated factorial.
!-
FUNCTION shortFactorial(z,N) result(out)
  INTEGER, INTENT(IN) :: z
  INTEGER, INTENT(IN) :: N
  INTEGER out
  INTEGER i

  out = 1.0
  DO i=1,N
    out = out*(z-i+1)
  ENDDO
END FUNCTION shortFactorial

END MODULE windowLS_mod












