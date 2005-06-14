!###########################################################################
!+
! Subroutine gfit3D(phase_coords,parameters)
!
! Subroutine to fit a 3 Gaussians to a set of spatial coordinates
!
! Modules needed:
!   use precision_def
!   use bmad_struct
!
! Input:
!   phase_coords  -- Type(coord_struct),Dimension(0:ndata): Data points (x,x_angle,y,y_angle,z,z_angle)
!   parameters    -- Real(RP),Dimension(3,3): Parameters from previous cycle are 
!                       used as the guess values to fit parameters for this cycle
! Output: 
!   parameters    -- Real(RP),Dimension(3,3): Parameters of fitted Gaussian (A,mu,sigma by x,y,z)
!
! Calls      :
!
! Author     :
!
! Modified   :
!-
!........................................................................
!
! $Id$
!
! $Log$
! Revision 1.1  2005/06/14 14:59:02  cesrulib
! Initial revision
!
!
!........................................................................
!
#include "CESR_platform.h"

SUBROUTINE gfit3D(phase_coords,parameters)
  use precision_def
  use bmad_struct

  IMPLICIT NONE
  TYPE(coord_struct), DIMENSION(1:), INTENT(IN) :: phase_coords
  REAL(RP), DIMENSION(1:,1:), INTENT(INOUT):: parameters
  REAL(RP), DIMENSION(1:size(phase_coords)) :: xreal
  REAL(RP), DIMENSION(1:size(parameters,1)) :: params
  INTEGER :: i, ndata, coord
  INTEGER, PARAMETER :: MinDataVals = 40

  INTERFACE
     SUBROUTINE gfit1D(xreal,params)
       use precision_def
       IMPLICIT NONE
       REAL(RP), DIMENSION(1:), INTENT(IN)    :: xreal
       REAL(RP), DIMENSION(1:), INTENT(INOUT) :: params
     END SUBROUTINE gfit1D
  END INTERFACE

  ! Provide initial guess for amplitudes if not already specified
  ! Only necessary on first turn
  IF (parameters(1,1) == 0.0 .AND. parameters(1,2) ==0.0 .AND. parameters(1,3) == 0.0) THEN
     parameters(1,1) = 1.0E-002
     parameters(1,2) = 1.0E-004
     parameters(1,3) = 1.0E-006
  END IF

  ! if phase_coords is too small, we won't be able to fit the data
  IF(size(phase_coords,1) <= MinDataVals) PRINT *,"Too few data values to fit successfully. ", &
       "Using statistical means and variances." 
  
  open(12,file='gfit3d.chisq',status='old',access='append')

  ! fit data for x,y,z separately
  DO coord=1,5,2                                ! 1=x,2=y,3=z
     ! fill xreal with x, y, or z coordinates
     DO i=1,size(phase_coords)
        xreal(i) = phase_coords(i)%vec(coord)
     END DO
     params(:) = parameters(:,(coord+1)/2)      ! guess for params based on last calculation
     CALL gfit1D(xreal,params)
     parameters(:,(coord+1)/2) = params(:)      ! save calculated parameters
  END DO

  close(12)

! DON'T WRITE gfit3d.stat
!  open(36,file='gfit3d.stat',status='old',access='append')
!  write(36,'()')
!  close(36)


END SUBROUTINE gfit3D

!###########################################################################
! Subroutine gfit1D(xreal,params)
!
! Subroutine to fit a Gaussian to a set of data
!
! Modules needed:
!   use precision_def
!
! Input:
!   xreal  -- Real(RP),Dimension(ndata): Data points
!   params -- Real(RP),Dimension(3): Parameters from previous cycle are 
!                 used as the guess values to fit parameters for this cycle
!
! Output: 
!   params -- Real(RP),Dimension(3): Parameters of fitted Gaussian (A,mu,sigma)
!

SUBROUTINE gfit1D(xreal,params)
  use precision_def
  IMPLICIT NONE
  REAL(RP), DIMENSION(1:), INTENT(IN)    :: xreal   ! data points
  REAL(RP), DIMENSION(1:), INTENT(INOUT) :: params  ! parameters of gaussian func

  INTEGER  :: nbin                               ! number of bins
  INTEGER, PARAMETER :: MinDataVals = 40         ! minimum # of data values we can fit
  REAL(RP) :: xmin, xmax, rmean, rvar            ! statistical quantities of real data
  REAL(RP), DIMENSION(:), ALLOCATABLE :: xbin    ! position of bins
  REAL(RP), DIMENSION(:), ALLOCATABLE :: ybin    ! contents of bins
  INTEGER  :: AllocateStatus, i
  REAL(RP) :: width                              ! dummy variable

  INTERFACE
     SUBROUTINE statistics(xreal,xmin,xmax,rmean,rvar)
       use precision_def
       IMPLICIT NONE
       REAL(RP), DIMENSION(1:), INTENT(IN) :: xreal
       REAL(RP), INTENT(OUT) :: xmin, xmax, rmean, rvar
     END SUBROUTINE statistics
     SUBROUTINE makehisto(xreal,xmin,xmax,xbin,ybin)
       use precision_def
       IMPLICIT NONE
       REAL(RP), DIMENSION(1:), INTENT(IN) :: xreal
       REAL(RP), INTENT(IN) :: xmin, xmax
       REAL(RP), DIMENSION(1:), INTENT(OUT) :: xbin, ybin
     END SUBROUTINE makehisto
     SUBROUTINE fitter(xbin,ybin,params,rmean,rvar)
       use nr, ONLY : mrqmin
       use precision_def
       IMPLICIT NONE
       REAL(RP), DIMENSION(1:), INTENT(IN)    :: xbin, ybin
       REAL(RP), DIMENSION(1:), INTENT(INOUT) :: params
       REAL(RP), INTENT(IN) :: rmean,rvar
     END SUBROUTINE fitter
  END INTERFACE

  CALL statistics(xreal,xmin,xmax,rmean,rvar)

  ! we can not fit too few data values, so return now with statistical quantites and guess for A
  IF(size(xreal) <= MinDataVals) THEN
     DO i=1,size(params),3
        params(i+1) = rmean
        params(i+2) = rvar
     END DO
     RETURN
  END IF

  ! if no initial guess provided for gaussian mean & width, use statistical values for guesses
  DO i=1,size(params),3
     IF(params(i+1) == 0.) params(i+1) = rmean
     IF(params(i+2) == 0.) params(i+2) = rvar
  END DO

  ! calculate nbin using eqn from "Probability Density Estimation" by Don Johnson, 8-20-03
  width = ((9.0*(rvar**5.0)*SQRT(2*3.14159))/(2.0*(size(xreal))*exp(-(1.0/2.0)*(rmean/rvar)**2.0)*(1-(rmean/rvar)**2.0)**2.0))**(1.0/5.0)
  nbin = INT((xmax-xmin)/width)
  IF(nbin < 4) nbin = 4

  ALLOCATE(xbin(1:nbin), ybin(1:nbin), STAT = AllocateStatus)
     IF(AllocateStatus /= 0) STOP "*****Not Enough Memory*****"
 
  CALL makehisto(xreal,xmin,xmax,xbin,ybin)
  CALL fitter(xbin,ybin,params,rmean,rvar)

  DEALLOCATE(xbin, ybin)

END SUBROUTINE gfit1D

!###########################################################################
! Recursive Subroutine statistics(xreal,xmin,xmax,rmean,rvar)
!
! Subroutine to calculate statistical quantities of a set of data
!
! Modules needed:
!   use precision_def
!
! Input:
!   xreal -- Real(rp), Dimension(ndata): A set of data values to be evaluated
! Output: 
!   xreal -- Real(rp), Dimension(ndata): Values lying far away from the mean are marked 
!                                     by assigning them a value BadData=999.0
!   xmin  -- Real(rp): Smallest value in the data set
!   xmax  -- Real(rp): Largest value in the data set
!   rmean -- Real(rp): Mean of the data set
!   rvar  -- Real(rp): Variance of the data set
!
SUBROUTINE statistics(xreal,xmin,xmax,rmean,rvar)
  use precision_def
  
  IMPLICIT NONE
  REAL(RP), DIMENSION(1:), INTENT(INOUT) :: xreal
  REAL(RP), INTENT(OUT) :: xmin, xmax, rmean, rvar
  REAL(RP) :: xsum, x2sum
  REAL(RP), PARAMETER :: BadData = 999.0, CutoffSig = 4.0
  INTEGER  :: i,num_bad_data,redo_count
  INTEGER, SAVE :: num
  LOGICAL  :: redo, redo2


  redo_count = -1
  redo = .true.
  do while(redo)
     redo_count = redo_count + 1
     
     ! find xmin, xmax
     DO i=1,size(xreal)
        xmin = xreal(i)
        xmax = xreal(i)
        IF(xmin /= BadData .AND. xmax /= BadData) EXIT
     END DO
     DO i=1,size(xreal)
        IF (xreal(i) <  xmin .AND. xreal(i) /= BadData)  xmin = xreal(i)
        IF (xreal(i) >  xmax .AND. xreal(i) /= BadData)  xmax = xreal(i)
     END DO
     
     ! sum data points
     xsum  = 0
     x2sum = 0
     num_bad_data = 0
     DO i=1,size(xreal)
        IF(xreal(i) /= BadData) xsum = xsum + xreal(i)
        IF(xreal(i) /= BadData) x2sum = x2sum + xreal(i)**2
        IF(xreal(i) == BadData) num_bad_data = num_bad_data + 1
     END DO
     
     ! do statistical calculations on sums
     rmean = xsum/(size(xreal)-num_bad_data)
     rvar = SQRT((x2sum/(size(xreal)-num_bad_data)) - rmean**2)
     
     ! eliminate outliers from xreal & recalculate mean,variance  
     redo = .false.
     DO i=1,size(xreal)
        IF ( (xreal(i) /= BadData)         .AND. &
             ( xreal(i) >= (rmean+CutoffSig*rvar)   .OR.  &
             xreal(i) <= (rmean-CutoffSig*rvar) ) ) THEN
           xreal(i) = BadData
           redo = .true.
           redo2= .true.
        END IF
     END DO
  end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  redo = .false.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!     IF (redo) THEN
!        CAll statistics(xreal,xmin,xmax,rmean,rvar)
!        RETURN
!     ELSE
 
  num_bad_data = 0
  DO i=1,size(xreal)
     IF(xreal(i) == BadData) num_bad_data = num_bad_data + 1
  END DO


! DON'T WRITE gfit3d.stat
     ! write number of data points used to gfit3d.stat file
!     open(unit=36, file='gfit3d.stat', access="append")
!     if(redo2) then
!        write(36,'(3i6,f6.2,A10)') size(xreal), num_bad_data, size(xreal)-num_bad_data,CutoffSig, '.true.'
!     else
!        write(36,'(3i6,f6.2,A10)') size(xreal), num_bad_data, size(xreal)-num_bad_data,CutoffSig, '.false.'
!     end if
!     close(36)
!  RETURN
!END IF

END SUBROUTINE statistics

!###########################################################################
! Subroutine makehisto(xreal,xmin,xmax,xbin,ybin)
!
! Subroutine to sort a set of data values into bins
!
! Modules needed:
!   use precision_def
!
! Input:
!   xreal -- Real(rp), Dimension(ndata): Set of data to be sorted
!   xmin  -- Real(rp): Smallest value in the data set 
!   xmax  -- Real(rp): Largest value in the data set
! Output:
!   xbin  -- Real(rp), Dimension(nbin): Positions of bins
!   ybin  -- Real(rp), Dimension(nbin): Number of data points falling into a bin
!
SUBROUTINE makehisto(xreal,xmin,xmax,xbin,ybin)
  use precision_def

  IMPLICIT NONE
  REAL(RP), DIMENSION(1:), INTENT(IN) :: xreal
  REAL(RP), INTENT(IN) :: xmin, xmax
  REAL(RP), DIMENSION(1:), INTENT(OUT) :: xbin, ybin
  REAL(RP) :: binwidth
  REAL(RP), PARAMETER :: BadData = 999.0
  INTEGER  :: i, j, k

  binwidth = (xmax-xmin)/size(xbin)               ! width of a bin
  DO k = 1,size(xbin)
     xbin(k) = xmin + (k-(1/2.0))*binwidth
     ybin(k) = 0.0
  END DO

  ! sort each real data point into a bin
  DO i=1,size(xreal)                                ! scroll through data points
     DO j=1,size(xbin)                              ! scroll through bins
        IF ((xreal(i) >=  (xmin + (j-1)*binwidth)) .AND. &
            (xreal(i) <=  (xmin +  (j) *binwidth)) .AND. &
            (xreal(i) /= BadData)) THEN            ! outliers were given value BadData=999.0 
           ybin(j) = ybin(j) + 1.0
           EXIT
           ! account for min & max points lost by roundoff
        ELSE IF (xreal(i) >= xmax .AND. xreal(i) /= BadData) THEN
           ybin(size(xbin)) = ybin(size(xbin)) + 1.0
           EXIT
        ELSE IF (xreal(i) <= xmin .AND. xreal(i) /= BadData) THEN
           ybin(1) = ybin(1) + 1.0
           EXIT
        END IF
     END DO
  END DO

! DON'T WRITE histo.dat
!  open(unit=3,file='histo.dat',access="append")
!  do i=1,size(xbin)
!     write(3,"(e12.4,f6.1)") xbin(i),ybin(i)
!  end do
!  write(3,"(A)") ''
!  close(3)
  
END SUBROUTINE makehisto

!###########################################################################
! Subroutine fitter(xbin,ybin,params,rmean,rvar)
!
! Subroutine to fit a Gaussian curve to data coordinates
!
! Modules needed:
!   use precision_def
!
! Input:
!   xbin   -- Real(rp), Dimension(nbin): Position of bins (x-coordinate)
!   ybin   -- Real(rp), Dimension(nbin): Contents of bins (y-coordinate)
!   rmean  -- Real(rp): Mean of data set; Used as initial guess
!   rvar   -- Real(rp): Variance of data set; used as initial guess
!   params -- Real(rp), Dimension(:): Parameters of previous fit are used as an initial guess
! Output: 
!   params -- Real(rp), Dimension(:): Parameters of fitted Gaussian (A,mu,sigma)
!


SUBROUTINE fitter(xbin,ybin,params,rmean,rvar)
  use nr, ONLY : mrqmin
  use precision_def
  
  IMPLICIT NONE
  REAL(RP), DIMENSION(1:), INTENT(IN) :: xbin, ybin
  REAL(RP), INTENT(IN) :: rmean,rvar
  REAL(RP), DIMENSION(1:), INTENT(INOUT) :: params

  REAL(RP), PARAMETER :: Eps = 0.01, MaxParam=1.0E+006
  REAL(RP), PARAMETER :: PI = 3.1415926535897932
  INTEGER,  PARAMETER :: MaxCount = 20
  INTEGER :: MaxCounter
  REAL(RP), DIMENSION(1:size(params)) :: a,ainit  ! parameters that mrqmin solves for
  REAL(RP), DIMENSION(1:size(xbin)) :: sig        ! variance of each data point
  REAL(RP), DIMENSION(1:size(params),1:size(params)) :: covar, alpha
  REAL(RP) :: chisq, chisqprev                    ! parameters ready when Chi^2 converges
  REAL(RP) :: alamda                              ! controls mrqmin
  REAL(RP) :: factor
  LOGICAL,  DIMENSION(1:size(params)) :: maska    ! true=fit param, false=frozen
  LOGICAL  :: blowup                              ! tells if parameters are blowing up
  LOGICAL  :: stall                               ! indicates when mrqmin has stalled
  INTEGER  :: i,j,count,counter

  INTERFACE
     SUBROUTINE funcs(x,a,y,dyda)
       use precision_def
       REAL(RP), DIMENSION(:), INTENT(IN):: x, a
       REAL(RP), DIMENSION(:), INTENT(OUT) :: y
       REAL(RP), DIMENSION(:,:), INTENT(OUT) :: dyda
     END SUBROUTINE funcs
  END INTERFACE

  ! initialization
  alamda = -1.0                         ! for initialization of mrqmin
  chisqprev = 1.0                       ! guess for chisq prev
  maska(:) = .true.                     ! fit all params
  DO i=1,size(params),3                 ! inital guess for parameters
     ainit(i)   = params(i)/(rvar*SQRT(2.0*PI))
     ainit(i+1) = rmean
     ainit(i+2) = rvar*SQRT(2.0)
  END DO
  a(:) = ainit(:)
  factor  = 10.0
  sig(:)  = 1                   ! guess for sig
  count   = 0
  counter = 0
  MaxCounter = 500*size(xbin)

  ! call mrqmin until parameters converge to the correct values (Chi^2 stops changing)
  DO

! DON'T WRITE gfit3d.iter     
!     open(unit=12,file='gfit3d.iter',access='append')
!     write(12,'(i8,2e20.8,3e20.8,3e20.8)') counter,chisq,chisqprev,a(1:3),ainit(1:3)
!     close(12)

     ! if convergence is impossible (usually due to too few data points), stop iterating
     counter = counter + 1
     IF(counter >= MaxCounter) THEN
        PRINT *,"Fitter could not converge.  Using statistical mean and variance."
        EXIT
     END IF
        
     CALL mrqmin(xbin,ybin,sig,a,maska,covar,alpha,chisq,funcs,alamda)

     ! count when chisq has not changed over the past iteration
     IF (chisq == chisqprev) THEN
        count = count + 1
     ELSE 
        count = 0
     END IF

     ! check for errors
     blowup = .false.
     DO i=1,size(a)
        IF (ABS(a(i)) > MaxParam) blowup = .true.
     END DO
     stall = .false.
     IF (count >= MaxCount) stall = .true.

     ! let params converge, then check if they've blown up; if not, exit
     IF ((chisq < chisqprev) .AND. (ABS(chisq - chisqprev) < Eps)) THEN
        IF (.not. blowup) EXIT
     END IF
 
     ! if params have blown up, clear mrqmin & reset it using new initial conditions
     IF (blowup .or. stall) THEN
        alamda = 0.0                         ! clear mrqmin (deallocate saved arrays)
           CALL mrqmin(xbin,ybin,sig,a,maska,covar,alpha,chisq,funcs,alamda)
        alamda = -1.0                        ! start over w/ new inital guess
           ! make inital conditions wander around until we find some that work
           DO i=1,size(ainit),3              ! 1,4,7... all the A's
              IF (ABS(ainit(i)) < (params(i)/(rvar*SQRT(2.0*PI)))*1.0E-003)  factor=5.0
              IF (ABS(ainit(i)) > (params(i)/(rvar*SQRT(2.0*PI)))*1.0E+003)  factor=0.2
              ainit(i) = ainit(i) * factor
           END DO
           a(:) = ainit(:)
           count = 0
           chisqprev = 1.0
     ! if params have not blown up, iterate as usual
     ELSE
        chisqprev = chisq
     END IF
  END DO

! WRITE gfit3d.chisq
  write(12,'(8x,e20.8)',ADVANCE='NO') chisq
  
  ! call mrqmin w/ alamda=0 to deallocate saved arrays
  alamda = 0
  CALL mrqmin(xbin,ybin,sig,a,maska,covar,alpha,chisq,funcs,alamda)

  ! if we were unable to fit a curve, use statistical values
  IF(counter >= MaxCounter) THEN
     DO i=1,size(params),3
        params(i+1) = rmean
        params(i+2) = rvar
     END DO
  ! otherwise store the fitted parameters (A,MU,SIGMA) to be passed out
  ELSE
     ! make sure A,sigma > 0
     DO i=1,size(a),3
        IF (a(i+2) < 0) a(i+2) = -a(i+2)
     END DO
     DO i=1,size(params),3
        params(i) = a(i)*a(i+2)*SQRT(PI)
        params(i+1) = a(i+1)
        params(i+2) = a(i+2)/SQRT(2.0_RP)
     END DO
  END IF

END SUBROUTINE fitter

! same as fgauss from NR
SUBROUTINE funcs(x,a,y,dyda)
  USE nrtype; USE nrutil, ONLY : assert_eq
  USE precision_def
  REAL(RP), DIMENSION(:),   INTENT(IN)  :: x, a
  REAL(RP), DIMENSION(:),   INTENT(OUT) :: y
  REAL(RP), DIMENSION(:,:), INTENT(OUT) :: dyda
  INTEGER(I4B) :: i, na, nx
  REAL(RP), DIMENSION(size(x)) :: arg, ex, fac
  nx=assert_eq(size(x),size(y),size(dyda,1),'fgauss: nx')
  na=assert_eq(size(a),size(dyda,2),'fgauss: na')
  y(:)=0.0
  DO i=1,na-1,3
     arg(:)=(x(:)-a(i+1))/a(i+2)
     ex(:)=exp(-arg(:)**2)
     fac(:)=a(i)*ex(:)*2.0_rp*arg(:)
     y(:)=y(:)+a(i)*ex(:)
     dyda(:,i)=ex(:)
     dyda(:,i+1)=fac(:)/a(i+2)
     dyda(:,i+2)=fac(:)*arg(:)/a(i+2)
  END DO
END SUBROUTINE funcs
