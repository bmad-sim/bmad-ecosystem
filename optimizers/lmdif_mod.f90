module lmdif_mod

use precision_def

contains

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!+
! function opti_lmdif (vec, n, merit, eps) result (this_opti)
!
! Function which tries to get the merit function/s as close to zero as possible
! by changing the values in vec. Multiple merit functions can be used.
!
! Initially call Initial_lmdif to initialize everything.
!
! Input:
!   vec(:)  -- real(dp): Array of variables to vary
!   n       -- integer: Max number of iterations
!   merit   -- Function to minimize. The interface is:
!                   function merit (vec, flag) result (mrslt)
!                     use precision_def
!                     real(dp) vec(:)
!                     real(dp) mrslt
!                     integer flag
!                   end function
!   eps     -- real(dp): Desired deviation from 0 for your merit value/s
!
! Output:
!   this_opti  -- real(dp): Result/s of merit function at optimized values.
!-

function opti_lmdif (vec, n, merit, eps) result (this_opti)

  implicit none

  interface
    function merit (vec, flag) result (mrslt)
      use precision_def
      real(dp) vec(:)
      real(dp) mrslt
      integer flag
    end function
  end interface

  real(dp) vec(:), this_opti
  integer flag
  real(dp) eps, fct(1)
  integer n, nfv
  logical at_end

!

  flag = 0.D+0
  nfv  = 1
  
  do 
    fct(1) = merit(vec, flag)
    call suggest_lmdif (vec, fct, eps, n, at_end)
    if (at_end) exit
  enddo

  this_opti = merit(vec, flag)

end function opti_lmdif

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!+
! subroutine initial_lmdif
!
! Subroutine that clears out previous saved values of the optimizer.
!
! Input: none
!
! Output: none
!-

subroutine initial_lmdif
  real(dp) xv(1), fv(1), eps
  integer itermx
  logical at_end

  XV = 0.D+0
  FV = 0.D+0
  EPS = 0.D+0
  ITERMX = 0
  call suggest_lmdif (XV, FV, EPS, ITERMX, at_end, .TRUE.)

end subroutine initial_lmdif

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!+
! subroutine suggest_lmdif (xv, fv, eps, itermx, at_end, reset_flag)
!
! Reverse communication subroutine. 
! It suggests values for your input variables based on
! the previous value of your merit function.
!
! Use initial_lmdif to initialize internal variables
!
! Input:
!   xv(:)      -- real(dp): Array of variables
!   fv(:)      -- real(dp): Array of function value/s that should be optimized to zero
!   eps        -- real(dp): Desired accuracy with which the optimum should be found.
!   itermx     -- integer: Max number of iterations
!   reset_flag -- logical: Optional. Used by initial_lmdif to clear 
!                           previous saved values
!
! Output:
!   xv(:)      -- real(dp): Suggested new values
!   fv(:)      -- real(dp): After the last optimization this returns the best values ever.
!   at_end     -- logical: Set to False if more optimization is recommended.
!                   If set to True then xv(:) will be the minimum found.
!-

subroutine suggest_lmdif (XV, FV, EPS, ITERMX, at_end, reset_flag)

      IMPLICIT none

      integer itermx
      integer NV, NFV, i, nf
      integer, save :: ipos, ipos1, info, nfev, iter
      real(dp) eps
      real(dp), save :: s, smin
      real(dp) XV(:)
      integer, allocatable, save :: IPVT(:)
      real(dp), allocatable, save :: XMIN(:), DIAG(:), QTF(:), FJAC(:,:)
      real(dp), allocatable, save :: WA4(:), WA1(:),WA2(:),WA3(:)
      real(dp), allocatable, save :: fvv(:)

      integer, save :: old_nv = -1

      real(dp) FV(:)
      real(dp), allocatable, save :: FMIN(:)

      integer, save :: old_nfv = -1

      real(dp) ftol, xtol, gtol, epsfcn, factor
      integer mode, nprint, nend, ldfjac

      logical, optional :: reset_flag
      logical at_end, error

!

      nv = size (xv)
      nfv = size(fv)
      nf = max(nv,nfv)

      if(present(reset_flag).or.(old_nv.ne.nv).or.(old_nfv.ne.nfv)) then
         iter = 0
         ipos=0
         ipos1=0
         info=0
         nfev=0
         s=0.D+0
         smin=1.D+99
         old_nv=-1
         old_nfv=-1
         ldfjac = nfv
         CALL LMDIF(NFV,NV,XV,FV,FTOL,XTOL,GTOL,EPS,ITERMX,EPSFCN, &
              DIAG,MODE,FACTOR,NPRINT,INFO,NFEV,FJAC,LDFJAC, &
              IPVT,QTF,WA1,WA2,WA3,WA3,NEND,IPOS,IPOS1, error, .TRUE.)
         if (error) then
           at_end = .true.
           return
         endif
         if (present(reset_flag)) return
      endif

      if((old_nv.ne.nv).or.(old_nfv.ne.nfv)) then
         old_nv=nv
         old_nfv=nfv
         if(allocated(ipvt)) deallocate(ipvt)
         if(allocated(xmin)) deallocate(xmin)
         if(allocated(diag)) deallocate(diag)
         if(allocated(qtf)) deallocate(qtf)
         if(allocated(wa1)) deallocate(wa1)
         if(allocated(wa2)) deallocate(wa2)
         if(allocated(wa3)) deallocate(wa3)
         if(allocated(wa4)) deallocate(wa4)
         if(allocated(fjac)) deallocate(fjac)
         if(allocated(fmin)) deallocate(fmin)
         if(allocated(fvv)) deallocate(fvv)
         allocate(ipvt(nv),xmin(nv),diag(nv),qtf(nv),wa1(nv),wa2(nv),wa3(nv))
         allocate(fjac(nf,nv))
         allocate(wa4(nf),fmin(nf),fvv(nf))
         ipvt=0
         xmin=0
         diag=0
         qtf=0
         wa1=0
         wa2=0
         wa3=0
         wa4=0
         fjac=0
         fmin=0
      endif


!
!     SPECIAL CASES
!     *************
!

      S = 0.D0
      DO I=1, NFV
        S = S + FV(I)**2
      enddo

      IF(ITERMX.EQ.0) THEN
         at_end = .true.
         IPOS = 0
         IPOS1 = 0
         RETURN
      ELSEIF(S.LE.SMIN) THEN
         SMIN = S
         FMIN(1:NFV) = FV(1:nfv)
         XMIN(1:nv) = XV(1:nv)
      ENDIF

!     LMDIF
!     *****
!
      FTOL = 1.0D-15
      XTOL = 1.0D-15
      GTOL = 1.0D-15
      EPSFCN = 1.0D-15
      MODE = 1
      FACTOR = 100.0D0
      NPRINT = 0
      fvv = 0
      fvv(1:nfv) = fv

      LDFJAC = nf
      CALL LMDIF(nf,NV,XV,fvv,FTOL,XTOL,GTOL,EPS,ITERMX,EPSFCN, &
                DIAG,MODE,FACTOR,NPRINT,INFO,NFEV,FJAC,LDFJAC, &
                IPVT,QTF,WA1,WA2,WA3,WA4,NEND,IPOS,IPOS1, error)  
      if (error) then
        at_end = .true.
        return
      endif

      iter = iter + 1
!
      GO TO 200
!v
!     
!     AFTER OPTIMIZERS
!     ****************
!
 200  CONTINUE
!
      IF((iter.GE.itermx) .OR.(NEND.EQ.1)) THEN
         DO I=1,NV
            XV(I) = XMIN(I)
         ENDDO
         DO I=1,NFV
            FV(I) = FMIN(I)
         ENDDO
         at_end = .true.
         ITER = 0
         smin = 1.D+99
         IPOS = 0
         IPOS1 = 0
         RETURN
      ELSE
         at_end = .false.
         RETURN
      ENDIF
      END SUBROUTINE
!

!*****************
! LMDIF PACKAGE ROUTINES
!********************
!     ***************************************************************
      SUBROUTINE LMDIF(M,N,X,FVEC,FTOL,XTOL,GTOL,FATOL,MAXFEV,EPSFCN, &
                      DIAG,MODE,FACTOR,NPRINT,INFO,NFEV,FJAC,LDFJAC, &
                      IPVT,QTF,WA1,WA2,WA3,WA4,IEND,IPOS,IPOS1, error, reset_flag)

      INTEGER M,N,MAXFEV,MODE,NPRINT,INFO,NFEV,LDFJAC,IEND,IPOS,IPOS1
      INTEGER IPVT(N)
      real(dp) FTOL,XTOL,GTOL,FATOL,EPSFCN,FACTOR
      real(dp) X(N),FVEC(M),DIAG(N),FJAC(LDFJAC,N),QTF(N), &
                      WA1(N),WA2(N),WA3(N),WA4(M)
!      EXTERNAL FCN
!     ***************************************************************
!
!     SUBROUTINE LMDIF
!
!     THE PURPOSE OF LMDIF IS TO MINIMIZE THE SUM OF THE SQUARES OF
!     M NONLINEAR FUNCTIONS IN N VARIABLES BY A MODIFICATION OF
!     THE LEVENBERG-MARQUARDT ALGORITHM. THE USER MUST PROVIDE A
!     SUBROUTINE WHICH CALCULATES THE FUNCTIONS. THE JACOBIAN IS
!     THEN CALCULATED BY A FORWARD-DIFFERENCE APPROXIMATION.
!
!     THE SUBROUTINE STATEMENT IS
!
!       SUBROUTINE LMDIF(M,N,X,FVEC,FTOL,XTOL,GTOL,FATOL,MAXFEV,EPSFCN,
!                        DIAG,MODE,FACTOR,NPRINT,INFO,NFEV,FJAC,
!                        LDFJAC,IPVT,QTF,WA1,WA2,WA3,WA4)
!
!     WHERE
!
!       M IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!         OF FUNCTIONS.
!
!       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!         OF VARIABLES. N MUST NOT EXCEED M.
!
!       X IS AN ARRAY OF LENGTH N. ON INPUT X MUST CONTAIN
!         AN INITIAL ESTIMATE OF THE SOLUTION VECTOR. ON OUTPUT X
!         CONTAINS THE FINAL ESTIMATE OF THE SOLUTION VECTOR.
!
!       FVEC IS AN OUTPUT ARRAY OF LENGTH M WHICH CONTAINS
!         THE FUNCTIONS EVALUATED AT THE OUTPUT X.
!
!       FTOL IS A NONNEGATIVE INPUT VARIABLE. TERMINATION
!         OCCURS WHEN BOTH THE ACTUAL AND PREDICTED RELATIVE
!         REDUCTIONS IN THE SUM OF SQUARES ARE AT MOST FTOL.
!         THEREFORE, FTOL MEASURES THE RELATIVE ERROR DESIRED
!         IN THE SUM OF SQUARES.
!
!       XTOL IS A NONNEGATIVE INPUT VARIABLE. TERMINATION
!         OCCURS WHEN THE RELATIVE ERROR BETWEEN TWO CONSECUTIVE
!         ITERATES IS AT MOST XTOL. THEREFORE, XTOL MEASURES THE
!         RELATIVE ERROR DESIRED IN THE APPROXIMATE SOLUTION.
!
!       GTOL IS A NONNEGATIVE INPUT VARIABLE. TERMINATION
!         OCCURS WHEN THE COSINE OF THE ANGLE BETWEEN FVEC AND
!         ANY COLUMN OF THE JACOBIAN IS AT MOST GTOL IN ABSOLUTE
!         VALUE. THEREFORE, GTOL MEASURES THE ORTHOGONALITY
!         DESIRED BETWEEN THE FUNCTION VECTOR AND THE COLUMNS
!         OF THE JACOBIAN.
!
!       MAXFEV IS A POSITIVE INTEGER INPUT VARIABLE. TERMINATION
!         OCCURS WHEN THE NUMBER OF CALLS TO FCN IS AT LEAST
!         MAXFEV BY THE END OF AN ITERATION.
!
!       EPSFCN IS AN INPUT VARIABLE USED IN DETERMINING A SUITABLE
!         STEP LENGTH FOR THE FORWARD-DIFFERENCE APPROXIMATION. THIS
!         APPROXIMATION ASSUMES THAT THE RELATIVE ERRORS IN THE
!         FUNCTIONS ARE OF THE ORDER OF EPSFCN. IF EPSFCN IS LESS
!         THAN THE MACHINE PRECISION, IT IS ASSUMED THAT THE RELATIVE
!         ERRORS IN THE FUNCTIONS ARE OF THE ORDER OF THE MACHINE
!         PRECISION.
!
!       DIAG IS AN ARRAY OF LENGTH N. IF MODE = 1 (SEE
!         BELOW), DIAG IS INTERNALLY SET. IF MODE = 2, DIAG
!         MUST CONTAIN POSITIVE ENTRIES THAT SERVE AS
!         MULTIPLICATIVE SCALE FACTORS FOR THE VARIABLES.
!
!       MODE IS AN INTEGER INPUT VARIABLE. IF MODE = 1, THE
!         VARIABLES WILL BE SCALED INTERNALLY. IF MODE = 2,
!         THE SCALING IS SPECIFIED BY THE INPUT DIAG. OTHER
!        VALUES OF MODE ARE EQUIVALENT TO MODE = 1.
!
!       FACTOR IS A POSITIVE INPUT VARIABLE USED IN DETERMINING THE
!         INITIAL STEP BOUND. THIS BOUND IS SET TO THE PRODUCT OF
!         FACTOR AND THE EUCLIDEAN NORM OF DIAG*X IF NONZERO, OR ELSE
!         TO FACTOR ITSELF. IN MOST CASES FACTOR SHOULD LIE IN THE
!         INTERVAL (.1,100.). 100. IS A GENERALLY RECOMMENDED VALUE.
!
!       NPRINT IS AN INTEGER INPUT VARIABLE THAT ENABLES CONTROLLED
!         PRINTING OF ITERATES IF IT IS POSITIVE. IN THIS CASE,
!         FCN IS CALLED WITH IFLAG = 0 AT THE BEGINNING OF THE FIRST
!         ITERATION AND EVERY NPRINT ITERATIONS THEREAFTER AND
!         IMMEDIATELY PRIOR TO RETURN, WITH X AND FVEC AVAILABLE
!         FOR PRINTING. IF NPRINT IS NOT POSITIVE, NO SPECIAL CALLS
!         OF FCN WITH IFLAG = 0 ARE MADE.
!
!       INFO IS AN INTEGER OUTPUT VARIABLE. IF THE USER HAS
!         TERMINATED EXECUTION, INFO IS SET TO THE (NEGATIVE)
!         VALUE OF IFLAG. SEE DESCRIPTION OF FCN. OTHERWISE,
!         INFO IS SET AS FOLLOWS.
!
!         INFO = 0  IMPROPER INPUT PARAMETERS.
!
!         INFO = 1  BOTH ACTUAL AND PREDICTED RELATIVE REDUCTIONS
!                   IN THE SUM OF SQUARES ARE AT MOST FTOL.
!
!         INFO = 2  RELATIVE ERROR BETWEEN TWO CONSECUTIVE ITERATES
!                   IS AT MOST XTOL.
!
!         INFO = 3  CONDITIONS FOR INFO = 1 AND INFO = 2 BOTH HOLD.
!
!         INFO = 4  THE COSINE OF THE ANGLE BETWEEN FVEC AND ANY
!                   COLUMN OF THE JACOBIAN IS AT MOST GTOL IN
!                   ABSOLUTE VALUE.
!
!         INFO = 5  NUMBER OF CALLS TO FCN HAS REACHED OR
!                   EXCEEDED MAXFEV.
!
!         INFO = 6  FTOL IS TOO SMALL. NO FURTHER REDUCTION IN
!                   THE SUM OF SQUARES IS POSSIBLE.
!
!         INFO = 7  XTOL IS TOO SMALL. NO FURTHER IMPROVEMENT IN
!                   THE APPROXIMATE SOLUTION X IS POSSIBLE.
!
!         INFO = 8  GTOL IS TOO SMALL. FVEC IS ORTHOGONAL TO THE
!                   COLUMNS OF THE JACOBIAN TO MACHINE PRECISION.
!
!         INFO = 9  THE TWO-NORM OF FVEC IS SMALLER THAN FATOL.
!
!       NFEV IS AN INTEGER OUTPUT VARIABLE SET TO THE NUMBER OF
!         CALLS TO FCN.
!
!       FJAC IS AN OUTPUT M BY N ARRAY. THE UPPER N BY N SUBMATRIX
!         OF FJAC CONTAINS AN UPPER TRIANGULAR MATRIX R WITH
!         DIAGONAL ELEMENTS OF NONINCREASING MAGNITUDE SUCH THAT
!
!                T     T           T
!               P *(JAC *JAC)*P = R *R,
!
!         WHERE P IS A PERMUTATION MATRIX AND JAC IS THE FINAL
!         CALCULATED JACOBIAN. COLUMN J OF P IS COLUMN IPVT(J)
!         (SEE BELOW) OF THE IDENTITY MATRIX. THE LOWER TRAPEZOIDAL
!         PART OF FJAC CONTAINS INFORMATION GENERATED DURING
!         THE COMPUTATION OF R.
!
!       LDFJAC IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN M
!         WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY FJAC.
!
!       IPVT IS AN INTEGER OUTPUT ARRAY OF LENGTH N. IPVT
!         DEFINES A PERMUTATION MATRIX P SUCH THAT JAC*P = Q*R,
!         WHERE JAC IS THE FINAL CALCULATED JACOBIAN, Q IS
!         ORTHOGONAL (NOT STORED), AND R IS UPPER TRIANGULAR
!         WITH DIAGONAL ELEMENTS OF NONINCREASING MAGNITUDE.
!         COLUMN J OF P IS COLUMN IPVT(J) OF THE IDENTITY MATRIX.
!
!       QTF IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS
!         THE FIRST N ELEMENTS OF THE VECTOR (Q TRANSPOSE)*FVEC.
!
!       WA1, WA2, AND WA3 ARE WORK ARRAYS OF LENGTH N.
!
!       WA4 IS A WORK ARRAY OF LENGTH M.
!
!     SUBPROGRAMS CALLED
!
!       USER-SUPPLIED ...... FCN
!
!       MINPACK-SUPPLIED ... DPMPAR,ENORM,FDJAC2,LMPAR,QRFAC
!
!       FORTRAN-SUPPLIED ... DABS,DMAX1,DMIN1,DSQRT,MOD
!
!     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
!
!     **********
      INTEGER I,IFLAG,ITER,J,L,ID, NV, NFV
      real(dp) ACTRED,DELTA,DIRDER,EPSMCH,FNORM,FNORM1,GNORM, &
           ONE,PAR,PNORM,PRERED,P1,P5,P25,P75,P0001,RATIO, &
           SUM,TEMP,TEMP1,TEMP2,XNORM,ZERO
!      real(dp) DPMPAR,ENORM
      real(dp), allocatable, save :: TEMPX(:),TEMPF(:)
      integer, save :: old_nv = -1
      integer, save :: old_nfv = -1

      logical, optional :: reset_flag
      logical error

      SAVE I,IFLAG,ITER,J,L,ID, &
           ACTRED,DELTA,DIRDER,EPSMCH,FNORM,FNORM1,GNORM, &
           ONE,PAR,PNORM,PRERED,P1,P5,P25,P75,P0001,RATIO, &
           SUM,TEMP,TEMP1,TEMP2,XNORM,ZERO
      DATA ONE,P1,P5,P25,P75,P0001,ZERO &
          /1.0D+0,1.0D-1,5.0D-1,2.5D-1,7.5D-1,1.0D-4,0.0D+0/
!
!     EPSMCH IS THE MACHINE PRECISION.
!
!-----------------------------------------------------------------
!-----------------------------------------------------------------

      error = .false.

!     INITIALIZE SAVED VARIABLES
      if(present(reset_flag)) then
         I=0
         IFLAG=0
         ITER=0
         J=0
         L=0
         ID=0
         ACTRED=0.D+0
         DELTA=0.D+0
         DIRDER=0.D+0
         EPSMCH=0.D+0
         FNORM=0.D+0
         FNORM1=0.D+0
         GNORM=0.D+0
         ONE=1.0D+0
         PAR=0.D+0
         PNORM=0.D+0
         PRERED=0.D+0
         P1=1.D-1
         P5=5.D-1
         P25=2.5D-1
         P75=7.5D-1
         P0001=1.D-4
         RATIO=0.D+0
         SUM=0.D+0
         TEMP=0.D+0
         TEMP1=0.D+0
         TEMP2=0.D+0
         XNORM=0.D+0
         ZERO=0.D+0
         old_nv=-1
         old_nfv=-1
         CALL FDJAC2(M,N,X,FVEC,FJAC,LDFJAC,IFLAG,EPSFCN,WA4,ID,IPOS1,.TRUE.)
         return
      endif

      nv = size (x)
      nfv = size(fvec)
      if((old_nv.ne.nv).or.(old_nfv.ne.nfv)) then
         old_nv=nv
         old_nfv=nfv
         if(allocated(tempx)) deallocate(tempx)
         if(allocated(tempf)) deallocate(tempf)
         allocate(tempx(nv))
         allocate(tempf(nfv))
         tempx=0.D+0
         tempf=0.D+0
      endif

      EPSMCH = DPMPAR(1)
      select case (ipos)
      case (1); goto 1001
      case (2); goto 1002
      case (3); goto 1003
      case (4); goto 1004
      case (5); goto 1005
      end select
!
      INFO = 0
      IFLAG = 0
      NFEV = 0
!
!     CHECK THE INPUT PARAMETERS FOR ERRORS.
!
      IF (N .LE. 0 .OR. M .LT. N .OR. LDFJAC .LT. M &
         .OR. FTOL .LT. ZERO .OR. XTOL .LT. ZERO .OR. GTOL .LT. ZERO &
         .OR. MAXFEV .LE. 0 .OR. FACTOR .LE. ZERO) GO TO 300
      IF (MODE .NE. 2) GO TO 20
      DO 10 J = 1, N
         IF (DIAG(J) .LE. ZERO) GO TO 300
   10    CONTINUE
   20 CONTINUE
!
!     EVALUATE THE FUNCTION AT THE STARTING POINT
!     AND CALCULATE ITS NORM.
!
      IFLAG = 1
      IPOS = 1
      IEND = 0
      RETURN
1001  NFEV = 1
      IF (IFLAG .LT. 0) GO TO 300
      FNORM = ENORM(M,FVEC)
!
!     INITIALIZE LEVENBERG-MARQUARDT PARAMETER AND ITERATION COUNT
!
      PAR = ZERO
      ITER = 1
!
!     BEGINNING OF THE OUTER LOOP.
!
   30 CONTINUE
!
!        CALCULATE THE JACOBIAN MATRIX.
!
         IFLAG = 2
         IPOS = 4
   IEND = 0
 1004    CALL FDJAC2(M,N,X,FVEC,FJAC,LDFJAC,IFLAG,EPSFCN,WA4,ID,IPOS1)
   IF (ID .EQ. 1) THEN
      IEND = 0
      RETURN
   ENDIF
         NFEV = NFEV + N
         IF (IFLAG .LT. 0) GO TO 300
!
!        IF REQUESTED, CALL FCN TO ENABLE PRINTING OF ITERATES.
!
         IF (NPRINT .LE. 0) GO TO 40
         IFLAG = 0
         IF (MOD(ITER-1,NPRINT) .EQ. 0) THEN
            IPOS = 2
            IEND = 0
            NFEV = NFEV + 1
            RETURN
         ENDIF
 1002    IF (IFLAG .LT. 0) GO TO 300
   40    CONTINUE
!
!        COMPUTE THE QR FACTORIZATION OF THE JACOBIAN.
!
         CALL QRFAC(M,N,FJAC,LDFJAC,.TRUE.,IPVT,N,WA1,WA2,WA3)
!
!        ON THE FIRST ITERATION AND IF MODE IS 1, SCALE ACCORDING
!        TO THE NORMS OF THE COLUMNS OF THE INITIAL JACOBIAN.
!
         IF (ITER .NE. 1) GO TO 80
         IF (MODE .EQ. 2) GO TO 60
         DO 50 J = 1, N
            DIAG(J) = WA2(J)
            IF (WA2(J) .EQ. ZERO) DIAG(J) = ONE
   50       CONTINUE
   60    CONTINUE
!
!     ON THE FIRST ITERATION, CALCULATE THE NORM OF THE SCALED X
!     AND INITIALIZE THE STEP BOUND DELTA.
!     
         DO 70 J = 1, N
            WA3(J) = DIAG(J)*X(J)
   70       CONTINUE
         XNORM = ENORM(N,WA3)
         DELTA = FACTOR*XNORM
         IF (DELTA .EQ. ZERO) DELTA = FACTOR
   80    CONTINUE
!
!        FORM (Q TRANSPOSE)*FVEC AND STORE THE FIRST N COMPONENTS IN
!        QTF.
!
         DO 90 I = 1, M
            WA4(I) = FVEC(I)
   90       CONTINUE
         DO 130 J = 1, N
            IF (FJAC(J,J) .EQ. ZERO) GO TO 120
            SUM = ZERO
            DO 100 I = J, M
               SUM = SUM + FJAC(I,J)*WA4(I)
  100          CONTINUE
            TEMP = -SUM/FJAC(J,J)
            DO 110 I = J, M
               WA4(I) = WA4(I) + FJAC(I,J)*TEMP
  110          CONTINUE
  120       CONTINUE
            FJAC(J,J) = WA1(J)
            QTF(J) = WA4(J)
  130       CONTINUE
!
!        COMPUTE THE NORM OF THE SCALED GRADIENT.
!
         GNORM = ZERO
         IF (FNORM .EQ. ZERO) GO TO 170
         DO 160 J = 1, N
            L = IPVT(J)
            IF (WA2(L) .EQ. ZERO) GO TO 150
            SUM = ZERO
            DO 140 I = 1, J
               SUM = SUM + FJAC(I,J)*(QTF(I)/FNORM)
  140          CONTINUE
            GNORM = DMAX1(GNORM,DABS(SUM/WA2(L)))
  150       CONTINUE
  160       CONTINUE
  170    CONTINUE
!
!        TEST FOR CONVERGENCE OF THE GRADIENT NORM.
!
!         IF (GNORM .LE. GTOL) INFO = 4
!         IF (INFO .NE. 0) GO TO 300
!       
!        RESCALE IF NECESSARY.
!
         IF (MODE .EQ. 2) GO TO 190
         DO 180 J = 1, N
            DIAG(J) = DMAX1(DIAG(J),WA2(J))
  180       CONTINUE
  190    CONTINUE
!
!        BEGINNING OF THE INNER LOOP.
!
  200    CONTINUE
!
!           DETERMINE THE LEVENBERG-MARQUARDT PARAMETER.
!       
            CALL LMPAR(N,FJAC,LDFJAC,IPVT,DIAG,QTF,DELTA,PAR,WA1,WA2, &
                      WA3,WA4, error)
            if (error) return
!
!           STORE THE DIRECTION P AND X + P. CALCULATE THE NORM OF P.
!
            DO 210 J = 1, N
               WA1(J) = -WA1(J)
               WA2(J) = X(J) + WA1(J)
               WA3(J) = DIAG(J)*WA1(J)
  210          CONTINUE
            PNORM = ENORM(N,WA3)
!
!           ON THE FIRST ITERATION, ADJUST THE INITIAL STEP BOUND.
!
            IF (ITER .EQ. 1) DELTA = DMIN1(DELTA,PNORM)
!
!           EVALUATE THE FUNCTION AT X + P AND CALCULATE ITS NORM.
!
            IFLAG = 1
      DO 310 J = 1, N
      TEMPX(J) = X(J)
      X(J) = WA2(J)
  310      CONTINUE
      DO 320 I = 1, M
      TEMPF(I) = FVEC(I)
  320      CONTINUE         
      IPOS = 3
      IEND = 0
      RETURN
 1003      DO 330 J = 1, N
      X(J) = TEMPX(J)
  330      CONTINUE
      DO 340 I = 1, M
      WA4(I) = FVEC(I)
      FVEC(I) = TEMPF(I)
  340      CONTINUE
            NFEV = NFEV + 1
            IF (IFLAG .LT. 0) GO TO 300
            FNORM1 = ENORM(M,WA4)
!
!           COMPUTE THE SCALED ACTUAL REDUCTION.
!
            ACTRED = -ONE
            IF (P1*FNORM1 .LT. FNORM) ACTRED = ONE - (FNORM1/FNORM)**2
!
!           COMPUTE THE SCALED PREDICTED REDUCTION AND
!           THE SCALED DIRECTIONAL DERIVATIVE.
!
            DO 230 J = 1, N
               WA3(J) = ZERO
               L = IPVT(J)
               TEMP = WA1(L)
               DO 220 I = 1, J
                  WA3(I) = WA3(I) + FJAC(I,J)*TEMP
  220             CONTINUE
  230          CONTINUE
            TEMP1 = ENORM(N,WA3)/FNORM
            TEMP2 = (DSQRT(PAR)*PNORM)/FNORM
            PRERED = TEMP1**2 + TEMP2**2/P5
            DIRDER = -(TEMP1**2 + TEMP2**2)
!
!           COMPUTE THE RATIO OF THE ACTUAL TO THE PREDICTED
!           REDUCTION.
!
            RATIO = ZERO
            IF (PRERED .NE. ZERO) RATIO = ACTRED/PRERED
!
!           UPDATE THE STEP BOUND.
!
            IF (RATIO .GT. P25) GO TO 240
               IF (ACTRED .GE. ZERO) TEMP = P5
               IF (ACTRED .LT. ZERO) &
                 TEMP = P5*DIRDER/(DIRDER + P5*ACTRED)
               IF (P1*FNORM1 .GE. FNORM .OR. TEMP .LT. P1) TEMP = P1
               DELTA = TEMP*DMIN1(DELTA,PNORM/P1)
               PAR = PAR/TEMP
               GO TO 260
  240       CONTINUE
               IF (PAR .NE. ZERO .AND. RATIO .LT. P75) GO TO 250
               DELTA = PNORM/P5
               PAR = P5*PAR
  250          CONTINUE
  260       CONTINUE
!
!           TEST FOR SUCCESSFUL ITERATION.
!
            IF (RATIO .LT. P0001) GO TO 290
!
!           SUCCESSFUL ITERATION. UPDATE X, FVEC, AND THEIR NORMS.
!
            DO 270 J = 1, N
               X(J) = WA2(J)
               WA2(J) = DIAG(J)*X(J)
  270          CONTINUE
            DO 280 I = 1, M
               FVEC(I) = WA4(I)
  280          CONTINUE
            XNORM = ENORM(N,WA2)
            FNORM = FNORM1
            ITER = ITER + 1
  290       CONTINUE
!
!           TESTS FOR CONVERGENCE.
!
            IF (DABS(ACTRED) .LE. FTOL .AND. PRERED .LE. FTOL &
               .AND. P5*RATIO .LE. ONE) INFO = 1
!            IF (DELTA .LE. XTOL*XNORM) INFO = 2
!            IF (DABS(ACTRED) .LE. FTOL .AND. PRERED .LE. FTOL
!     *          .AND. P5*RATIO .LE. ONE .AND. INFO .EQ. 2) INFO = 3
            IF (FNORM .LE. FATOL) INFO = 9
            IF (INFO .NE. 0) GO TO 300
!
!           TESTS FOR TERMINATION AND STRINGENT TOLERANCES.
!
            IF (NFEV .GE. MAXFEV) INFO = 5
            IF (DABS(ACTRED) .LE. EPSMCH .AND. PRERED .LE. EPSMCH &
               .AND. P5*RATIO .LE. ONE) INFO = 6
!            IF (DELTA .LE. EPSMCH*XNORM) INFO = 7
!            IF (GNORM .LE. EPSMCH) INFO = 8
            IF (INFO .NE. 0) GO TO 300
!
!           END OF THE INNER LOOP. REPEAT IF ITERATION UNSUCCESSFUL.
!
            IF (RATIO .LT. P0001) GO TO 200
!
!        END OF THE OUTER LOOP.
!
         GO TO 30
  300 CONTINUE
!
!     TERMINATION, EITHER NORMAL OR USER IMPOSED.
!
      IF (IFLAG .LT. 0) INFO = IFLAG
      IFLAG = 0
      IF (NPRINT .GT. 0) THEN
         IPOS = 5
         IEND = 0
         NFEV = NFEV + 1
         RETURN
      ENDIF
 1005 IEND = 1          
      IPOS = 0
      IPOS1 = 0
      RETURN
!
!     LAST CARD OF SUBROUTINE LMDIF.
!
      END SUBROUTINE

!     **************************************
      real(dp) FUNCTION ENORM(N,X)
      INTEGER N
      real(dp) X(N)
!     ***************************************
!
!     FUNCTION ENORM
!
!     GIVEN AN N-VECTOR X, THIS FUNCTION CALCULATES THE
!     EUCLIDEAN NORM OF X.
!
!     THE EUCLIDEAN NORM IS COMPUTED BY ACCUMULATING THE SUM OF
!     SQUARES IN THREE DIFFERENT SUMS. THE SUMS OF SQUARES FOR THE
!     SMALL AND LARGE COMPONENTS ARE SCALED SO THAT NO OVERFLOWS
!     OCCUR. NON-DESTRUCTIVE UNDERFLOWS ARE PERMITTED. UNDERFLOWS
!     AND OVERFLOWS DO NOT OCCUR IN THE COMPUTATION OF THE UNSCALED
!     SUM OF SQUARES FOR THE INTERMEDIATE COMPONENTS.
!     THE DEFINITIONS OF SMALL, INTERMEDIATE AND LARGE COMPONENTS
!     DEPEND ON TWO CONSTANTS, RDWARF AND RGIANT. THE MAIN
!     RESTRICTIONS ON THESE CONSTANTS ARE THAT RDWARF**2 NOT
!     UNDERFLOW AND RGIANT**2 NOT OVERFLOW. THE CONSTANTS
!     GIVEN HERE ARE SUITABLE FOR EVERY KNOWN COMPUTER.
!
!     THE FUNCTION STATEMENT IS
!
!       real(dp) FUNCTION ENORM(N,X)
!
!     WHERE
!
!       N IS A POSITIVE INTEGER INPUT VARIABLE.
!
!       X IS AN INPUT ARRAY OF LENGTH N.
!
!     SUBPROGRAMS CALLED
!
!       FORTRAN-SUPPLIED ... DABS,DSQRT
!
!     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
!
!     **********
      INTEGER I
      real(dp) AGIANT,FLOATN,ONE,RDWARF,RGIANT,S1,S2,S3,XABS, &
                      X1MAX,X3MAX,ZERO
      DATA ONE,ZERO,RDWARF,RGIANT /1.0D0,0.0D0,3.834D-20,1.304D19/
      S1 = ZERO
      S2 = ZERO
      S3 = ZERO
      X1MAX = ZERO
      X3MAX = ZERO
      FLOATN = N
      AGIANT = RGIANT/FLOATN
      DO 90 I = 1, N
         XABS = DABS(X(I))
         IF (XABS .GT. RDWARF .AND. XABS .LT. AGIANT) GO TO 70
            IF (XABS .LE. RDWARF) GO TO 30
!
!              SUM FOR LARGE COMPONENTS.
!
               IF (XABS .LE. X1MAX) GO TO 10
                  S1 = ONE + S1*(X1MAX/XABS)**2
                  X1MAX = XABS
                  GO TO 20
   10          CONTINUE
                  S1 = S1 + (XABS/X1MAX)**2
   20          CONTINUE
               GO TO 60
   30       CONTINUE
!
!              SUM FOR SMALL COMPONENTS.
!
               IF (XABS .LE. X3MAX) GO TO 40
                  S3 = ONE + S3*(X3MAX/XABS)**2
                  X3MAX = XABS
                  GO TO 50
   40          CONTINUE
                  IF (XABS .NE. ZERO) S3 = S3 + (XABS/X3MAX)**2
   50          CONTINUE
   60       CONTINUE
            GO TO 80
   70    CONTINUE
!
!           SUM FOR INTERMEDIATE COMPONENTS.
!
            S2 = S2 + XABS**2
   80    CONTINUE
   90    CONTINUE
!
!     CALCULATION OF NORM.
!
      IF (S1 .EQ. ZERO) GO TO 100
         ENORM = X1MAX*DSQRT(S1+(S2/X1MAX)/X1MAX)
         GO TO 130
  100 CONTINUE
         IF (S2 .EQ. ZERO) GO TO 110
            IF (S2 .GE. X3MAX) &
              ENORM = DSQRT(S2*(ONE+(X3MAX/S2)*(X3MAX*S3)))
            IF (S2 .LT. X3MAX) &
              ENORM = DSQRT(X3MAX*((S2/X3MAX)+(X3MAX*S3)))
            GO TO 120
  110    CONTINUE
            ENORM = X3MAX*DSQRT(S3)
  120    CONTINUE
  130 CONTINUE
      RETURN
!
!     LAST CARD OF FUNCTION ENORM.
!
      END FUNCTION
!     **************************************************************
  SUBROUTINE FDJAC2(M,N,X,FVEC,FJAC,LDFJAC,IFLAG,EPSFCN,WA,ID,IPOS,reset_flag)
      INTEGER M,N,LDFJAC,IFLAG,ID,IPOS
      real(dp) EPSFCN
      real(dp) X(N),FVEC(M),FJAC(LDFJAC,N),WA(M)
!     ***************************************************************
!
!     SUBROUTINE FDJAC2
!
!     THIS SUBROUTINE COMPUTES A FORWARD-DIFFERENCE APPROXIMATION
!     TO THE M BY N JACOBIAN MATRIX ASSOCIATED WITH A SPECIFIED
!     PROBLEM OF M FUNCTIONS IN N VARIABLES.
!
!     THE SUBROUTINE STATEMENT IS
!
!       SUBROUTINE FDJAC2(FCN,M,N,X,FVEC,FJAC,LDFJAC,IFLAG,EPSFCN,WA)
!
!     WHERE
!
!       FCN IS THE NAME OF THE USER-SUPPLIED SUBROUTINE WHICH
!         CALCULATES THE FUNCTIONS. FCN MUST BE DECLARED
!         IN AN EXTERNAL STATEMENT IN THE USER CALLING
!         PROGRAM, AND SHOULD BE WRITTEN AS FOLLOWS.
!
!         SUBROUTINE FCN(M,N,X,FVEC,IFLAG)
!         INTEGER M,N,IFLAG
!         real(dp) X(N),FVEC(M)
!         ----------
!         CALCULATE THE FUNCTIONS AT X AND
!         RETURN THIS VECTOR IN FVEC.
!         ----------
!         RETURN
!         END
!
!         THE VALUE OF IFLAG SHOULD NOT BE CHANGED BY FCN UNLESS
!         THE USER WANTS TO TERMINATE EXECUTION OF FDJAC2.
!         IN THIS CASE SET IFLAG TO A NEGATIVE INTEGER.
!
!       M IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!         OF FUNCTIONS.
!
!       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!         OF VARIABLES. N MUST NOT EXCEED M.
!
!       X IS AN INPUT ARRAY OF LENGTH N.
!
!       FVEC IS AN INPUT ARRAY OF LENGTH M WHICH MUST CONTAIN THE
!         FUNCTIONS EVALUATED AT X.
!
!       FJAC IS AN OUTPUT M BY N ARRAY WHICH CONTAINS THE
!         APPROXIMATION TO THE JACOBIAN MATRIX EVALUATED AT X.
!
!       LDFJAC IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN M
!         WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY FJAC.
!
!       IFLAG IS AN INTEGER VARIABLE WHICH CAN BE USED TO TERMINATE
!         THE EXECUTION OF FDJAC2. SEE DESCRIPTION OF FCN.
!
!       EPSFCN IS AN INPUT VARIABLE USED IN DETERMINING A SUITABLE
!         STEP LENGTH FOR THE FORWARD-DIFFERENCE APPROXIMATION. THIS
!         APPROXIMATION ASSUMES THAT THE RELATIVE ERRORS IN THE
!         FUNCTIONS ARE OF THE ORDER OF EPSFCN. IF EPSFCN IS LESS
!         THAN THE MACHINE PRECISION, IT IS ASSUMED THAT THE RELATIVE
!         ERRORS IN THE FUNCTIONS ARE OF THE ORDER OF THE MACHINE
!         PRECISION.
!
!       WA IS A WORK ARRAY OF LENGTH M.
!
!     SUBPROGRAMS CALLED
!
!       USER-SUPPLIED ...... FCN
!
!       MINPACK-SUPPLIED ... DPMPAR
!
!       FORTRAN-SUPPLIED ... DABS,DMAX1,DSQRT
!
!     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
!
!     **********
!     EXTERNAL FCN
      INTEGER I,J, NFV
      real(dp) EPS,EPSMCH,H,TEMP,ZERO
!      real(dp) DPMPAR
      real(dp), allocatable, save :: TEMPF(:)
      integer, save :: old_nfv = -1

      logical, optional :: reset_flag

      SAVE I,J,EPS,EPSMCH,H,TEMP,ZERO
      DATA ZERO /0.0D0/

      if(present(reset_flag)) then
        I=0
        J=0
        EPS=0.D+0
        EPSMCH=0.D+0
        H=0.D+0
        TEMP=0.D+0
        ZERO=0.D+0
        old_nfv=-1
        return
      endif

      nfv = size(fvec)
      if(old_nfv.ne.nfv) then
         old_nfv=nfv
         if(allocated(tempf)) deallocate(tempf)
         allocate(tempf(nfv))
         tempf=0.D+0
      endif
!
!     EPSMCH IS THE MACHINE PRECISION.
!
      EPSMCH = DPMPAR(1)
      IF (IPOS .EQ. 4) GO TO 104  
!
      EPS = DSQRT(DMAX1(EPSFCN,EPSMCH))
    J = 0
  20   J = J + 1
         TEMP = X(J)
         H = EPS*DABS(TEMP)
         IF (H .EQ. ZERO) H = EPS
         X(J) = TEMP + H
   DO 40 I = 1, M
   TEMPF(I) = FVEC(I)
  40   CONTINUE
   IPOS = 4
   ID = 1
   RETURN
 104    DO 50 I = 1, M
   WA(I) = FVEC(I)
   FVEC(I) = TEMPF(I)
  50    CONTINUE  
         IF (IFLAG .LT. 0) GO TO 30
         X(J) = TEMP
         DO 10 I = 1, M
            FJAC(I,J) = (WA(I) - FVEC(I))/H
   10       CONTINUE
  IF (J .LT. N) GO TO 20
   30 CONTINUE
      ID = 0  
      IPOS = 0
      RETURN
!
!     LAST CARD OF SUBROUTINE FDJAC2.
!
      END SUBROUTINE
!     *************************************
      real(dp) FUNCTION DPMPAR(I)
      INTEGER I
!     **********
!
!     FUNCTION DPMPAR
!    
!     THIS FUNCTION PROVIDES real(dp) MACHINE PARAMETERS
!     WHEN THE APPROPRIATE SET OF DATA STATEMENTS IS ACTIVATED (BY
!     REMOVING THE C FROM COLUMN 1) AND ALL OTHER DATA STATEMENTS ARE
!     RENDERED INACTIVE. MOST OF THE PARAMETER VALUES WERE OBTAINED
!     FROM THE CORRESPONDING BELL LABORATORIES PORT LIBRARY FUNCTION.
!
!     THE FUNCTION STATEMENT IS
!
!       real(dp) FUNCTION DPMPAR(I)
!
!     WHERE
!
!       I IS AN INTEGER INPUT VARIABLE SET TO 1, 2, OR 3 WHICH
!         SELECTS THE DESIRED MACHINE PARAMETER. IF THE MACHINE HAS
!         T BASE B DIGITS AND ITS SMALLEST AND LARGEST EXPONENTS ARE
!         EMIN AND EMAX, RESPECTIVELY, THEN THESE PARAMETERS ARE
!
!         DPMPAR(1) = B**(1 - T), THE MACHINE PRECISION,
!
!         DPMPAR(2) = B**(EMIN - 1), THE SMALLEST MAGNITUDE,
!
!         DPMPAR(3) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
!
!     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983.
!     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
!
!     **********
      INTEGER MCHEPS(4)
      INTEGER MINMAG(4)
      INTEGER MAXMAG(4)
      real(dp) DMACH(3)
      EQUIVALENCE (DMACH(1),MCHEPS(1))
      EQUIVALENCE (DMACH(2),MINMAG(1))
      EQUIVALENCE (DMACH(3),MAXMAG(1))
!
!     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
!     THE AMDAHL 470/V6, THE ICL 2900, THE ITEL AS/6,
!     THE XEROX SIGMA 5/7/9 AND THE SEL SYSTEMS 85/86.
!
!        DATA MCHEPS(1),MCHEPS(2) / Z'34100000', Z'00000000' /
!        DATA MINMAG(1),MINMAG(2) / Z'00100000', Z'00000000' /
!        DATA MAXMAG(1),MAXMAG(2) / Z'7FFFFFFF', Z'FFFFFFFF' /
!
!    My IBM constants as of July 1989
!      DATA DMACH /1.0D-15,1.0D-78,1.0D75/
!    My VAX constants as of July 1989
!      DATA DMACH /1.0D-15,1.0D-38,1.0D38/
!    My DEC/ULTRIX constants as of July 1989
      DATA DMACH /1.0D-15,1.0D-300,1.0D300/
!
!     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES.
!
!     DATA MCHEPS(1),MCHEPS(2) / O606400000000, O000000000000 /
!     DATA MINMAG(1),MINMAG(2) / O402400000000, O000000000000 /
!     DATA MAXMAG(1),MAXMAG(2) / O376777777777, O777777777777 /
!
!     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES.
!
!     DATA MCHEPS(1) / 15614000000000000000B /
!     DATA MCHEPS(2) / 15010000000000000000B /
!
!     DATA MINMAG(1) / 00604000000000000000B /
!     DATA MINMAG(2) / 00000000000000000000B /
!
!     DATA MAXMAG(1) / 37767777777777777777B /
!     DATA MAXMAG(2) / 37167777777777777777B /
!
!     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR).
!
!     DATA MCHEPS(1),MCHEPS(2) / "114400000000, "000000000000 /
!     DATA MINMAG(1),MINMAG(2) / "033400000000, "000000000000 /
!     DATA MAXMAG(1),MAXMAG(2) / "377777777777, "344777777777 /
!
!     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR).
!
!    DATA MCHEPS(1),MCHEPS(2) / "104400000000, "000000000000 /
!     DATA MINMAG(1),MINMAG(2) / "000400000000, "000000000000 /
!     DATA MAXMAG(1),MAXMAG(2) / "377777777777, "377777777777 /
!
!     MACHINE CONSTANTS FOR THE PDP-11.
!
!     DATA MCHEPS(1),MCHEPS(2) /   9472,      0 /
!     DATA MCHEPS(3),MCHEPS(4) /      0,      0 /
!
!     DATA MINMAG(1),MINMAG(2) /    128,      0 /
!     DATA MINMAG(3),MINMAG(4) /      0,      0 /
!
!     DATA MAXMAG(1),MAXMAG(2) /  32767,     -1 /
!     DATA MAXMAG(3),MAXMAG(4) /     -1,     -1 /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS.
!
!     DATA MCHEPS(1) / O1451000000000000 /
!     DATA MCHEPS(2) / O0000000000000000 /
!
!     DATA MINMAG(1) / O1771000000000000 /
!     DATA MINMAG(2) / O7770000000000000 /
!
!     DATA MAXMAG(1) / O0777777777777777 /
!     DATA MAXMAG(2) / O7777777777777777 /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM.
!
!     DATA MCHEPS(1) / O1451000000000000 /
!     DATA MCHEPS(2) / O0000000000000000 /
!
!     DATA MINMAG(1) / O1771000000000000 /
!     DATA MINMAG(2) / O0000000000000000 /
!
!     DATA MAXMAG(1) / O0777777777777777 /
!     DATA MAXMAG(2) / O0007777777777777 /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
!
!     DATA MCHEPS(1) / ZCC6800000 /
!     DATA MCHEPS(2) / Z000000000 /
!
!     DATA MINMAG(1) / ZC00800000 /
!     DATA MINMAG(2) / Z000000000 /
!
!     DATA MAXMAG(1) / ZDFFFFFFFF /
!     DATA MAXMAG(2) / ZFFFFFFFFF /
!
!     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
!
!     DATA MCHEPS(1),MCHEPS(2) / O170640000000, O000000000000 /
!     DATA MINMAG(1),MINMAG(2) / O000040000000, O000000000000 /
!     DATA MAXMAG(1),MAXMAG(2) / O377777777777, O777777777777 /
!
!     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200.
!
!     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING CARD -
!     STATIC DMACH(3)
!
!     DATA MINMAG/20K,3*0/,MAXMAG/77777K,3*177777K/
!     DATA MCHEPS/32020K,3*0/
!
!     MACHINE CONSTANTS FOR THE HARRIS 220.
!
!     DATA MCHEPS(1),MCHEPS(2) / '20000000, '00000334 /
!     DATA MINMAG(1),MINMAG(2) / '20000000, '00000201 /
!     DATA MAXMAG(1),MAXMAG(2) / '37777777, '37777577 /
!
!     MACHINE CONSTANTS FOR THE CRAY-1.
!
!     DATA MCHEPS(1) / 0376424000000000000000B /
!     DATA MCHEPS(2) / 0000000000000000000000B /
!
!     DATA MINMAG(1) / 0200034000000000000000B /
!     DATA MINMAG(2) / 0000000000000000000000B /
!
!     DATA MAXMAG(1) / 0577777777777777777777B /
!     DATA MAXMAG(2) / 0000007777777777777776B /
!
!     MACHINE CONSTANTS FOR THE PRIME 400.
!
!     DATA MCHEPS(1),MCHEPS(2) / :10000000000, :00000000123 /
!     DATA MINMAG(1),MINMAG(2) / :10000000000, :00000100000 /
!     DATA MAXMAG(1),MAXMAG(2) / :17777777777, :37777677776 /
!
!     MACHINE CONSTANTS FOR THE VAX-11.
!
!     DATA MCHEPS(1),MCHEPS(2) /   9472,  0 /
!     DATA MINMAG(1),MINMAG(2) /    128,  0 /
!     DATA MAXMAG(1),MAXMAG(2) / -32769, -1 /
!
      DPMPAR = DMACH(I)
      RETURN
!
!     LAST CARD OF FUNCTION DPMPAR.
!
      END FUNCTION
!     **************************************************************
      SUBROUTINE LMPAR(N,R,LDR,IPVT,DIAG,QTB,DELTA,PAR,X,SDIAG,WA1, &
                      WA2, error)
      INTEGER N,LDR
      INTEGER IPVT(N)
      real(dp) DELTA,PAR
      real(dp) R(LDR,N),DIAG(N),QTB(N),X(N),SDIAG(N),WA1(N), &
                      WA2(N)
!     ***************************************************************
!
!     SUBROUTINE LMPAR
!
!     GIVEN AN M BY N MATRIX A, AN N BY N NONSINGULAR DIAGONAL
!     MATRIX D, AN M-VECTOR B, AND A POSITIVE NUMBER DELTA,
!     THE PROBLEM IS TO DETERMINE A VALUE FOR THE PARAMETER
!     PAR SUCH THAT IF X SOLVES THE SYSTEM
!
!           A*X = B ,     SQRT(PAR)*D*X = 0 ,
!
!     IN THE LEAST SQUARES SENSE, AND DXNORM IS THE EUCLIDEAN
!     NORM OF D*X, THEN EITHER PAR IS ZERO AND
!
!           (DXNORM-DELTA) .LE. 0.1*DELTA ,
!
!     OR PAR IS POSITIVE AND
!
!           ABS(DXNORM-DELTA) .LE. 0.1*DELTA .
!
!     THIS SUBROUTINE COMPLETES THE SOLUTION OF THE PROBLEM
!     IF IT IS PROVIDED WITH THE NECESSARY INFORMATION FROM THE
!     QR FACTORIZATION, WITH COLUMN PIVOTING, OF A. THAT IS, IF
!     A*P = Q*R, WHERE P IS A PERMUTATION MATRIX, Q HAS ORTHOGONAL
!     COLUMNS, AND R IS AN UPPER TRIANGULAR MATRIX WITH DIAGONAL
!     ELEMENTS OF NONINCREASING MAGNITUDE, THEN LMPAR EXPECTS
!     THE FULL UPPER TRIANGLE OF R, THE PERMUTATION MATRIX P,
!     AND THE FIRST N COMPONENTS OF (Q TRANSPOSE)*B. ON OUTPUT
!     LMPAR ALSO PROVIDES AN UPPER TRIANGULAR MATRIX S SUCH THAT
!
!            T   T                   T
!           P *(A *A + PAR*D*D)*P = S *S .
!
!     S IS EMPLOYED WITHIN LMPAR AND MAY BE OF SEPARATE INTEREST.
!
!     ONLY A FEW ITERATIONS ARE GENERALLY NEEDED FOR CONVERGENCE
!     OF THE ALGORITHM. IF, HOWEVER, THE LIMIT OF 10 ITERATIONS
!     IS REACHED, THEN THE OUTPUT PAR WILL CONTAIN THE BEST
!     VALUE OBTAINED SO FAR.
!
!     THE SUBROUTINE STATEMENT IS
!
!       SUBROUTINE LMPAR(N,R,LDR,IPVT,DIAG,QTB,DELTA,PAR,X,SDIAG,
!                        WA1,WA2, error)
!
!     WHERE
!
!       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE ORDER OF R.
!
!       R IS AN N BY N ARRAY. ON INPUT THE FULL UPPER TRIANGLE
!         MUST CONTAIN THE FULL UPPER TRIANGLE OF THE MATRIX R.
!         ON OUTPUT THE FULL UPPER TRIANGLE IS UNALTERED, AND THE
!         STRICT LOWER TRIANGLE CONTAINS THE STRICT UPPER TRIANGLE
!         (TRANSPOSED) OF THE UPPER TRIANGULAR MATRIX S.
!
!       LDR IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN N
!         WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY R.
!
!       IPVT IS AN INTEGER INPUT ARRAY OF LENGTH N WHICH DEFINES THE
!         PERMUTATION MATRIX P SUCH THAT A*P = Q*R. COLUMN J OF P
!         IS COLUMN IPVT(J) OF THE IDENTITY MATRIX.
!
!       DIAG IS AN INPUT ARRAY OF LENGTH N WHICH MUST CONTAIN THE
!         DIAGONAL ELEMENTS OF THE MATRIX D.
!
!       QTB IS AN INPUT ARRAY OF LENGTH N WHICH MUST CONTAIN THE FIRST
!         N ELEMENTS OF THE VECTOR (Q TRANSPOSE)*B.
!
!       DELTA IS A POSITIVE INPUT VARIABLE WHICH SPECIFIES AN UPPER
!         BOUND ON THE EUCLIDEAN NORM OF D*X.
!
!       PAR IS A NONNEGATIVE VARIABLE. ON INPUT PAR CONTAINS AN
!         INITIAL ESTIMATE OF THE LEVENBERG-MARQUARDT PARAMETER.
!         ON OUTPUT PAR CONTAINS THE FINAL ESTIMATE.
!
!       X IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE LEAST
!         SQUARES SOLUTION OF THE SYSTEM A*X = B, SQRT(PAR)*D*X = 0,
!         FOR THE OUTPUT PAR.
!
!       SDIAG IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE
!         DIAGONAL ELEMENTS OF THE UPPER TRIANGULAR MATRIX S.
!
!       WA1 AND WA2 ARE WORK ARRAYS OF LENGTH N.
!
!       error is a logical set True if there is an error and flase otherwise
!
!     SUBPROGRAMS CALLED
!
!       MINPACK-SUPPLIED ... DPMPAR,ENORM,QRSOLV
!
!       FORTRAN-SUPPLIED ... DABS,DMAX1,DMIN1,DSQRT
!
!     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
!
!     **********
      INTEGER I,ITER,J,JM1,JP1,K,L,NSING
      real(dp) DXNORM,DWARF,FP,GNORM,PARC,PARL,PARU,P1,P001, &
                      SUM,TEMP,ZERO
      logical error
!      real(dp) DPMPAR,ENORM
      DATA P1,P001,ZERO /1.0D-1,1.0D-3,0.0D0/
!
!     DWARF IS THE SMALLEST POSITIVE MAGNITUDE.
!
      DWARF = DPMPAR(2)
      error = .false.
!
!     COMPUTE AND STORE IN X THE GAUSS-NEWTON DIRECTION. IF THE
!     JACOBIAN IS RANK-DEFICIENT, OBTAIN A LEAST SQUARES SOLUTION.
!
      NSING = N
      DO 10 J = 1, N
         WA1(J) = QTB(J)
         IF (R(J,J) .EQ. ZERO .AND. NSING .EQ. N) NSING = J - 1
         IF (NSING .LT. N) WA1(J) = ZERO
   10    CONTINUE
      IF (NSING .LT. 1) GO TO 50
      DO 40 K = 1, NSING
         J = NSING - K + 1
         WA1(J) = WA1(J)/R(J,J)
         TEMP = WA1(J)
         JM1 = J - 1
         IF (JM1 .LT. 1) GO TO 30
         DO 20 I = 1, JM1
            WA1(I) = WA1(I) - R(I,J)*TEMP
   20       CONTINUE
   30    CONTINUE
   40    CONTINUE
   50 CONTINUE
      DO 60 J = 1, N
         L = IPVT(J)
         X(L) = WA1(J)
   60    CONTINUE
!
!     INITIALIZE THE ITERATION COUNTER.
!     EVALUATE THE FUNCTION AT THE ORIGIN, AND TEST
!     FOR ACCEPTANCE OF THE GAUSS-NEWTON DIRECTION.
!
      ITER = 0
      DO 70 J = 1, N
         WA2(J) = DIAG(J)*X(J)
   70    CONTINUE
      DXNORM = ENORM(N,WA2)
      FP = DXNORM - DELTA
      IF (FP .LE. P1*DELTA) GO TO 220
!
!     IF THE JACOBIAN IS NOT RANK DEFICIENT, THE NEWTON
!     STEP PROVIDES A LOWER BOUND, PARL, FOR THE ZERO OF
!     THE FUNCTION. OTHERWISE SET THIS BOUND TO ZERO.
!
      PARL = ZERO
      IF (NSING .LT. N) GO TO 120
      DO 80 J = 1, N
         L = IPVT(J)
         WA1(J) = DIAG(L)*(WA2(L)/DXNORM)
   80    CONTINUE
      DO 110 J = 1, N
         SUM = ZERO
         JM1 = J - 1
         IF (JM1 .LT. 1) GO TO 100
         DO 90 I = 1, JM1
            SUM = SUM + R(I,J)*WA1(I)
   90       CONTINUE
  100    CONTINUE
         WA1(J) = (WA1(J) - SUM)/R(J,J)
  110    CONTINUE
      TEMP = ENORM(N,WA1)
      PARL = ((FP/DELTA)/TEMP)/TEMP
  120 CONTINUE
!
!     CALCULATE AN UPPER BOUND, PARU, FOR THE ZERO OF THE FUNCTION.
!
      DO 140 J = 1, N
         SUM = ZERO
         DO 130 I = 1, J
            SUM = SUM + R(I,J)*QTB(I)
  130       CONTINUE
         L = IPVT(J)
         WA1(J) = SUM/DIAG(L)
  140    CONTINUE
      GNORM = ENORM(N,WA1)
      PARU = GNORM/DELTA
      IF (PARU .EQ. ZERO) PARU = DWARF/DMIN1(DELTA,P1)
!
!     IF THE INPUT PAR LIES OUTSIDE OF THE INTERVAL (PARL,PARU),
!     SET PAR TO THE CLOSER ENDPOINT.
!
      PAR = DMAX1(PAR,PARL)
      PAR = DMIN1(PAR,PARU)
      IF (PAR .EQ. ZERO) PAR = GNORM/DXNORM
!
!     BEGINNING OF AN ITERATION.
!
  150 CONTINUE
         ITER = ITER + 1
!
!        EVALUATE THE FUNCTION AT THE CURRENT VALUE OF PAR.
!
         IF (PAR .EQ. ZERO) PAR = DMAX1(DWARF,P001*PARU)
         TEMP = DSQRT(PAR)
         DO 160 J = 1, N
            WA1(J) = TEMP*DIAG(J)
  160       CONTINUE
         CALL QRSOLV(N,R,LDR,IPVT,WA1,QTB,X,SDIAG,WA2)
         DO 170 J = 1, N
            WA2(J) = DIAG(J)*X(J)
  170       CONTINUE
         DXNORM = ENORM(N,WA2)
         TEMP = FP
         FP = DXNORM - DELTA
!
!        IF THE FUNCTION IS SMALL ENOUGH, ACCEPT THE CURRENT VALUE
!        OF PAR. ALSO TEST FOR THE EXCEPTIONAL CASES WHERE PARL
!        IS ZERO OR THE NUMBER OF ITERATIONS HAS REACHED 10.
!
         IF (DABS(FP) .LE. P1*DELTA &
            .OR. PARL .EQ. ZERO .AND. FP .LE. TEMP &
                 .AND. TEMP .LT. ZERO .OR. ITER .EQ. 10) GO TO 220
!
!        COMPUTE THE NEWTON CORRECTION.
!
          if (dxnorm == 0) then
            error = .true.
            print *, 'ERROR IN LMDIF: FLOATING UNDERFLOW. WILL STOP HERE.'
            return
          endif

         DO 180 J = 1, N
            L = IPVT(J)
            WA1(J) = DIAG(L)*(WA2(L)/DXNORM)
  180       CONTINUE
         DO 210 J = 1, N
            WA1(J) = WA1(J)/SDIAG(J)
            TEMP = WA1(J)
            JP1 = J + 1
            IF (N .LT. JP1) GO TO 200
            DO 190 I = JP1, N
               WA1(I) = WA1(I) - R(I,J)*TEMP
  190          CONTINUE
  200       CONTINUE
  210       CONTINUE
         TEMP = ENORM(N,WA1)
         PARC = ((FP/DELTA)/TEMP)/TEMP
!
!        DEPENDING ON THE SIGN OF THE FUNCTION, UPDATE PARL OR PARU.
!
         IF (FP .GT. ZERO) PARL = DMAX1(PARL,PAR)
         IF (FP .LT. ZERO) PARU = DMIN1(PARU,PAR)
!
!        COMPUTE AN IMPROVED ESTIMATE FOR PAR.
!
         PAR = DMAX1(PARL,PAR+PARC)
!
!        END OF AN ITERATION.
!
         GO TO 150
  220 CONTINUE
!
!     TERMINATION.
!
      IF (ITER .EQ. 0) PAR = ZERO
      RETURN
!
!     LAST CARD OF SUBROUTINE LMPAR.
!
      END SUBROUTINE
!     ***********************
      SUBROUTINE QRFAC(M,N,A,LDA,PIVOT,IPVT,LIPVT,RDIAG,ACNORM,WA)
      INTEGER M,N,LDA,LIPVT
      INTEGER IPVT(LIPVT)
      LOGICAL PIVOT
      real(dp) A(LDA,N),RDIAG(N),ACNORM(N),WA(N)
!     **************************************************************
!
!     SUBROUTINE QRFAC
!
!     THIS SUBROUTINE USES HOUSEHOLDER TRANSFORMATIONS WITH COLUMN
!     PIVOTING (OPTIONAL) TO COMPUTE A QR FACTORIZATION OF THE
!     M BY N MATRIX A. THAT IS, QRFAC DETERMINES AN ORTHOGONAL
!     MATRIX Q, A PERMUTATION MATRIX P, AND AN UPPER TRAPEZOIDAL
!     MATRIX R WITH DIAGONAL ELEMENTS OF NONINCREASING MAGNITUDE,
!     SUCH THAT A*P = Q*R. THE HOUSEHOLDER TRANSFORMATION FOR
!     COLUMN K, K = 1,2,...,MIN(M,N), IS OF THE FORM
!
!                           T
!           I - (1/U(K))*U*U
!
!     WHERE U HAS ZEROS IN THE FIRST K-1 POSITIONS. THE FORM OF
!     THIS TRANSFORMATION AND THE METHOD OF PIVOTING FIRST
!     APPEARED IN THE CORRESPONDING LINPACK SUBROUTINE.
!
!     THE SUBROUTINE STATEMENT IS
!
!       SUBROUTINE QRFAC(M,N,A,LDA,PIVOT,IPVT,LIPVT,RDIAG,ACNORM,WA)
!
!     WHERE
!
!       M IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!         OF ROWS OF A.
!
!       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!         OF COLUMNS OF A.
!
!       A IS AN M BY N ARRAY. ON INPUT A CONTAINS THE MATRIX FOR
!         WHICH THE QR FACTORIZATION IS TO BE COMPUTED. ON OUTPUT
!         THE STRICT UPPER TRAPEZOIDAL PART OF A CONTAINS THE STRICT
!         UPPER TRAPEZOIDAL PART OF R, AND THE LOWER TRAPEZOIDAL
!         PART OF A CONTAINS A FACTORED FORM OF Q (THE NON-TRIVIAL
!         ELEMENTS OF THE U VECTORS DESCRIBED ABOVE).
!
!       LDA IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN M
!         WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY A.
!
!       PIVOT IS A LOGICAL INPUT VARIABLE. IF PIVOT IS SET TRUE,
!         THEN COLUMN PIVOTING IS ENFORCED. IF PIVOT IS SET FALSE,
!         THEN NO COLUMN PIVOTING IS DONE.
!
!       IPVT IS AN INTEGER OUTPUT ARRAY OF LENGTH LIPVT. IPVT
!         DEFINES THE PERMUTATION MATRIX P SUCH THAT A*P = Q*R.
!         COLUMN J OF P IS COLUMN IPVT(J) OF THE IDENTITY MATRIX.
!         IF PIVOT IS FALSE, IPVT IS NOT REFERENCED.
!
!       LIPVT IS A POSITIVE INTEGER INPUT VARIABLE. IF PIVOT IS FALSE,
!         THEN LIPVT MAY BE AS SMALL AS 1. IF PIVOT IS TRUE, THEN
!         LIPVT MUST BE AT LEAST N.
!
!       RDIAG IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE
!         DIAGONAL ELEMENTS OF R.
!
!       ACNORM IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE
!         NORMS OF THE CORRESPONDING COLUMNS OF THE INPUT MATRIX A.
!         IF THIS INFORMATION IS NOT NEEDED, THEN ACNORM CAN COINCIDE
!         WITH RDIAG.
!
!       WA IS A WORK ARRAY OF LENGTH N. IF PIVOT IS FALSE, THEN WA
!         CAN COINCIDE WITH RDIAG.
!
!     SUBPROGRAMS CALLED
!
!       MINPACK-SUPPLIED ... DPMPAR,ENORM
!
!       FORTRAN-SUPPLIED ... DMAX1,DSQRT,MIN0
!
!     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
!
!     **********
      INTEGER I,J,JP1,K,KMAX,MINMN
      real(dp) AJNORM,EPSMCH,ONE,P05,SUM,TEMP,ZERO
!      real(dp) DPMPAR,ENORM
      DATA ONE,P05,ZERO /1.0D0,5.0D-2,0.0D0/
!
!     EPSMCH IS THE MACHINE PRECISION.
!
      EPSMCH = DPMPAR(1)
!
!     COMPUTE THE INITIAL COLUMN NORMS AND INITIALIZE SEVERAL ARRAYS.
!
      DO 10 J = 1, N
         ACNORM(J) = ENORM(M,A(1,J))
         RDIAG(J) = ACNORM(J)
         WA(J) = RDIAG(J)
         IF (PIVOT) IPVT(J) = J
   10    CONTINUE
!
!     REDUCE A TO R WITH HOUSEHOLDER TRANSFORMATIONS.
!
      MINMN = MIN0(M,N)
      DO 110 J = 1, MINMN
         IF (.NOT.PIVOT) GO TO 40
!
!        BRING THE COLUMN OF LARGEST NORM INTO THE PIVOT POSITION.
!
         KMAX = J
         DO 20 K = J, N
            IF (RDIAG(K) .GT. RDIAG(KMAX)) KMAX = K
   20       CONTINUE
         IF (KMAX .EQ. J) GO TO 40
         DO 30 I = 1, M
            TEMP = A(I,J)
            A(I,J) = A(I,KMAX)
            A(I,KMAX) = TEMP
   30       CONTINUE
         RDIAG(KMAX) = RDIAG(J)
         WA(KMAX) = WA(J)
         K = IPVT(J)
         IPVT(J) = IPVT(KMAX)
         IPVT(KMAX) = K
   40    CONTINUE
!
!        COMPUTE THE HOUSEHOLDER TRANSFORMATION TO REDUCE THE
!        J-TH COLUMN OF A TO A MULTIPLE OF THE J-TH UNIT VECTOR.
!
         AJNORM = ENORM(M-J+1,A(J,J))
         IF (AJNORM .EQ. ZERO) GO TO 100
         IF (A(J,J) .LT. ZERO) AJNORM = -AJNORM
         DO 50 I = J, M
            A(I,J) = A(I,J)/AJNORM
   50       CONTINUE
         A(J,J) = A(J,J) + ONE
!
!        APPLY THE TRANSFORMATION TO THE REMAINING COLUMNS
!        AND UPDATE THE NORMS.
!
         JP1 = J + 1
         IF (N .LT. JP1) GO TO 100
         DO 90 K = JP1, N
            SUM = ZERO
            DO 60 I = J, M
               SUM = SUM + A(I,J)*A(I,K)
   60          CONTINUE
            TEMP = SUM/A(J,J)
            DO 70 I = J, M
               A(I,K) = A(I,K) - TEMP*A(I,J)
   70          CONTINUE
            IF (.NOT.PIVOT .OR. RDIAG(K) .EQ. ZERO) GO TO 80
            TEMP = A(J,K)/RDIAG(K)
            RDIAG(K) = RDIAG(K)*DSQRT(DMAX1(ZERO,ONE-TEMP**2))
            IF (P05*(RDIAG(K)/WA(K))**2 .GT. EPSMCH) GO TO 80
            RDIAG(K) = ENORM(M-J,A(JP1,K))
            WA(K) = RDIAG(K)
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE
         RDIAG(J) = -AJNORM
  110    CONTINUE
      RETURN
!
!     LAST CARD OF SUBROUTINE QRFAC.
!
      END SUBROUTINE
!     **************************************************************
      SUBROUTINE QRSOLV(N,R,LDR,IPVT,DIAG,QTB,X,SDIAG,WA)
      INTEGER N,LDR
      INTEGER IPVT(N)
      real(dp) R(LDR,N),DIAG(N),QTB(N),X(N),SDIAG(N),WA(N)
!     ***************************************************************
!
!     SUBROUTINE QRSOLV
!
!     GIVEN AN M BY N MATRIX A, AN N BY N DIAGONAL MATRIX D,
!     AND AN M-VECTOR B, THE PROBLEM IS TO DETERMINE AN X WHICH
!     SOLVES THE SYSTEM
!
!           A*X = B ,     D*X = 0 ,
!
!     IN THE LEAST SQUARES SENSE.
!
!     THIS SUBROUTINE COMPLETES THE SOLUTION OF THE PROBLEM
!     IF IT IS PROVIDED WITH THE NECESSARY INFORMATION FROM THE
!     QR FACTORIZATION, WITH COLUMN PIVOTING, OF A. THAT IS, IF
!     A*P = Q*R, WHERE P IS A PERMUTATION MATRIX, Q HAS ORTHOGONAL
!     COLUMNS, AND R IS AN UPPER TRIANGULAR MATRIX WITH DIAGONAL
!     ELEMENTS OF NONINCREASING MAGNITUDE, THEN QRSOLV EXPECTS
!     THE FULL UPPER TRIANGLE OF R, THE PERMUTATION MATRIX P,
!     AND THE FIRST N COMPONENTS OF (Q TRANSPOSE)*B. THE SYSTEM
!     A*X = B, D*X = 0, IS THEN EQUIVALENT TO
!
!                  T       T
!           R*Z = Q *B ,  P *D*P*Z = 0 ,
!
!     WHERE X = P*Z. IF THIS SYSTEM DOES NOT HAVE FULL RANK,
!     THEN A LEAST SQUARES SOLUTION IS OBTAINED. ON OUTPUT QRSOLV
!     ALSO PROVIDES AN UPPER TRIANGULAR MATRIX S SUCH THAT
!
!            T   T               T
!           P *(A *A + D*D)*P = S *S .
!
!     S IS COMPUTED WITHIN QRSOLV AND MAY BE OF SEPARATE INTEREST.
!
!  THE SUBROUTINE STATEMENT IS
!
!       SUBROUTINE QRSOLV(N,R,LDR,IPVT,DIAG,QTB,X,SDIAG,WA)
!
!     WHERE
!
!       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE ORDER OF R.
!
!       R IS AN N BY N ARRAY. ON INPUT THE FULL UPPER TRIANGLE
!         MUST CONTAIN THE FULL UPPER TRIANGLE OF THE MATRIX R.
!         ON OUTPUT THE FULL UPPER TRIANGLE IS UNALTERED, AND THE
!         STRICT LOWER TRIANGLE CONTAINS THE STRICT UPPER TRIANGLE
!         (TRANSPOSED) OF THE UPPER TRIANGULAR MATRIX S.
!
!       LDR IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN N
!         WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY R.
!
!       IPVT IS AN INTEGER INPUT ARRAY OF LENGTH N WHICH DEFINES THE
!         PERMUTATION MATRIX P SUCH THAT A*P = Q*R. COLUMN J OF P
!         IS COLUMN IPVT(J) OF THE IDENTITY MATRIX.
!
!       DIAG IS AN INPUT ARRAY OF LENGTH N WHICH MUST CONTAIN THE
!         DIAGONAL ELEMENTS OF THE MATRIX D.
!
!       QTB IS AN INPUT ARRAY OF LENGTH N WHICH MUST CONTAIN THE FIRST
!         N ELEMENTS OF THE VECTOR (Q TRANSPOSE)*B.
!
!       X IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE LEAST
!         SQUARES SOLUTION OF THE SYSTEM A*X = B, D*X = 0.
!
!       SDIAG IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE
!         DIAGONAL ELEMENTS OF THE UPPER TRIANGULAR MATRIX S.
!
!       WA IS A WORK ARRAY OF LENGTH N.
!
!     SUBPROGRAMS CALLED
!
!       FORTRAN-SUPPLIED ... DABS,DSQRT
!
!     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
!
!     **********
      INTEGER I,J,JP1,K,KP1,L,NSING
      real(dp) COS,COTAN,P5,P25,QTBPJ,SIN,SUM,TAN,TEMP,ZERO
      DATA P5,P25,ZERO /5.0D-1,2.5D-1,0.0D0/
!
!     COPY R AND (Q TRANSPOSE)*B TO PRESERVE INPUT AND INITIALIZE S.
!     IN PARTICULAR, SAVE THE DIAGONAL ELEMENTS OF R IN X.
!
      DO 20 J = 1, N
         DO 10 I = J, N
            R(I,J) = R(J,I)
   10       CONTINUE
         X(J) = R(J,J)
         WA(J) = QTB(J)
   20    CONTINUE
!
!     ELIMINATE THE DIAGONAL MATRIX D USING A GIVENS ROTATION.
!
      DO 100 J = 1, N
!
!        PREPARE THE ROW OF D TO BE ELIMINATED, LOCATING THE
!        DIAGONAL ELEMENT USING P FROM THE QR FACTORIZATION.
!
         L = IPVT(J)
         IF (DIAG(L) .EQ. ZERO) GO TO 90
         DO 30 K = J, N
            SDIAG(K) = ZERO
   30       CONTINUE
         SDIAG(J) = DIAG(L)
!
!        THE TRANSFORMATIONS TO ELIMINATE THE ROW OF D
!        MODIFY ONLY A SINGLE ELEMENT OF (Q TRANSPOSE)*B
!        BEYOND THE FIRST N, WHICH IS INITIALLY ZERO.
!
         QTBPJ = ZERO
         DO 80 K = J, N
!
!           DETERMINE A GIVENS ROTATION WHICH ELIMINATES THE
!           APPROPRIATE ELEMENT IN THE CURRENT ROW OF D.
!
            IF (SDIAG(K) .EQ. ZERO) GO TO 70
            IF (DABS(R(K,K)) .GE. DABS(SDIAG(K))) GO TO 40
               COTAN = R(K,K)/SDIAG(K)
               SIN = P5/DSQRT(P25+P25*COTAN**2)
               COS = SIN*COTAN
               GO TO 50
   40       CONTINUE
               TAN = SDIAG(K)/R(K,K)
               COS = P5/DSQRT(P25+P25*TAN**2)
               SIN = COS*TAN
   50       CONTINUE
!
!           COMPUTE THE MODIFIED DIAGONAL ELEMENT OF R AND
!           THE MODIFIED ELEMENT OF ((Q TRANSPOSE)*B,0).
!
            R(K,K) = COS*R(K,K) + SIN*SDIAG(K)
            TEMP = COS*WA(K) + SIN*QTBPJ
            QTBPJ = -SIN*WA(K) + COS*QTBPJ
            WA(K) = TEMP
!
!           ACCUMULATE THE TRANFORMATION IN THE ROW OF S.
!
            KP1 = K + 1
            IF (N .LT. KP1) GO TO 70
            DO 60 I = KP1, N
               TEMP = COS*R(I,K) + SIN*SDIAG(I)
               SDIAG(I) = -SIN*R(I,K) + COS*SDIAG(I)
               R(I,K) = TEMP
   60          CONTINUE
   70       CONTINUE
   80       CONTINUE
   90    CONTINUE
!
!        STORE THE DIAGONAL ELEMENT OF S AND RESTORE
!        THE CORRESPONDING DIAGONAL ELEMENT OF R.
!
         SDIAG(J) = R(J,J)
         R(J,J) = X(J)
  100    CONTINUE
!
!     SOLVE THE TRIANGULAR SYSTEM FOR Z. IF THE SYSTEM IS
!     SINGULAR, THEN OBTAIN A LEAST SQUARES SOLUTION.
!
      NSING = N
      DO 110 J = 1, N
         IF (SDIAG(J) .EQ. ZERO .AND. NSING .EQ. N) NSING = J - 1
         IF (NSING .LT. N) WA(J) = ZERO
  110    CONTINUE
      IF (NSING .LT. 1) GO TO 150
      DO 140 K = 1, NSING
         J = NSING - K + 1
         SUM = ZERO
         JP1 = J + 1
         IF (NSING .LT. JP1) GO TO 130
         DO 120 I = JP1, NSING
            SUM = SUM + R(I,J)*WA(I)
  120       CONTINUE
  130    CONTINUE
         WA(J) = (WA(J) - SUM)/SDIAG(J)
  140    CONTINUE
  150 CONTINUE
!
!     PERMUTE THE COMPONENTS OF Z BACK TO COMPONENTS OF X.
!
      DO 160 J = 1, N
         L = IPVT(J)
         X(L) = WA(J)
  160    CONTINUE
      RETURN
!
!     LAST CARD OF SUBROUTINE QRSOLV.
!
      END SUBROUTINE

endmodule lmdif_mod
