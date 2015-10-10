!+
! Function inverse (funct, y, x1, x2, tol) result (x)
!
! Function to invert a function Y = FUNCT(X) to return X given Y
! This function is a slight modification of the function ZBRENT given
!     in Numerical Recipes.
!
! Inputs:
!   funct(x) -- Function: Function to be inverted.
!   y        -- Real(rp): Value to be inverted.
!   x1, x2   -- Real(rp): X must lie between: X1 < X < X2.
!   tol      -- Real(rp): Absolute accuracy for X.
!
! Output:
!   x -- Real(rp): Inverted function.
!-

function inverse(funct, y, x1, x2, tol) result (x)

  use precision_def

  implicit none

  real(rp) x, y, tol
  real(rp) x1, x2, a, b, c, d, e, fa, fb, fc
  real(rp) p, q, r, s, xm, tol1
  integer iter

  integer, parameter :: itmax = 100
  real(rp), parameter :: eps = 3.d-8

  interface
     function funct(x)
       use precision_def
       implicit none
       real(rp) funct, x
     end function funct
  end interface

!-----------------------------------------------------

  A=X1
  B=X2
  FA=FUNCT(A) - y
  FB=FUNCT(B) - y

  IF(FB*FA>0.) then
    print *, 'Root must be bracketed for INVERSE.'
    if (global_com%exit_on_error) call err_exit
  endif

  C=B
  FC=FB
  DO 11 ITER=1,ITMAX
    IF(FB*FC>0.) THEN
      C=A
      FC=FA
      D=B-A
      E=D
    ENDIF
    IF(ABS(FC)<ABS(FB)) THEN
      A=B
      B=C
      C=A
      FA=FB
      FB=FC
      FC=FA
    ENDIF
    TOL1=2.*EPS*ABS(B)+0.5*TOL
    XM=.5*(C-B)
    IF(ABS(XM)<=TOL1 .OR. FB==0.)THEN
      X=B
      RETURN
    ENDIF
    IF(ABS(E)>=TOL1 .AND. ABS(FA)>ABS(FB)) THEN
      S=FB/FA
      IF(A==C) THEN
        P=2.*XM*S
        Q=1.-S
      ELSE
        Q=FA/FC
        R=FB/FC
        P=S*(2.*XM*Q*(Q-R)-(B-A)*(R-1.))
        Q=(Q-1.)*(R-1.)*(S-1.)
      ENDIF
      IF(P>0.) Q=-Q
      P=ABS(P)
      IF(2.*P < MIN(3.*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
        E=D
        D=P/Q
      ELSE
        D=XM
        E=D
      ENDIF
    ELSE
      D=XM
      E=D
    ENDIF
    A=B
    FA=FB
    IF(ABS(D) > TOL1) THEN
      B=B+D
    ELSE
      B=B+SIGN(TOL1,XM)
    ENDIF
    FB=FUNCT(B) - y
11    CONTINUE
  print *, 'INVERSE: exceeding maximum iterations.'
  X = B
END function
