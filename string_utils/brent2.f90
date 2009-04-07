!+
! Function BRENT2
!
! This is the subroutine BRENT from Numerical Recipes with the addition
! of a test to exit if the function F is not varying significantly.
!-

#include "CESR_platform.inc"

FUNCTION BRENT2 (AX,BX,CX,F,TOL,XMIN)


  use precision_def

  implicit none

  real(rp) brent2
  integer, PARAMETER :: ITMAX=100
  real(rp), parameter :: CGOLD=.3819660, ZEPS=1.0E-10
  real(rp) ax, bx, cx, tol, xmin
  real(rp) a, b, d, e, p, q, r, u, v, w, x, fu, fv, fw, fx
  real(rp) xm, tol1, tol2, etemp
  integer iter

  interface
     function f(x)
       use precision_def
       real(rp) f, x
     end function f
  end interface

  A=MIN(AX,CX)
  B=MAX(AX,CX)
  V=BX
  W=V
  X=V
  E=0.
  FX=F(X)
  FV=FX
  FW=FX
  DO 11 ITER=1,ITMAX
    XM=0.5*(A+B)
    TOL1=TOL*ABS(X)+ZEPS
    TOL2=2.*TOL1
    IF(ABS(X-XM)<=(TOL2-.5*(B-A))) GOTO 3
    IF(ABS(E)>TOL1) THEN
      R=(X-W)*(FX-FV)
      Q=(X-V)*(FX-FW)
      P=(X-V)*Q-(X-W)*R
      Q=2.*(Q-R)
      IF(Q>0.) P=-P
      Q=ABS(Q)
      ETEMP=E
      E=D
      IF(ABS(P)>=ABS(.5*Q*ETEMP).OR.P<=Q*(A-X).OR. &
          P>=Q*(B-X)) GOTO 1
      D=P/Q
      U=X+D
      IF(U-A<TOL2 .OR. B-U<TOL2) D=SIGN(TOL1,XM-X)
      GOTO 2
    ENDIF
1       IF(X>=XM) THEN
      E=A-X
    ELSE
      E=B-X
    ENDIF
    D=CGOLD*E
2       IF(ABS(D)>=TOL1) THEN
      U=X+D
    ELSE
      U=X+SIGN(TOL1,D)
    ENDIF
    FU=F(U)
    IF(FU<=FX) THEN
      IF(U>=X) THEN
        A=X
      ELSE
        B=X
      ENDIF
      V=W
      FV=FW
      W=X
      FW=FX
      X=U
      FX=FU
    ELSE
      IF(U<X) THEN
        A=U
      ELSE
        B=U
      ENDIF
      IF(FU<=FW .OR. W==X) THEN
        V=W
        FV=FW
        W=U
        FW=FU
      ELSE IF(FU<=FV .OR. V==X .OR. V==W) THEN
        V=U
        FV=FU
      ENDIF
    ENDIF

! Test to see if FX is not significantly different

    if ((iter > 3) .and. (fx >= max(fu, fv, fw) * 0.999999)) goto 3

11    CONTINUE
  print *,  'BRENT2: WARNING, MAXIMUM ITERATIONS EXCEEDED'
3     XMIN=X
  BRENT2=FX
  RETURN
  END
