!+
! This file originally was pgsubs.f from the demo section of the PGPLOT library.
!-

! This #if def wraps the entire source

#if defined (CESR_PGPLOT)

SUBROUTINE PALETT(TYPE, CONTRA, BRIGHT)
!-----------------------------------------------------------------------
! Set a "palette" of colors in the range of color indices used by
! PGIMAG.
!-----------------------------------------------------------------------
INTEGER TYPE
REAL CONTRA, BRIGHT
!
REAL GL(2), GR(2), GG(2), GB(2)
REAL RL(9), RR(9), RG(9), RB(9)
REAL HL(5), HR(5), HG(5), HB(5)
REAL WL(10), WR(10), WG(10), WB(10)
REAL AL(20), AR(20), AG(20), AB(20)
!
DATA GL /0.0, 1.0/
DATA GR /0.0, 1.0/
DATA GG /0.0, 1.0/
DATA GB /0.0, 1.0/
!
DATA RL /-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7/
DATA RR / 0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0/
DATA RG / 0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0/
DATA RB / 0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0/
!
DATA HL /0.0, 0.2, 0.4, 0.6, 1.0/
DATA HR /0.0, 0.5, 1.0, 1.0, 1.0/
DATA HG /0.0, 0.0, 0.5, 1.0, 1.0/
DATA HB /0.0, 0.0, 0.0, 0.3, 1.0/
!
DATA WL /0.0, 0.5, 0.5, 0.7, 0.7, 0.85, 0.85, 0.95, 0.95, 1.0/
DATA WR /0.0, 1.0, 0.0, 0.0, 0.3,  0.8,  0.3,  1.0,  1.0, 1.0/
DATA WG /0.0, 0.5, 0.4, 1.0, 0.0,  0.0,  0.2,  0.7,  1.0, 1.0/
DATA WB /0.0, 0.0, 0.0, 0.0, 0.4,  1.0,  0.0,  0.0, 0.95, 1.0/
!
DATA AL /0.0, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5, &
                 0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0/
DATA AR /0.0, 0.0, 0.3, 0.3, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0, &
                 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/
DATA AG /0.0, 0.0, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8, &
                 0.6, 0.6, 1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.0, 0.0/
DATA AB /0.0, 0.0, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.9, 0.9, &
                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/
!
IF (TYPE == 1) THEN
!        -- gray scale
   CALL PGCTAB(GL, GR, GG, GB, 2, CONTRA, BRIGHT)
ELSE IF (TYPE == 2) THEN
!        -- rainbow
   CALL PGCTAB(RL, RR, RG, RB, 9, CONTRA, BRIGHT)
ELSE IF (TYPE == 3) THEN
!        -- heat
   CALL PGCTAB(HL, HR, HG, HB, 5, CONTRA, BRIGHT)
ELSE IF (TYPE == 4) THEN
!        -- weird IRAF
   CALL PGCTAB(WL, WR, WG, WB, 10, CONTRA, BRIGHT)
ELSE IF (TYPE == 5) THEN
!        -- AIPS
   CALL PGCTAB(AL, AR, AG, AB, 20, CONTRA, BRIGHT)
END IF
END

SUBROUTINE SETVP
!-----------------------------------------------------------------------
! Set the viewport, allowing margins around the edge for annotation.
! (This is similar in effect to PGVSTD, but has different margins.)
! The routine determines the view-surface size and allocates margins
! as fractions of the minimum of width and height.
!-----------------------------------------------------------------------
REAL D, VPX1, VPX2, VPY1, VPY2
!
CALL PGSVP(0.0, 1.0, 0.0, 1.0)
CALL PGQVP(1, VPX1, VPX2, VPY1, VPY2)
D = MIN(VPX2-VPX1, VPY2-VPY1)/40.0
VPX1 = VPX1 + 5.0*D
VPX2 = VPX2 - 2.0*D
VPY1 = VPY1 + 8.0*D
VPY2 = VPY2 - 2.0*D
CALL PGVSIZ(VPX1, VPX2, VPY1, VPY2)
END

SUBROUTINE FIDDLE
!
INTEGER P, IER, PGCURS
REAL CONTRA, BRIGHT, X, Y, SIGN
REAL X1, Y1, X2, Y2, B1, B2, C1, C2
CHARACTER CH
!
WRITE (*,*) 'Use cursor to adjust color table:'
WRITE (*,*) ' Keys 1,2,3,4,5 select different palettes'
WRITE (*,*) ' Key P cycles through available palettes'
WRITE (*,*) ' Key F adjusts contrast and brightness, with'
WRITE (*,*) '  cursor x position setting brightness [0.0 - 1.0]'
WRITE (*,*) '   and y position setting contrast [0.0 - 10.0]'
WRITE (*,*) '  (Hold down F key while moving cursor to change'
WRITE (*,*) '  contrast and brightness continuously)'
WRITE (*,*) ' Key C resets contrast=1.0, brightness=0.5'
WRITE (*,*) ' Key - reverses color palette'
WRITE (*,*) ' Key X or right mouse button exits program'
!
P = 2
CONTRA = 1.0
BRIGHT = 0.5
X = 0.5
Y = 1.0
SIGN = +1.0
!
CALL PGQWIN(X1, X2, Y1, Y2)
B1 = 0.0
B2 = 1.0
C1 = 0.0
C2 = 10.0
CALL PGSWIN(B1, B2, C1, C2)
 10   IER = PGCURS(X, Y, CH)
IF (CH == CHAR(0) .OR. CH == 'x' .OR. CH == 'X') THEN
   CALL PGSWIN(X1, X2, Y1, Y2)
   RETURN
ELSE IF (CH == 'F' .OR. CH == 'f') THEN
   BRIGHT = MAX(B1, MIN(B2,X))
   CONTRA = MAX(C1, MIN(C2,Y))
ELSE IF (CH == 'C' .OR. CH == 'c') THEN
   CONTRA = 1.0
   Y = 1.0
   BRIGHT = 0.5
   X = 0.5
ELSE IF (CH == '-') THEN
   SIGN = -SIGN
ELSE IF (CH == '1') THEN
   P = 1
ELSE IF (CH == '2') THEN
   P = 2
ELSE IF (CH == '3') THEN
   P = 3
ELSE IF (CH == '4') THEN
   P = 4
ELSE IF (CH == '5') THEN
   P = 5
ELSE IF (CH == 'P' .OR. CH == 'p') THEN
   P = 1 + MOD(P,5)
END IF
CALL PALETT(P, SIGN*CONTRA, BRIGHT)
GOTO 10
END

SUBROUTINE FUNC(F, M, N, FMIN, FMAX)
INTEGER M,N
REAL F(M,N), FMIN, FMAX
!
INTEGER I, J
REAL R
!
FMIN = 1E30
FMAX = -1E30
DO 20 I=1,M
   DO 10 J=1,N
      R = SQRT(REAL(I)**2 + REAL(J)**2)
      F(I,J) = COS(0.6*SQRT(I*80./M)-16.0*J/(3.*N))* &
                   COS(16.0*I/(3.*M))+(I/REAL(M)-J/REAL(N)) + &
                   0.05*SIN(R)
      FMIN = MIN(F(I,J),FMIN)
      FMAX = MAX(F(I,J),FMAX)
 10      CONTINUE
 20   CONTINUE
END

SUBROUTINE OUTLIN(I1,I2,J1,J2,TR)
INTEGER I1,I2,J1,J2
REAL TR(6)
!-----------------------------------------------------------------------
! Draw the enclosing rectangle of the subarray to be contoured,
! applying the transformation TR.
!
! For a contour map, the corners are (I1,J1) and (I2,J2); for
! a gray-scale map, they are (I1-0.5,J1-0.5), (I2+0.5, J2+0.5).
!-----------------------------------------------------------------------
INTEGER K
REAL XW(5), YW(5), T
!
XW(1) = I1
YW(1) = J1
XW(2) = I1
YW(2) = J2
XW(3) = I2
YW(3) = J2
XW(4) = I2
YW(4) = J1
XW(5) = I1
YW(5) = J1
DO 10 K=1,5
    T = XW(K)
    XW(K) = TR(1) + TR(2)*T + TR(3)*YW(K)
    YW(K) = TR(4) + TR(5)*T + TR(6)*YW(K)
   10 CONTINUE
CALL PGLINE(5,XW,YW)
END

#endif



