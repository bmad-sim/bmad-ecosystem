        SUBROUTINE MINOP(N,X,F,G,STEP,ACC,MAXFUN,NITYP,INIT)
        implicit none
        integer:: n,maxfun,nityp,i,j,k,nfcall,ngcall,iupd,iptest
        real:: step,acc,epssq,qsq,dss,f,fa,c,d,ggdiag,hdiag
        real:: gggg,hggg,dg,dga,dggd,wbsq,b0,c1,c2,c3
        real:: a,b,e,cc,ca,cb,theta,dsq,flag,sum,zwc,zgv,vgv,vsq
        real:: zv,phi,usq,gsq,uv,sigma,bound,wasq,zgz,vwc
        logical:: init
        integer,parameter:: mq=100
        real:: X(*), G(*),XA(mq), GA(mq),H(mq,mq), GG(mq,mq)
        real:: V(mq),U(mq), WA(mq), WB(mq), WC(mq), WT(mq),WS(mq)
        real:: Y(mq), Z(mq), GZ(mq), GV(mq), PZ(mq)
!     ******************************************************************
!     * INITIATLIZE THE FOLLOWING PARAMETERS                           *
!     * NFCALL & NGCALL = NO. OF FUNCTION & GRADIENT CALLS RESP..      *
!     ******************************************************************
	QSQ=1.			!FOR COMPILER
	NFCALL = 1
	NGCALL = 1
	IUPD = 1
	EPSSQ = ACC * ACC
	IPTEST = 1
	DSS = STEP*STEP
	CALL CALCF (N,X,F)
	CALL CALCG (N,X,G)
	print *,' MINOP 1ST EVAL ',F
	IF(INIT.NE.-1) GO TO 8 !SAVE HESSIAN IF NOT FIRST RUN
	GSQ = 0
	DO 3 I = 1,N
   3    GSQ = GSQ + G(I) * G(I)
!
!     ******************************************************************
!     * PREPARE FOR THE FIRST ITERATION BY INITIALIZING MATRICES H&GG  *
!     ******************************************************************
!
	C = -SQRT( DSS / GSQ )
	GGDIAG = 0.01 * SQRT( GSQ) / STEP
	HDIAG = 1. / GGDIAG
	DO 6 I= 1,N
	DO 5 J = 1,N
	H(I,J) = 0.0D0
	GG(I,J) = 0.0D0
    5   CONTINUE
	H(I,I) = HDIAG
	GG(I,I) = GGDIAG
    6   Z(I) = C * G(I)
    8   GSQ = 0
	DO 9 I = 1, N
    9   GSQ = GSQ + G(I) * G(I)
	IF(GSQ.LT.EPSSQ) print *,QSQ,' MINOP RETURN ON GRADIENT OK '
	IF(GSQ.LT.EPSSQ) RETURN
!     ******************************************************************
!     * TEST WHETHER MAXFUN CALLS OF CALCFG HAVE BEEN MADE             *
!     ******************************************************************
!
18	IF(NFCALL.EQ.1) CALL OOO(0,NFCALL,F)
 	IF(NFCALL.GT.MAXFUN) print *,' MINOP , FUNCTIONS EXHAUSTED'
	IF(NFCALL.GT.MAXFUN) RETURN !FUNC CALLS EXHAUSTED
	BOUND = SQRT (DSS)    !WAS FOR PRINTOUT
!     ******************************************************************
!     * CALCUALTE NEWTON POINT                                         *
!     ******************************************************************
   26   WASQ = 0.0D0
	DO 28 I = 1,N
	WA(I) = 0.0D0
	DO 27 J= 1,N
   27   WA(I) = WA(I) - H(I,J) * G(J)
   28   WASQ = WASQ + WA(I) * WA(I)
	IF ( WASQ .LE. DSS) GO TO 39
	GGGG = 0
	HGGG = 0.0D0
	DO 30 I = 1, N
	WB(I) = 0.0D0
	DO 29 J = 1,N
   29   WB(I) = WB(I) + GG(I,J) * G(J)
	HGGG = HGGG - WA(I) * G(I)
   30   GGGG = GGGG + WB(I) * G(I)
	D = HGGG * GGGG
	D = GSQ * GSQ / D
	D = 0.2 + 0.8 * D
	IF ( D * D * WASQ .GT. DSS ) GO TO 32
	D = SQRT ( DSS / WASQ )
	DO 31 I = 1,N
   31   WA(I) = D * WA(I)
	GO TO 41
32	CONTINUE
!	IF ((GGGG *(GGGG * DSS)).GT. GSQ ** 3 ) GO TO 34
	IF (((GGGG/GSQ) *(GGGG/GSQ) * DSS).GT. GSQ  ) GO TO 34
!
!     ******************************************************************
!     * THE CAUCHY POINT IS OUTSIDE CIRCLE, SO PICK THE POINT ON THE   *
!     * BOUNDARY                                                       *
!     ******************************************************************
!	if((gsq*1.0e-25).lt.dss) then
	 C = -SQRT(DSS)/sqrt(GSQ)
!	else  			!avoid float underflow
!	 c=-1.0e-15
!	endif
	DO 33 I = 1,N
   33   WA(I) = C * G(I)
	GO TO 41
!
!     ******************************************************************
!     * THE CAUCHY POINT IS INSIDE CIRCLE, USE DOG LEG STEP.           *
!     ******************************************************************
   34	DO 35 I = 1,N
   35   WA(I) = D * WA(I)
	CA = 0.0
	CB = 0
	C = -GSQ / GGGG
!     ******************************************************************
!     * SET THE STEEPEST DESCENT CORRECTION IN WB AND THE DIFFERENCE   *
!     * BETWEEN WA AND WB IN WC.                                       *
!     ******************************************************************
	DO 36 I = 1, N
	WB(I) = C * G(I)
	WC(I)= WA(I) - WB(I)
	CA = CA + WB(I) * WC(I)
   36   CB = CB + WC(I) * WC(I)
!
!     * FIND CORRECTION VECTOR AT THE INTERSECTION OF WA-WB AND CIRCLE *
	C = DSS -C * C * GSQ
	THETA = SIGN (C / (ABS(CA) + SQRT( CA * CA + C * CB)), CA)
!     * TEST WHETHER TO USE THE GENERALIZED NEWTON STEP                *
	IF (THETA-1.) 37, 39, 39
   37   DO 38 I = 1, N
   38   WA(I) = WB(I) + THETA * WC(I)
	GO TO 41
   39   CONTINUE
!
   41   DO 42 I = 1,N
   42   XA(I) = X(I) + WA(I)
	NFCALL = NFCALL + 1
	CALL CALCF(N,XA,FA)
   43   DSQ = 0.0D0
	DO 44 I = 1,N
   44   DSQ = DSQ + WA(I) * WA(I)
!     ******************************************************************
!     * TEST IF FUNCTION VALUE IS DECREASED                            *
!     ******************************************************************
	IF(MOD(NFCALL,NITYP).EQ.0) CALL OOO(0,NFCALL,FA)
	IF (FA .LT. F) GO TO 46
	DSS = 0.25 * DSQ
	GO TO 18   !SBP MOD
!
!     * CALCULATE SOME NUMBERS FOR REVISING STEP BOUND                 *
   46   NGCALL = NGCALL + 1
	CALL CALCG(N,XA,GA)
	DG = 0.0
	DGA = 0.0
	DGGD = 0.0
	WBSQ = 0.0
	DO 47 I = 1, N
	WC(I) = 0.0D0
	DO 47 J = 1,N
   47   WC(I) = WC(I) + GG(I,J) * WA(J)
!
!     * SET  (Y-GS) IN WB                                              *
	DO 48 I = 1,N
	U(I) = WC(I)
	WB(I) = GA(I) -G(I) - WC(I)
	DG = DG + G(I) * WA(I)
	DGA = DGA + GA(I) * WA(I)
	DGGD = DGGD + WC(I) * WA(I)
   48   WBSQ = WBSQ + WB(I) *WB(I)
!     ******************************************************************
!     * TEST WHETHER TO DECREASE THE STEP-BOUND                        *
!     ******************************************************************
	IF (FA-F-0.1 * DG-0.05 * DGGD) 51, 51, 50
   50   DSS = 0.25 * DSQ
	IF(DSS.LT.1E-15) print *,' STEP-BOUND<E-15  MINOP EXITING'
	IF(DSS.LT.1E-15) RETURN
	GO TO 54    !NOT TOO SMALL
!     ******************************************************************
!     * TEST WHETHER TO INCREASE THE STEP-BOUND                        *
!     ******************************************************************
   51   DSS = DSQ
	IF (WBSQ - 0.25 * GSQ) 53, 53, 52
   52   IF (DG -DGA -DGA) 54, 53, 53
   53   DSS = 4.0 * DSQ
!     * SET THE DIFFERENCE BETWEN GRADIENTS IN WC                      *
   54   DO 55 I = 1,N
   55   WC(I) = GA(I) - G(I)
	F = FA
	B0 = 0.0
	DO 56 I = 1,N
	X(I) = XA(I)
	G(I) = GA(I)
   56   B0 = B0 + WC(I) * WA(I)
!     ******************************************************************
!     * TEST WHETHER WE NEED TO UPDATE H AND G                         *
!     ******************************************************************
	IF (B0 .GE. 1.0D-30) GO TO 58
!   WRITE(6,57)  NFCALL
   57   FORMAT('FUNCTION DECREASED ON',I4,'TH ITERATION, NO UPDATE IS MADE &
       ,SINCE Y*Z IS NEG. OR TOO SMALL')
	GO TO 8
   58   IF (IUPD .EQ. 1) GO TO 63
	IF (FLAG .EQ. 0.0) GO TO 60
	DO  59 I = 1, N
   59   Z(I) = WA(I)
	GO TO 63
   60   C2=0.0
	C3 = 0.0
	SUM = 0.0
	DO 61 I = 1, N
	C2 = C2 + WS(I) * WA(I)
	C3 = C3 + WT(I) * WA(I)
   61   SUM = SUM +Y(I) * WA(I)
	DO 62 I = 1, N
   62   Z(I) = C1 *(C2 * Z(I) + C3 * V(I)) - (SUM/ZWC) * Z(I)
   63   IUPD = IUPD+1
	DO 64 I = 1, N
	GV(I) = U(I) - WC(I)
	V(I) = WA(I)
	GZ(I) = 0.0
	DO 64 J = 1, N
	GZ(I) = GZ(I) + GG(I,J) * Z(J)
   64   V(I) =V(I) - H(I,J) * WC(J)
	ZGV = 0.0
	VGV = 0.0
	ZWC = 0.0
	ZGZ = 0.0
	VWC = 0.0
	DO 65 I = 1, N
	ZWC = ZWC + Z(I) * WC(I)
	ZGZ = ZGZ + Z(I) * GZ(I)
	VWC = VWC + V(I) * WC(I)
	ZGV = ZGV + Z(I) * GV(I)
	VGV = VGV + V(I) * GV(I)
   65   CONTINUE
	C1 = 1. / (ZGZ * VGV - ZGV * ZGV)
	C2 = ZWC * VGV - VWC * ZGV
	C3 = -ZWC * ZGV + VWC * ZGZ
	VSQ =0.0
	ZV  =0.0
	A = 0.0
	E = 0.0
	DO 66 I = 1,N
	Y(I) = C1 * (C2 * GZ(I) + C3 * GV(I))
	WS(I) = VGV * GZ(I) - ZGV * GV(I)
	WT(I) = -ZGV * GZ(I) + ZGZ * GV(I)
	ZV  =  ZV + Z(I) * V(I)
	VSQ = VSQ + V(I) * V(I)
	A = A + WS(I) * WA(I)
   66   E = E + WT(I) * WA(I)
	B = 0.0
	DO 67 I = 1,N
	PZ(I) = C1 * ( A * Z(I) + E * V(I))
   67   B = B +PZ(I) * Y(I)
	IF (B .GE. 1.0D-30) GO TO 70
	B = B0
!      WRITE (6,68)
   68   FORMAT(15X,'IMAGE PRODUCT NEG.OR TOO SMALL,ORIGINAL STEP USED')
	DO 69 I = 1, N
	PZ(I) = WA(I)
   69   GA(I) = WC(I)
	GO TO 72
   70   DO 71 I = 1, N
   71   GA(I) = Y(I)
   72   CC = 0.0
	A  = 0.0
	DO 74 I = 1, N
	WB(I) = 0.0
	GZ(I) = 0.0
	DO  73 J = 1, N
	GZ(I) = GZ(I) + GG(I,J) *PZ(J)
   73   WB(I) = WB(I) + H(I,J) * GA(J)
	A = A + GA(I) * WB(I)
   74   CC = CC + PZ(I) * GZ(I)
!     ******************************************************************
!     *  OPTIMAL CONDITIONING                                          *
!     ******************************************************************
	IF ( B .LE. 2.0D0 * A * CC / (A + CC)) GO TO 75
	PHI = B / (B - A)
	GO TO 76
   75   PHI = B * (CC - B)/(A * CC - B * B)
   76   USQ = 0.0
	UV  = 0.0
!     ******************************************************************
!     * TEST WHETHER Z & V ARE LINEARLY INDEPENDENT                    *
!     ******************************************************************
	DO 77 I = 1, N
	U(I) = Z(I) - ( ZV / VSQ ) * V(I)
	USQ = USQ + U(I) * U(I)
   77   UV  = UV  + U(I) * V(I)
	FLAG = 0.0
	IF (1.0D6 * UV * UV  .GE. VSQ * USQ ) FLAG =1.0
!     ******************************************************************
!     *  IF FLAG = 1 THEN Z & V  ARE LINEARLY DEPENDENT                *
!     ******************************************************************
	SIGMA = 1.+ PHI * (A * CC/ (B * B) - 1.)
        IF ( SIGMA .LT. 1.0D-30) WRITE (6,78)
   78   FORMAT(20X,' G IS ALMOST SINGULAR')
	DO 79 I = 1, N
	WA(I) = (1./B) * PZ(I) - (1./A) * WB(I)
	WC(I) = (1./B) * GZ(I) -CC/(B * B) * GA(I)
   79   U(I) = GA(I)- GZ(I)
!     ******************************************************************
!     *  START TO UPDATE MATRICES H & GG                               *
!     ******************************************************************
	DO 80 I = 1, N
	DO 80 J = 1, N
	GG(I,J) = GG(I,J) + (1./B) * (U(I) * GA(J) + GA(I) * U(J)) &
               -(( B- CC ) / (B * B)) * GA(I) * GA(J) &
               -(PHI* A / SIGMA) * WC(I) * WC(J)
	H(I,J) = H(I,J) - (1./A) * WB(I) * WB(J) + (1./B) * PZ(I) * PZ(J) &
               + PHI * A * WA(I) * WA(J)
   80   CONTINUE
	GO TO 8
	END
!
