MODULE fast_fourier_am

! Robert D. Ryne 06 February 2018
! Perform a 1d complex-to-complex fft.
! Minor edits starting with Alan Miller's code from https://jblevins.org/mirror/amiller/fft.f90
! I have removed the routines fft,mfft and use only fft1,mfft1.
! I also added top level routine ccfftam.
!
! Code converted using TO_F90 by Alan Miller
! Date: 2003-05-19  Time: 21:30:06
!
!ryne---------
use, intrinsic :: iso_fortran_env

implicit none

! Fortran 2008
integer, parameter, private :: dp = REAL64
!ryne---------



CONTAINS

subroutine ccfftam(c,n,isn,ierr)
implicit none
complex(dp), dimension(*), intent(inout) :: c !note that I got rank mismatch when * replaced with :
integer :: n,isn,ierr
real(dp), dimension(n) :: a,b
a(1:n)=real(c(1:n))
b(1:n)=aimag(c(1:n))
call fft1(a,b,n,isn,ierr)
c(1:n)=cmplx(a(1:n),b(1:n), dp)
return
end subroutine ccfftam

SUBROUTINE fft1 (a, b, n, isn, ierr)

REAL(dp), INTENT(IN OUT)  :: a(:)
REAL(dp), INTENT(IN OUT)  :: b(:)
INTEGER, INTENT(IN)       :: n
INTEGER, INTENT(IN)       :: isn
INTEGER, INTENT(OUT)      :: ierr

!     ------------
IF (ABS(isn) /= 1) GO TO 10
CALL sfft (a, b, n, n, n, isn, ierr)
RETURN
10 ierr = 4
RETURN
END SUBROUTINE fft1

SUBROUTINE mfft1 (a, b, n, ndim, isn, ierr)

REAL(dp), INTENT(IN OUT)  :: a(:)
REAL(dp), INTENT(IN OUT)  :: b(:)
INTEGER, INTENT(IN)       :: n(:)
INTEGER, INTENT(IN)       :: ndim
INTEGER, INTENT(IN)       :: isn
INTEGER, INTENT(OUT)      :: ierr

! Local variables
INTEGER  :: i, nspan, ntot
!     ------------
IF (ABS(isn) /= 1) GO TO 40
IF (ndim <= 0) GO TO 50
ntot = 1
DO  i = 1,ndim
  ntot = n(i)*ntot
END DO
IF (ntot < 1) GO TO 30

nspan = 1
DO  i = 1,ndim
  nspan = n(i)*nspan
  CALL sfft (a, b, ntot, n(i), nspan, isn, ierr)
  IF (ierr /= 0) RETURN
END DO
RETURN

30 ierr = 1
RETURN
40 ierr = 4
RETURN
50 ierr = 5
RETURN
END SUBROUTINE mfft1


SUBROUTINE sfft(a, b, ntot, n, nspan, isn, ierr)
!  MULTIVARIATE COMPLEX FOURIER TRANSFORM, COMPUTED IN PLACE USING
!    MIXED-RADIX FAST FOURIER TRANSFORM ALGORITHM.
!  BY R. C. SINGLETON, STANFORD RESEARCH INSTITUTE, OCT. 1968
!    MODIFIED BY A. H. MORRIS, NSWC/DL, DAHLGREN VA
!  ARRAYS A AND B ORIGINALLY HOLD THE REAL AND IMAGINARY COMPONENTS OF
!    THE DATA, AND RETURN THE REAL AND IMAGINARY COMPONENTS OF THE
!    RESULTING FOURIER COEFFICIENTS.
!  MULTIVARIATE DATA IS INDEXED ACCORDING TO THE FORTRAN ARRAY ELEMENT
!    SUCCESSOR FUNCTION, WITHOUT LIMIT ON THE NUMBER OF IMPLIED MULTIPLE
!    SUBSCRIPTS.
!    THE SUBROUTINE IS CALLED ONCE FOR EACH VARIATE.
!    THE CALLS FOR A MULTIVARIATE TRANSFORM MAY BE IN ANY ORDER.
!  NTOT IS THE TOTAL NUMBER OF COMPLEX DATA VALUES.
!  N IS THE DIMENSION OF THE CURRENT VARIABLE.
!  NSPAN/N IS THE SPACING OF CONSECUTIVE DATA VALUES
!    WHILE INDEXING THE CURRENT VARIABLE.
!  THE SIGN OF ISN DETERMINES THE SIGN OF THE COMPLEX EXPONENTIAL,
!    AND THE MAGNITUDE OF ISN IS NORMALLY ONE.
!  A TRI-VARIATE TRANSFORM WITH A(N1,N2,N3), B(N1,N2,N3) IS COMPUTED BY
!      CALL SFFT(A, B, N1*N2*N3, N1, N1, 1, IERR)
!      CALL SFFT(A, B, N1*N2*N3, N2, N1*N2, 1, IERR)
!      CALL SFFT(A, B, N1*N2*N3, N3, N1*N2*N3, 1, IERR)
!  FOR A SINGLE-VARIATE TRANSFORM,
!    NTOT = N = NSPAN = (NUMBER OF COMPLEX DATA VALUES), E.G.
!      CALL SFFT(A, B, N, N, N, 1, IERR)
!  THE DATA MAY ALTERNATIVELY BE STORED IN A SINGLE COMPLEX ARRAY A,
!    THEN THE MAGNITUDE OF ISN CHANGED TO TWO TO GIVE THE CORRECT INDEXING
!    INCREMENT AND A(2) USED TO PASS THE INITIAL ADDRESS FOR THE SEQUENCE
!    OF IMAGINARY VALUES, E.G.
!      CALL SFFT(A, A(2), NTOT, N, NSPAN, 2, IERR)
!  ARRAYS NFAC(MAXN), NP(MAXP), AT(MAXF), CK(MAXF), BT(MAXF), SK(MAXF)
!    ARE USED FOR TEMPORARY STORAGE.
!    MAXN MUST BE >= THE NUMBER OF FACTORS OF N
!    MAXF MUST BE >= THE MAXIMUM PRIME FACTOR OF N.
!    MAXP MUST BE > THE NUMBER OF PRIME FACTORS OF N.
!    IN ADDITION, MAXN IS ASSUMED TO BE ODD.
!    IF THE SQUARE-FREE PORTION K OF N HAS TWO OR MORE PRIME FACTORS,
!    THEN MAXP MUST BE >= K-1.
!  IERR IS A VARIABLE. IERR IS SET TO 0 IF NO INPUT ERRORS ARE
!    DETECTED. OTHERWISE, IERR IS ASSIGNED ONE OF THE VALUES
!      IERR=1    N IS LESS THAN 1
!      IERR=2    N HAS MORE THAN MAXN FACTORS
!      IERR=3    N HAS A PRIME FACTOR GREATER THAN MAXF OR THE SQUARE-FREE
!                PORTION OF N IS GREATER THAN MAXP+1

REAL(dp), INTENT(IN OUT)  :: a(*)
REAL(dp), INTENT(IN OUT)  :: b(*)
INTEGER, INTENT(IN)       :: ntot
INTEGER, INTENT(IN)       :: n
INTEGER, INTENT(IN)       :: nspan
INTEGER, INTENT(IN)       :: isn
INTEGER, INTENT(OUT)      :: ierr

!  ARRAY STORAGE IN NFAC FOR A MAXIMUM OF 15 FACTORS OF N.
!  IF N HAS MORE THAN ONE SQUARE-FREE FACTOR, THE PRODUCT OF THE
!    SQUARE-FREE FACTORS MUST BE <= 210
INTEGER   :: nfac(15), np(209)

!  ARRAY STORAGE FOR MAXIMUM PRIME FACTOR OF 23
REAL(dp)  :: at(23), ck(23), bt(23), sk(23)

INTEGER   :: i, inc, j, jc, jf, jj, k, k1, k2, k3, k4, kk, ks, kspan, kspnn,  &
             kt, l, m, max, maxf, maxn, maxp, nn, nt, num
REAL(dp)  :: aa, aj, ajm, ajp, ak, akm, akp, bb, bj, bjm, bjp, bk, bkm, bkp,  &
             c1, c2, c3, c72, cd, rad, radf, s1, s2, s3, s72, s120, sd, u, v

!  EQUIVALENCE (i,ii)
!  THE FOLLOWING CONSTANTS SHOULD AGREE WITH THE ARRAY DIMENSIONS.
maxn = 15
maxf = 23
maxp = 209
!  SET THE FOLLOWING CONSTANTS
!     RAD = 2.0*PI
!     S72 = SIN(RAD/5.0)
!     C72 = COS(RAD/5.0)
!     S120 = SQRT(0.75)
rad = 6.2831853071796_dp
s72 = .951056516295154_dp
c72 = .309016994374947_dp
s120 = .86602540378444_dp

ierr = 0
IF(n-1 < 0) THEN
  GO TO  1000
ELSE IF (n-1 == 0) THEN
  GO TO   960
END IF
inc = isn
IF(isn >= 0) GO TO 10
s72 = -s72
s120 = -s120
rad = -rad
inc = -inc
10 nt = inc*ntot
ks = inc*nspan
kspan = ks
nn = nt-inc
jc = ks/n
radf = rad*jc*0.5
i = 0
jf = 0

!  DETERMINE THE FACTORS OF N
m = 0
k = n
MAX = maxn/2
GO TO 20

15 IF(m == MAX) GO TO 1001
m = m+1
nfac(m) = 4
k = l
20 l = k/16
IF(k == l*16) GO TO 15
j = 3
jj = 9
GO TO 30

25 IF(m == MAX) GO TO 1001
m = m+1
nfac(m) = j
k = k/jj
30 IF(MOD(k,jj) == 0) GO TO 25
j = j+2
jj = j**2
IF(j <= maxf .AND. jj <= k) GO TO 30
IF(k > 4) GO TO 40
kt = m
nfac(m+1) = k
IF(k /= 1) m = m+1
GO TO 80

40 l = k/4
IF(k /= l*4) GO TO 50
IF(m == MAX) GO TO 1001
m = m+1
nfac(m) = 2
k = l
kt = m
IF(k == 1) GO TO 85
50 kt = m
IF(k-1 > maxp) GO TO 1002
num = maxn-kt-kt
j = 2
60 IF(MOD(k,j) /= 0) GO TO 70
m = m+1
nfac(m) = j
num = num-1
k = k/j
IF(k == 1) GO TO 80
IF(num <= 0) GO TO 1001
70 l = (j+1)/2
j = l+l+1
IF(j <= maxf) GO TO 60
GO TO 1002

80 IF(kt == 0) GO TO 100
85 j = kt
90 m = m+1
nfac(m) = nfac(j)
j = j-1
IF(j /= 0) GO TO 90

!  COMPUTE FOURIER TRANSFORM
100 sd = radf/kspan
cd = 2.0*SIN(sd)**2
sd = SIN(sd+sd)
kk = 1
i = i+1
IF(nfac(i) /= 2) GO TO 400

!  TRANSFORM FOR FACTOR OF 2 (INCLUDING ROTATION FACTOR)
kspan = kspan/2
k1 = kspan+2
210 k2 = kk+kspan
ak = a(k2)
bk = b(k2)
a(k2) = a(kk)-ak
b(k2) = b(kk)-bk
a(kk) = a(kk)+ak
b(kk) = b(kk)+bk
kk = k2+kspan
IF(kk <= nn) GO TO 210
kk = kk-nn
IF(kk <= jc) GO TO 210
IF(kk > kspan) GO TO 800
220 c1 = 1.0-cd
s1 = sd
230 k2 = kk+kspan
ak = a(kk)-a(k2)
bk = b(kk)-b(k2)
a(kk) = a(kk)+a(k2)
b(kk) = b(kk)+b(k2)
a(k2) = c1*ak-s1*bk
b(k2) = s1*ak+c1*bk
kk = k2+kspan
IF(kk < nt) GO TO 230
k2 = kk-nt
c1 = -c1
kk = k1-k2
IF(kk > k2) GO TO 230
u = sd*s1+cd*c1
v = sd*c1-cd*s1
ak = c1-u
s1 = s1+v

!  THE FOLLOWING THREE STATEMENTS COMPENSATE FOR TRUNCATION ERROR.
!    IF ROUNDED ARITHMETIC IS USED THEN ONE MAY SUBSTITUTE
!     C1 = AK
c1 = 1.5-0.5*(ak*ak+s1*s1)
s1 = c1*s1
c1 = c1*ak
kk = kk+jc
IF(kk < k2) GO TO 230
k1 = k1+inc+inc
kk = (k1-kspan)/2+jc
IF(kk <= jc+jc) GO TO 220
GO TO 100

!  TRANSFORM FOR FACTOR OF 3 (OPTIONAL CODE)
320 k1 = kk+kspan
k2 = k1+kspan
ak = a(kk)
bk = b(kk)
aj = a(k1)+a(k2)
bj = b(k1)+b(k2)
a(kk) = ak+aj
b(kk) = bk+bj
ak = -0.5*aj+ak
bk = -0.5*bj+bk
aj = (a(k1)-a(k2))*s120
bj = (b(k1)-b(k2))*s120
a(k1) = ak-bj
b(k1) = bk+aj
a(k2) = ak+bj
b(k2) = bk-aj
kk = k2+kspan
IF(kk < nn) GO TO 320
kk = kk-nn
IF(kk <= kspan) GO TO 320
GO TO 700

!  TRANSFORM FOR FACTOR OF 4
400 IF(nfac(i) /= 4) GO TO 600
kspnn = kspan
kspan = kspan/4
410 c1 = 1.0
s1 = 0.0
420 k1 = kk+kspan
k2 = k1+kspan
k3 = k2+kspan
akp = a(kk)+a(k2)
akm = a(kk)-a(k2)
ajp = a(k1)+a(k3)
ajm = a(k1)-a(k3)
a(kk) = akp+ajp
ajp = akp-ajp
bkp = b(kk)+b(k2)
bkm = b(kk)-b(k2)
bjp = b(k1)+b(k3)
bjm = b(k1)-b(k3)
b(kk) = bkp+bjp
bjp = bkp-bjp
IF(isn < 0) GO TO 450
akp = akm-bjm
akm = akm+bjm
bkp = bkm+ajm
bkm = bkm-ajm
IF(s1 == 0.0) GO TO 460
430 a(k1) = akp*c1-bkp*s1
b(k1) = akp*s1+bkp*c1
a(k2) = ajp*c2-bjp*s2
b(k2) = ajp*s2+bjp*c2
a(k3) = akm*c3-bkm*s3
b(k3) = akm*s3+bkm*c3
kk = k3+kspan
IF(kk <= nt) GO TO 420

440 u = sd*s1+cd*c1
v = sd*c1-cd*s1
c2 = c1-u
s1 = s1+v
!  THE FOLLOWING THREE STATEMENTS COMPENSATE FOR TRUNCATION ERROR.
!    IF ROUNDED ARITHMETIC IS USED THEN ONE MAY SUBSTITUTE
!     C1 = C2
c1 = 1.5-0.5*(c2*c2+s1*s1)
s1 = c1*s1
c1 = c1*c2
c2 = c1*c1-s1*s1
s2 = 2.0*c1*s1
c3 = c2*c1-s2*s1
s3 = c2*s1+s2*c1
kk = kk-nt+jc
IF(kk <= kspan) GO TO 420
kk = kk-kspan+inc
IF(kk <= jc) GO TO 410
IF(kspan == jc) GO TO 800
GO TO 100

450 akp = akm+bjm
akm = akm-bjm
bkp = bkm-ajm
bkm = bkm+ajm
IF(s1 /= 0.0) GO TO 430
460 a(k1) = akp
b(k1) = bkp
a(k2) = ajp
b(k2) = bjp
a(k3) = akm
b(k3) = bkm
kk = k3+kspan
IF(kk <= nt) GO TO 420
GO TO 440

!  TRANSFORM FOR FACTOR OF 5 (OPTIONAL CODE)
510 c2 = c72**2-s72**2
s2 = 2.0*c72*s72
520 k1 = kk+kspan
k2 = k1+kspan
k3 = k2+kspan
k4 = k3+kspan
akp = a(k1)+a(k4)
akm = a(k1)-a(k4)
bkp = b(k1)+b(k4)
bkm = b(k1)-b(k4)
ajp = a(k2)+a(k3)
ajm = a(k2)-a(k3)
bjp = b(k2)+b(k3)
bjm = b(k2)-b(k3)
aa = a(kk)
bb = b(kk)
a(kk) = aa+akp+ajp
b(kk) = bb+bkp+bjp
ak = akp*c72+ajp*c2+aa
bk = bkp*c72+bjp*c2+bb
aj = akm*s72+ajm*s2
bj = bkm*s72+bjm*s2
a(k1) = ak-bj
a(k4) = ak+bj
b(k1) = bk+aj
b(k4) = bk-aj
ak = akp*c2+ajp*c72+aa
bk = bkp*c2+bjp*c72+bb
aj = akm*s2-ajm*s72
bj = bkm*s2-bjm*s72
a(k2) = ak-bj
a(k3) = ak+bj
b(k2) = bk+aj
b(k3) = bk-aj
kk = k4+kspan
IF(kk < nn) GO TO 520
kk = kk-nn
IF(kk <= kspan) GO TO 520
GO TO 700

!  TRANSFORM FOR ODD FACTORS
600 k = nfac(i)
kspnn = kspan
kspan = kspan/k
IF(k == 3) GO TO 320
IF(k == 5) GO TO 510
IF(k == jf) GO TO 640
jf = k
s1 = rad/k
c1 = COS(s1)
s1 = SIN(s1)
ck(jf) = 1.0
sk(jf) = 0.0
j = 1
630 ck(j) = ck(k)*c1+sk(k)*s1
sk(j) = ck(k)*s1-sk(k)*c1
k = k-1
ck(k) = ck(j)
sk(k) = -sk(j)
j = j+1
IF(j < k) GO TO 630
640 k1 = kk
k2 = kk+kspnn
aa = a(kk)
bb = b(kk)
ak = aa
bk = bb
j = 1
k1 = k1+kspan
650 k2 = k2-kspan
j = j+1
at(j) = a(k1)+a(k2)
ak = at(j)+ak
bt(j) = b(k1)+b(k2)
bk = bt(j)+bk
j = j+1
at(j) = a(k1)-a(k2)
bt(j) = b(k1)-b(k2)
k1 = k1+kspan
IF(k1 < k2) GO TO 650
a(kk) = ak
b(kk) = bk
k1 = kk
k2 = kk+kspnn
j = 1
660 k1 = k1+kspan
k2 = k2-kspan
jj = j
ak = aa
bk = bb
aj = 0.0
bj = 0.0
k = 1
670 k = k+1
ak = at(k)*ck(jj)+ak
bk = bt(k)*ck(jj)+bk
k = k+1
aj = at(k)*sk(jj)+aj
bj = bt(k)*sk(jj)+bj
jj = jj+j
IF(jj > jf) jj = jj-jf
IF(k < jf) GO TO 670
k = jf-j
a(k1) = ak-bj
b(k1) = bk+aj
a(k2) = ak+bj
b(k2) = bk-aj
j = j+1
IF(j < k) GO TO 660
kk = kk+kspnn
IF(kk <= nn) GO TO 640
kk = kk-nn
IF(kk <= kspan) GO TO 640

!  MULTIPLY BY ROTATION FACTOR (EXCEPT FOR FACTORS OF 2 AND 4)
700 IF(i == m) GO TO 800
kk = jc+1
710 c2 = 1.0-cd
s1 = sd
720 c1 = c2
s2 = s1
kk = kk+kspan
730 ak = a(kk)
a(kk) = c2*ak-s2*b(kk)
b(kk) = s2*ak+c2*b(kk)
kk = kk+kspnn
IF(kk <= nt) GO TO 730
ak = s1*s2
s2 = s1*c2+c1*s2
c2 = c1*c2-ak
kk = kk-nt+kspan
IF(kk <= kspnn) GO TO 730
u = sd*s1+cd*c1
v = sd*c1-cd*s1
c2 = c1-u
s1 = s1+v
!  THE FOLLOWING THREE STATEMENTS COMPENSATE FOR TRUNCATION ERROR.
!    IF ROUNDED ARITHMETIC IS USED THEN THEY MAY BE DELETED.
c1 = 1.5-0.5*(c2*c2+s1*s1)
s1 = c1*s1
c2 = c1*c2
kk = kk-kspnn+jc
IF(kk <= kspan) GO TO 720
kk = kk-kspan+jc+inc
IF(kk <= jc+jc) GO TO 710
GO TO 100

!  PERMUTE THE RESULTS TO NORMAL ORDER---DONE IN TWO STAGES
!  PERMUTATION FOR SQUARE FACTORS OF N
800 np(1) = ks
IF(kt == 0) GO TO 890
k = kt+kt+1
IF(m < k) k = k-1
j = 1
np(k+1) = jc
810 np(j+1) = np(j)/nfac(j)
np(k) = np(k+1)*nfac(j)
j = j+1
k = k-1
IF(j < k) GO TO 810
k3 = np(k+1)
kspan = np(2)
kk = jc+1
k2 = kspan+1
j = 1
IF(n /= ntot) GO TO 850

!  PERMUTATION FOR SINGLE-VARIATE TRANSFORM (OPTIONAL CODE)
820 ak = a(kk)
a(kk) = a(k2)
a(k2) = ak
bk = b(kk)
b(kk) = b(k2)
b(k2) = bk
kk = kk+inc
k2 = kspan+k2
IF(k2 < ks) GO TO 820
830 k2 = k2-np(j)
j = j+1
k2 = np(j+1)+k2
IF(k2 > np(j)) GO TO 830
j = 1
840 IF(kk < k2) GO TO 820
kk = kk+inc
k2 = kspan+k2
IF(k2 < ks) GO TO 840
IF(kk < ks) GO TO 830
jc = k3
GO TO 890

!  PERMUTATION FOR MULTIVARIATE TRANSFORM
850 k = kk+jc
860 ak = a(kk)
a(kk) = a(k2)
a(k2) = ak
bk = b(kk)
b(kk) = b(k2)
b(k2) = bk
kk = kk+inc
k2 = k2+inc
IF(kk < k) GO TO 860
kk = kk+ks-jc
k2 = k2+ks-jc
IF(kk < nt) GO TO 850
k2 = k2-nt+kspan
kk = kk-nt+jc
IF(k2 < ks) GO TO 850
870 k2 = k2-np(j)
j = j+1
k2 = np(j+1)+k2
IF(k2 > np(j)) GO TO 870
j = 1
880 IF(kk < k2) GO TO 850
kk = kk+jc
k2 = kspan+k2
IF(k2 < ks) GO TO 880
IF(kk < ks) GO TO 870
jc = k3
890 IF(2*kt+1 >= m) RETURN
kspnn = np(kt+1)

!  PERMUTATION FOR SQUARE-FREE FACTORS OF N
j = m-kt
nfac(j+1) = 1
900 nfac(j) = nfac(j)*nfac(j+1)
j = j-1
IF(j /= kt) GO TO 900
kt = kt+1
nn = nfac(kt)-1
jj = 0
j = 0
GO TO 906

902 jj = jj-k2
k2 = kk
k = k+1
kk = nfac(k)
904 jj = kk+jj
IF(jj >= k2) GO TO 902
np(j) = jj
906 k2 = nfac(kt)
k = kt+1
kk = nfac(k)
j = j+1
IF(j <= nn) GO TO 904

!  DETERMINE THE PERMUTATION CYCLES OF LENGTH GREATER THAN 1
j = 0
GO TO 914

910 k = kk
kk = np(k)
np(k) = -kk
IF(kk /= j) GO TO 910
k3 = kk
914 j = j+1
kk = np(j)
IF(kk < 0) GO TO 914
IF(kk /= j) GO TO 910
np(j) = -j
IF(j /= nn) GO TO 914
maxf = inc*maxf
!  REORDER A AND B, FOLLOWING THE PERMUTATION CYCLES
GO TO 950

924 j = j-1
IF(np(j) < 0) GO TO 924
jj = jc
926 kspan = jj
IF(jj > maxf) kspan = maxf
jj = jj-kspan
k = np(j)
kk = jc*k+i+jj
k1 = kk+kspan
k2 = 0
928 k2 = k2+1
at(k2) = a(k1)
bt(k2) = b(k1)
k1 = k1-inc
IF(k1 /= kk) GO TO 928
932 k1 = kk+kspan
k2 = k1-jc*(k+np(k))
k = -np(k)
936 a(k1) = a(k2)
b(k1) = b(k2)
k1 = k1-inc
k2 = k2-inc
IF(k1 /= kk) GO TO 936
kk = k2
IF(k /= j) GO TO 932
k1 = kk+kspan
k2 = 0
940 k2 = k2+1
a(k1) = at(k2)
b(k1) = bt(k2)
k1 = k1-inc
IF(k1 /= kk) GO TO 940
IF(jj /= 0) GO TO 926
IF(j /= 1) GO TO 924

950 j = k3+1
nt = nt-kspnn
i = nt-inc+1
IF(nt >= 0) GO TO 924
960 RETURN

!  ERROR FINISH - THERE IS AN INPUT ERROR
1000 ierr = 1
RETURN
1001 ierr = 2
RETURN
1002 ierr = 3
RETURN
END SUBROUTINE sfft

END MODULE fast_fourier_am
