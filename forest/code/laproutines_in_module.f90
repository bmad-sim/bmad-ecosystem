module la_constants
!  -- LAPACK auxiliary module --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

!  Standard constants for
   integer, parameter :: sp = kind(1.e0)

   real(sp), parameter :: szero = 0.0_sp
   real(sp), parameter :: shalf = 0.5_sp
   real(sp), parameter :: sone = 1.0_sp
   real(sp), parameter :: stwo = 2.0_sp
   real(sp), parameter :: sthree = 3.0_sp
   real(sp), parameter :: sfour = 4.0_sp
   real(sp), parameter :: seight = 8.0_sp
   real(sp), parameter :: sten = 10.0_sp
   complex(sp), parameter :: czero = ( 0.0_sp, 0.0_sp )
   complex(sp), parameter :: chalf = ( 0.5_sp, 0.0_sp )
   complex(sp), parameter :: cone = ( 1.0_sp, 0.0_sp )
   character*1, parameter :: sprefix = 'S'
   character*1, parameter :: cprefix = 'C'

!  Scaling constants
   real(sp), parameter :: sulp = epsilon(0._sp)
   real(sp), parameter :: seps = sulp * 0.5_sp

   real(sp), parameter :: ssafmin = real(radix(0._sp),sp)**max( &
      minexponent(0._sp)-1, &
      1-maxexponent(0._sp) &
   )

   real(sp), parameter :: ssafmax = sone / ssafmin
   real(sp), parameter :: ssmlnum = ssafmin / sulp
   real(sp), parameter :: sbignum = ssafmax * sulp
   real(sp), parameter :: srtmin = sqrt(ssmlnum)
   real(sp), parameter :: srtmax = sqrt(sbignum)

!  Blue's scaling constants

   real(sp), parameter :: stsml = real(radix(0._sp), sp)**ceiling( &
       (minexponent(0._sp) - 1) * 0.5_sp)

   real(sp), parameter :: stbig = real(radix(0._sp), sp)**floor( &
       (maxexponent(0._sp) - digits(0._sp) + 1) * 0.5_sp)

!  ssml >= 1/s, where s was defined in https://doi.org/10.1145/355769.355771
!  The correction was added in https://doi.org/10.1145/3061665 to scale denormalized numbers correctly

   real(sp), parameter :: sssml = real(radix(0._sp), sp)**( - floor( &
       (minexponent(0._sp) - digits(0._sp)) * 0.5_sp))

!  sbig = 1/S, where S was defined in https://doi.org/10.1145/355769.355771

   real(sp), parameter :: ssbig = real(radix(0._sp), sp)**( - ceiling( &
       (maxexponent(0._sp) + digits(0._sp) - 1) * 0.5_sp))

!  Standard constants for
   integer, parameter :: dp = kind(1.d0)

   real(dp), parameter :: dzero = 0.0_dp
   real(dp), parameter :: dhalf = 0.5_dp
   real(dp), parameter :: done = 1.0_dp
   real(dp), parameter :: dtwo = 2.0_dp
   real(dp), parameter :: dthree = 3.0_dp
   real(dp), parameter :: dfour = 4.0_dp
   real(dp), parameter :: deight = 8.0_dp
   real(dp), parameter :: dten = 10.0_dp
   complex(dp), parameter :: zzero = ( 0.0_dp, 0.0_dp )
   complex(dp), parameter :: zhalf = ( 0.5_dp, 0.0_dp )
   complex(dp), parameter :: zone = ( 1.0_dp, 0.0_dp )
   character*1, parameter :: dprefix = 'D'
   character*1, parameter :: zprefix = 'Z'

!  Scaling constants
   real(dp), parameter :: dulp = epsilon(0._dp)
   real(dp), parameter :: deps = dulp * 0.5_dp

   real(dp), parameter :: dsafmin = real(radix(0._dp),dp)**max( &
      minexponent(0._dp)-1, &
      1-maxexponent(0._dp) &
   )

   real(dp), parameter :: dsafmax = done / dsafmin
   real(dp), parameter :: dsmlnum = dsafmin / dulp
   real(dp), parameter :: dbignum = dsafmax * dulp
   real(dp), parameter :: drtmin = sqrt(dsmlnum)
   real(dp), parameter :: drtmax = sqrt(dbignum)

!  Blue's scaling constants

   real(dp), parameter :: dtsml = real(radix(0._dp), dp)**ceiling( &
       (minexponent(0._dp) - 1) * 0.5_dp)

   real(dp), parameter :: dtbig = real(radix(0._dp), dp)**floor( &
       (maxexponent(0._dp) - digits(0._dp) + 1) * 0.5_dp)

!  ssml >= 1/s, where s was defined in https://doi.org/10.1145/355769.355771
!  The correction was added in https://doi.org/10.1145/3061665 to scale denormalized numbers correctly

   real(dp), parameter :: dssml = real(radix(0._dp), dp)**( - floor( &
       (minexponent(0._dp) - digits(0._dp)) * 0.5_dp))

!  sbig = 1/S, where S was defined in https://doi.org/10.1145/355769.355771

   real(dp), parameter :: dsbig = real(radix(0._dp), dp)**( - ceiling( &
       (maxexponent(0._dp) + digits(0._dp) - 1) * 0.5_dp))

end module la_constants

!
module la_xisnan

   interface la_isnan

   module procedure sisnan
   module procedure disnan

   end interface

contains

   logical function sisnan( x )   ! etienne
   use la_constants, only: wp=>sp
!#ifdef USE_IEEE_INTRINSIC
!   use, intrinsic :: ieee_arithmetic
!#elif USE_ISNAN
!   intrinsic :: isnan
!#endif
   real(wp) :: x
!#ifdef USE_IEEE_INTRINSIC
!   sisnan = ieee_is_nan(x)
!#elif USE_ISNAN
!   sisnan = isnan(x)
!#else
   sisnan = slaisnan(x,x)

   contains

   logical function slaisnan( x, y )
   use la_constants, only: wp=>sp
   real(wp) :: x, y
   slaisnan = ( x.ne.y )

   end function slaisnan
!#endif

   end function sisnan

   logical function disnan( x )
   use la_constants, only: wp=>dp
!#ifdef USE_IEEE_INTRINSIC
!   use, intrinsic :: ieee_arithmetic
!#elif USE_ISNAN
!   intrinsic :: isnan
!#endif
   real(wp) :: x
!#ifdef USE_IEEE_INTRINSIC
!   disnan = ieee_is_nan(x)
!#elif USE_ISNAN
!   disnan = isnan(x)
!#else
   disnan = dlaisnan(x,x)

   contains

   logical function dlaisnan( x, y )
   use la_constants, only: wp=>dp
   real(wp) :: x, y
   dlaisnan = ( x.ne.y )

   end function dlaisnan
!#endif

   end function disnan

end module la_xisnan


!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION DNRM2(N,X,INCX)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION X(*)
!       ..
!
!
!  =============
!
!  Arguments:
!  ==========
!
!
!  Authors:
!  ========
!
!
!
!
!  ==================
!
!  =====================
!  =====================================================================

function dnrm2( n, x, incx ) 
   integer, parameter :: wp = kind(1.d0)
   real(wp) :: dnrm2
!
!  -- Reference BLAS level1 routine (version 3.9.1) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     March 2021
!
!  .. Constants ..
   real(wp), parameter :: zero = 0.0_wp
   real(wp), parameter :: one  = 1.0_wp
   real(wp), parameter :: maxn = huge(0.0_wp)
!  ..
!  .. Blue's scaling constants ..
   real(wp), parameter :: tsml = real(radix(0._wp), wp)**ceiling( &
       (minexponent(0._wp) - 1) * 0.5_wp)
   real(wp), parameter :: tbig = real(radix(0._wp), wp)**floor( &
       (maxexponent(0._wp) - digits(0._wp) + 1) * 0.5_wp)
   real(wp), parameter :: ssml = real(radix(0._wp), wp)**( - floor( &
       (minexponent(0._wp) - digits(0._wp)) * 0.5_wp))
   real(wp), parameter :: sbig = real(radix(0._wp), wp)**( - ceiling( &
       (maxexponent(0._wp) + digits(0._wp) - 1) * 0.5_wp))
!  ..
!  .. Scalar Arguments ..
   integer :: incx, n
!  ..
!  .. Array Arguments ..
   real(wp) :: x(*)
!  ..
!  .. Local Scalars ..
   integer :: i, ix
   logical :: notbig
   real(wp) :: abig, amed, asml, ax, scl, sumsq, ymax, ymin
!
!  Quick return if possible
!
   dnrm2 = zero
   if( n <= 0 ) return
!
   scl = one
   sumsq = zero
!
!  Compute the sum of squares in 3 accumulators:
!     abig -- sums of squares scaled down to avoid overflow
!     asml -- sums of squares scaled up to avoid underflow
!     amed -- sums of squares that do not require scaling
!  The thresholds and multipliers are
!     tbig -- values bigger than this are scaled down by sbig
!     tsml -- values smaller than this are scaled up by ssml
!
   notbig = .true.
   asml = zero
   amed = zero
   abig = zero
   ix = 1
   if( incx < 0 ) ix = 1 - (n-1)*incx
   do i = 1, n
      ax = abs(x(ix))
      if (ax > tbig) then
         abig = abig + (ax*sbig)**2
         notbig = .false.
      else if (ax < tsml) then
         if (notbig) asml = asml + (ax*ssml)**2
      else
         amed = amed + ax**2
      end if
      ix = ix + incx
   end do
!
!  Combine abig and amed or amed and asml if more than one
!  accumulator was used.
!
   if (abig > zero) then
!
!     Combine abig and amed if abig > 0.
!
      if ( (amed > zero) .or. (amed > maxn) .or. (amed /= amed) ) then
         abig = abig + (amed*sbig)*sbig
      end if
      scl = one / sbig
      sumsq = abig
   else if (asml > zero) then
!
!     Combine amed and asml if asml > 0.
!
      if ( (amed > zero) .or. (amed > maxn) .or. (amed /= amed) ) then
         amed = sqrt(amed)
         asml = sqrt(asml) / ssml
         if (asml > amed) then
            ymin = amed
            ymax = asml
         else
            ymin = asml
            ymax = amed
         end if
         scl = one
         sumsq = ymax**2*( one + (ymin/ymax)**2 )
      else
         scl = one / ssml
         sumsq = asml
      end if
   else
!
!     Otherwise all values are mid-range
!
      scl = one
      sumsq = amed
   end if
   dnrm2 = scl*sqrt( sumsq )
   return

end function

!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLARTG( F, G, C, S, R )
!
!       .. Scalar Arguments ..
!       REAL(wp)          C, F, G, R, S
!       ..
!
!  =============
!
!  Arguments:
!  ==========
!
!
!  Authors:
!  ========
!
!
!
!
!  ==================
!
!  =====================
!

subroutine dlartg( f, g, c, s, r )
   use la_constants, &
   only: wp=>dp, zero=>dzero, half=>dhalf, one=>done, &
         safmin=>dsafmin, safmax=>dsafmax
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     February 2021
!
!  .. Scalar Arguments ..
   real(wp) :: c, f, g, r, s
!  ..
!  .. Local Scalars ..
   real(wp) :: d, f1, fs, g1, gs, u, rtmin, rtmax
!  ..
!  .. Intrinsic Functions ..
   intrinsic :: abs, sign, sqrt
!  ..
!  .. Constants ..
   rtmin = sqrt( safmin )
   rtmax = sqrt( safmax/2 )
!  ..
!  .. Executable Statements ..
!
   f1 = abs( f )
   g1 = abs( g )
   if( g == zero ) then
      c = one
      s = zero
      r = f
   else if( f == zero ) then
      c = zero
      s = sign( one, g )
      r = g1
   else if( f1 > rtmin .and. f1 < rtmax .and. &
            g1 > rtmin .and. g1 < rtmax ) then
      d = sqrt( f*f + g*g )
      c = f1 / d
      r = sign( d, f )
      s = g / r
   else
      u = min( safmax, max( safmin, f1, g1 ) )
      fs = f / u
      gs = g / u
      d = sqrt( fs*fs + gs*gs )
      c = abs( fs ) / d
      r = sign( d, f )
      s = gs / r
      r = r*u
   end if
   return

end subroutine

!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Authors:
!  ========
!
!
!
!
!  ==================
!
!  =====================
!

!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLASSQ( N, X, INCX, SCALE, SUMSQ )
!
!       .. Scalar Arguments ..
!       INTEGER            INCX, N
!       DOUBLE PRECISION   SCALE, SUMSQ
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   X( * )
!       ..
!
!
!  =============
!
!  Arguments:
!  ==========
!
!
!  Authors:
!  ========
!
!
!  ==================
!
!  =====================
!
!
!  =====================================================================

subroutine dlassq( n, x, incx, scale, sumsq )
   use la_constants, &
      only: wp=>dp, zero=>dzero, one=>done, &
            sbig=>dsbig, ssml=>dssml, tbig=>dtbig, tsml=>dtsml
   use la_xisnan
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!  .. Scalar Arguments ..
   integer :: incx, n
   real(wp) :: scale, sumsq
!  ..
!  .. Array Arguments ..
   real(wp) :: x(*)
!  ..
!  .. Local Scalars ..
   integer :: i, ix
   logical :: notbig
   real(wp) :: abig, amed, asml, ax, ymax, ymin
!  ..
!
!  Quick return if possible
!
   if( la_isnan(scale) .or. la_isnan(sumsq) ) return
   if( sumsq == zero ) scale = one
   if( scale == zero ) then
      scale = one
      sumsq = zero
   end if
   if (n <= 0) then
      return
   end if
!
!  Compute the sum of squares in 3 accumulators:
!     abig -- sums of squares scaled down to avoid overflow
!     asml -- sums of squares scaled up to avoid underflow
!     amed -- sums of squares that do not require scaling
!  The thresholds and multipliers are
!     tbig -- values bigger than this are scaled down by sbig
!     tsml -- values smaller than this are scaled up by ssml
!
   notbig = .true.
   asml = zero
   amed = zero
   abig = zero
   ix = 1
   if( incx < 0 ) ix = 1 - (n-1)*incx
   do i = 1, n
      ax = abs(x(ix))
      if (ax > tbig) then
         abig = abig + (ax*sbig)**2
         notbig = .false.
      else if (ax < tsml) then
         if (notbig) asml = asml + (ax*ssml)**2
      else
         amed = amed + ax**2
      end if
      ix = ix + incx
   end do
!
!  Put the existing sum of squares into one of the accumulators
!
   if( sumsq > zero ) then
      ax = scale*sqrt( sumsq )
      if (ax > tbig) then
         if (scale > one) then
            scale = scale * sbig
            abig = abig + scale * (scale * sumsq)
         else
            ! sumsq > tbig^2 => (sbig * (sbig * sumsq)) is representable
            abig = abig + scale * (scale * (sbig * (sbig * sumsq)))
         end if
      else if (ax < tsml) then
         if (notbig) then
            if (scale < one) then
               scale = scale * ssml
               asml = asml + scale * (scale * sumsq)
            else
               ! sumsq < tsml^2 => (ssml * (ssml * sumsq)) is representable
               asml = asml + scale * (scale * (ssml * (ssml * sumsq)))
            end if
         end if
      else
         amed = amed + scale * (scale * sumsq)
      end if
   end if
!
!  Combine abig and amed or amed and asml if more than one
!  accumulator was used.
!
   if (abig > zero) then
!
!     Combine abig and amed if abig > 0.
!
      if (amed > zero .or. la_isnan(amed)) then
         abig = abig + (amed*sbig)*sbig
      end if
      scale = one / sbig
      sumsq = abig
   else if (asml > zero) then
!
!     Combine amed and asml if asml > 0.
!
      if (amed > zero .or. la_isnan(amed)) then
         amed = sqrt(amed)
         asml = sqrt(asml) / ssml
         if (asml > amed) then
            ymin = amed
            ymax = asml
         else
            ymin = asml
            ymax = amed
         end if
         scale = one
         sumsq = ymax**2*( one + (ymin/ymax)**2 )
      else
         scale = one / ssml
         sumsq = asml
      end if
   else
!
!     Otherwise all values are mid-range or zero
!
      scale = one
      sumsq = amed
   end if
   return

end subroutine

module lapack_in_fpp
use la_xisnan, disnanmod=>disnan
implicit none
private
public DGEEV



contains


!> \brief <b> DGEEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download DGEEV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgeev.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgeev.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgeev.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR,
!                         LDVR, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBVL, JOBVR
!       INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ),
!      $                   WI( * ), WORK( * ), WR( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGEEV computes for an N-by-N real nonsymmetric matrix A, the
!> eigenvalues and, optionally, the left and/or right eigenvectors.
!>
!> The right eigenvector v(j) of A satisfies
!>                  A * v(j) = lambda(j) * v(j)
!> where lambda(j) is its eigenvalue.
!> The left eigenvector u(j) of A satisfies
!>               u(j)**H * A = lambda(j) * u(j)**H
!> where u(j)**H denotes the conjugate-transpose of u(j).
!>
!> The computed eigenvectors are normalized to have Euclidean norm
!> equal to 1 and largest component real.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOBVL
!> \verbatim
!>          JOBVL is CHARACTER*1
!>          = 'N': left eigenvectors of A are not computed;
!>          = 'V': left eigenvectors of A are computed.
!> \endverbatim
!>
!> \param[in] JOBVR
!> \verbatim
!>          JOBVR is CHARACTER*1
!>          = 'N': right eigenvectors of A are not computed;
!>          = 'V': right eigenvectors of A are computed.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A. N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the N-by-N matrix A.
!>          On exit, A has been overwritten.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] WR
!> \verbatim
!>          WR is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] WI
!> \verbatim
!>          WI is DOUBLE PRECISION array, dimension (N)
!>          WR and WI contain the real and imaginary parts,
!>          respectively, of the computed eigenvalues.  Complex
!>          conjugate pairs of eigenvalues appear consecutively
!>          with the eigenvalue having the positive imaginary part
!>          first.
!> \endverbatim
!>
!> \param[out] VL
!> \verbatim
!>          VL is DOUBLE PRECISION array, dimension (LDVL,N)
!>          If JOBVL = 'V', the left eigenvectors u(j) are stored one
!>          after another in the columns of VL, in the same order
!>          as their eigenvalues.
!>          If JOBVL = 'N', VL is not referenced.
!>          If the j-th eigenvalue is real, then u(j) = VL(:,j),
!>          the j-th column of VL.
!>          If the j-th and (j+1)-st eigenvalues form a complex
!>          conjugate pair, then u(j) = VL(:,j) + i*VL(:,j+1) and
!>          u(j+1) = VL(:,j) - i*VL(:,j+1).
!> \endverbatim
!>
!> \param[in] LDVL
!> \verbatim
!>          LDVL is INTEGER
!>          The leading dimension of the array VL.  LDVL >= 1; if
!>          JOBVL = 'V', LDVL >= N.
!> \endverbatim
!>
!> \param[out] VR
!> \verbatim
!>          VR is DOUBLE PRECISION array, dimension (LDVR,N)
!>          If JOBVR = 'V', the right eigenvectors v(j) are stored one
!>          after another in the columns of VR, in the same order
!>          as their eigenvalues.
!>          If JOBVR = 'N', VR is not referenced.
!>          If the j-th eigenvalue is real, then v(j) = VR(:,j),
!>          the j-th column of VR.
!>          If the j-th and (j+1)-st eigenvalues form a complex
!>          conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and
!>          v(j+1) = VR(:,j) - i*VR(:,j+1).
!> \endverbatim
!>
!> \param[in] LDVR
!> \verbatim
!>          LDVR is INTEGER
!>          The leading dimension of the array VR.  LDVR >= 1; if
!>          JOBVR = 'V', LDVR >= N.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.  LWORK >= max(1,3*N), and
!>          if JOBVL = 'V' or JOBVR = 'V', LWORK >= 4*N.  For good
!>          performance, LWORK must generally be larger.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  if INFO = i, the QR algorithm failed to compute all the
!>                eigenvalues, and no eigenvectors have been computed;
!>                elements i+1:N of WR and WI contain eigenvalues which
!>                have converged.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!
!  @precisions fortran d -> s
!
!> \ingroup geev
!
!  =====================================================================

      SUBROUTINE dgeev( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, &
     &                  VR, &
     &                  LDVR, WORK, LWORK, INFO )
      implicit none
!
!  -- LAPACK driver routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          JOBVL, JOBVR
      INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ), &
     &                   WI( * ), WORK( * ), WR( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0d0, one = 1.0d0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY, SCALEA, WANTVL, WANTVR
      CHARACTER          SIDE
      INTEGER            HSWORK, I, IBAL, IERR, IHI, ILO, ITAU, IWRK, K, &
     &                   lwork_trevc, maxwrk, minwrk, nout
      DOUBLE PRECISION   ANRM, BIGNUM, CS, CSCALE, EPS, R, SCL, SMLNUM, &
     &                   SN
!     ..
!     .. Local Arrays ..
      LOGICAL            SELECT( 1 )
      DOUBLE PRECISION   DUM( 1 )
!     ..
!     .. External Subroutines ..
!      EXTERNAL           dgebak, dgebal, dgehrd, dhseqr, dlacpy, &
!     &                   dlartg, &
!     &                   dlascl, dorghr, drot, dscal, dtrevc3, xerbla
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      INTEGER            IDAMAX, ILAENV
!      DOUBLE PRECISION   DLAMCH, DLANGE, DLAPY2, DNRM2
!      EXTERNAL           lsame, idamax, ilaenv, dlamch, dlange, &
!     &                   dlapy2, &
!     &                   dnrm2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          max, sqrt
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      info = 0
      lquery = ( lwork.EQ.-1 )
      wantvl = lsame( jobvl, 'V' )
      wantvr = lsame( jobvr, 'V' )
      IF( ( .NOT.wantvl ) .AND. ( .NOT.lsame( jobvl, 'N' ) ) ) THEN
         info = -1
      ELSE IF( ( .NOT.wantvr ) .AND. &
     &         ( .NOT.lsame( jobvr, 'N' ) ) ) THEN
         info = -2
      ELSE IF( n.LT.0 ) THEN
         info = -3
      ELSE IF( lda.LT.max( 1, n ) ) THEN
         info = -5
      ELSE IF( ldvl.LT.1 .OR. ( wantvl .AND. ldvl.LT.n ) ) THEN
         info = -9
      ELSE IF( ldvr.LT.1 .OR. ( wantvr .AND. ldvr.LT.n ) ) THEN
         info = -11
      END IF
!
!     Compute workspace
!      (Note: Comments in the code beginning "Workspace:" describe the
!       minimal amount of workspace needed at that point in the code,
!       as well as the preferred amount for good performance.
!       NB refers to the optimal block size for the immediately
!       following subroutine, as returned by ILAENV.
!       HSWORK refers to the workspace preferred by DHSEQR, as
!       calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
!       the worst case.)
!
      IF( info.EQ.0 ) THEN
         IF( n.EQ.0 ) THEN
            minwrk = 1
            maxwrk = 1
         ELSE
            maxwrk = 2*n + n*ilaenv( 1, 'DGEHRD', ' ', n, 1, n, 0 )
            IF( wantvl ) THEN
               minwrk = 4*n
               maxwrk = max( maxwrk, 2*n + ( n - 1 )*ilaenv( 1, &
     &                       'DORGHR', ' ', n, 1, n, -1 ) )
               CALL dhseqr( 'S', 'V', n, 1, n, a, lda, wr, wi, vl, &
     &                      ldvl, &
     &                      work, -1, info )
               hswork = int( work(1) )
               maxwrk = max( maxwrk, n + 1, n + hswork )
               CALL dtrevc3( 'L', 'B', SELECT, n, a, lda, &
     &                       vl, ldvl, vr, ldvr, n, nout, &
     &                       work, -1, ierr )
               lwork_trevc = int( work(1) )
               maxwrk = max( maxwrk, n + lwork_trevc )
               maxwrk = max( maxwrk, 4*n )
            ELSE IF( wantvr ) THEN
               minwrk = 4*n
               maxwrk = max( maxwrk, 2*n + ( n - 1 )*ilaenv( 1, &
     &                       'DORGHR', ' ', n, 1, n, -1 ) )
               CALL dhseqr( 'S', 'V', n, 1, n, a, lda, wr, wi, vr, &
     &                      ldvr, &
     &                      work, -1, info )
               hswork = int( work(1) )
               maxwrk = max( maxwrk, n + 1, n + hswork )
               CALL dtrevc3( 'R', 'B', SELECT, n, a, lda, &
     &                       vl, ldvl, vr, ldvr, n, nout, &
     &                       work, -1, ierr )
               lwork_trevc = int( work(1) )
               maxwrk = max( maxwrk, n + lwork_trevc )
               maxwrk = max( maxwrk, 4*n )
            ELSE
               minwrk = 3*n
               CALL dhseqr( 'E', 'N', n, 1, n, a, lda, wr, wi, vr, &
     &                      ldvr, &
     &                      work, -1, info )
               hswork = int( work(1) )
               maxwrk = max( maxwrk, n + 1, n + hswork )
            END IF
            maxwrk = max( maxwrk, minwrk )
         END IF
         work( 1 ) = maxwrk
!
         IF( lwork.LT.minwrk .AND. .NOT.lquery ) THEN
            info = -13
         END IF
      END IF
!
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DGEEV ', -info )
         RETURN
      ELSE IF( lquery ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( n.EQ.0 ) &
     &   RETURN
!
!     Get machine constants
!
      eps = dlamch( 'P' )
      smlnum = dlamch( 'S' )
      bignum = one / smlnum
      smlnum = sqrt( smlnum ) / eps
      bignum = one / smlnum
!
!     Scale A if max element outside range [SMLNUM,BIGNUM]
!
      anrm = dlange( 'M', n, n, a, lda, dum )
      scalea = .false.
      IF( anrm.GT.zero .AND. anrm.LT.smlnum ) THEN
         scalea = .true.
         cscale = smlnum
      ELSE IF( anrm.GT.bignum ) THEN
         scalea = .true.
         cscale = bignum
      END IF
      IF( scalea ) &
     &   CALL dlascl( 'G', 0, 0, anrm, cscale, n, n, a, lda, ierr )
!
!     Balance the matrix
!     (Workspace: need N)
!
      ibal = 1
      CALL dgebal( 'B', n, a, lda, ilo, ihi, work( ibal ), ierr )
!
!     Reduce to upper Hessenberg form
!     (Workspace: need 3*N, prefer 2*N+N*NB)
!
      itau = ibal + n
      iwrk = itau + n
      CALL dgehrd( n, ilo, ihi, a, lda, work( itau ), work( iwrk ), &
     &             lwork-iwrk+1, ierr )
!
      IF( wantvl ) THEN
!
!        Want left eigenvectors
!        Copy Householder vectors to VL
!
         side = 'L'
         CALL dlacpy( 'L', n, n, a, lda, vl, ldvl )
!
!        Generate orthogonal matrix in VL
!        (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB)
!
         CALL dorghr( n, ilo, ihi, vl, ldvl, work( itau ), &
     &                work( iwrk ), &
     &                lwork-iwrk+1, ierr )
!
!        Perform QR iteration, accumulating Schur vectors in VL
!        (Workspace: need N+1, prefer N+HSWORK (see comments) )
!
         iwrk = itau
         CALL dhseqr( 'S', 'V', n, ilo, ihi, a, lda, wr, wi, vl, &
     &                ldvl, &
     &                work( iwrk ), lwork-iwrk+1, info )
!
         IF( wantvr ) THEN
!
!           Want left and right eigenvectors
!           Copy Schur vectors to VR
!
            side = 'B'
            CALL dlacpy( 'F', n, n, vl, ldvl, vr, ldvr )
         END IF
!
      ELSE IF( wantvr ) THEN
!
!        Want right eigenvectors
!        Copy Householder vectors to VR
!
         side = 'R'
         CALL dlacpy( 'L', n, n, a, lda, vr, ldvr )
!
!        Generate orthogonal matrix in VR
!        (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB)
!
         CALL dorghr( n, ilo, ihi, vr, ldvr, work( itau ), &
     &                work( iwrk ), &
     &                lwork-iwrk+1, ierr )
!
!        Perform QR iteration, accumulating Schur vectors in VR
!        (Workspace: need N+1, prefer N+HSWORK (see comments) )
!
         iwrk = itau
         CALL dhseqr( 'S', 'V', n, ilo, ihi, a, lda, wr, wi, vr, &
     &                ldvr, &
     &                work( iwrk ), lwork-iwrk+1, info )
!
      ELSE
!
!        Compute eigenvalues only
!        (Workspace: need N+1, prefer N+HSWORK (see comments) )
!
         iwrk = itau
         CALL dhseqr( 'E', 'N', n, ilo, ihi, a, lda, wr, wi, vr, &
     &                ldvr, &
     &                work( iwrk ), lwork-iwrk+1, info )
      END IF
!
!     If INFO .NE. 0 from DHSEQR, then quit
!
      IF( info.NE.0 ) &
     &   GO TO 50
!
      IF( wantvl .OR. wantvr ) THEN
!
!        Compute left and/or right eigenvectors
!        (Workspace: need 4*N, prefer N + N + 2*N*NB)
!
         CALL dtrevc3( side, 'B', SELECT, n, a, lda, vl, ldvl, vr, &
     &                 ldvr, &
     &                 n, nout, work( iwrk ), lwork-iwrk+1, ierr )
      END IF
!
      IF( wantvl ) THEN
!
!        Undo balancing of left eigenvectors
!        (Workspace: need N)
!
         CALL dgebak( 'B', 'L', n, ilo, ihi, work( ibal ), n, vl, &
     &                ldvl, &
     &                ierr )
!
!        Normalize left eigenvectors and make largest component real
!
         DO 20 i = 1, n
            IF( wi( i ).EQ.zero ) THEN
               scl = one / dnrm2( n, vl( 1, i ), 1 )
               CALL dscal( n, scl, vl( 1, i ), 1 )
            ELSE IF( wi( i ).GT.zero ) THEN
               scl = one / dlapy2( dnrm2( n, vl( 1, i ), 1 ), &
     &               dnrm2( n, vl( 1, i+1 ), 1 ) )
               CALL dscal( n, scl, vl( 1, i ), 1 )
               CALL dscal( n, scl, vl( 1, i+1 ), 1 )
               DO 10 k = 1, n
                  work( iwrk+k-1 ) = vl( k, i )**2 + vl( k, i+1 )**2
10             CONTINUE
               k = idamax( n, work( iwrk ), 1 )
               CALL dlartg( vl( k, i ), vl( k, i+1 ), cs, sn, r )
               CALL drot( n, vl( 1, i ), 1, vl( 1, i+1 ), 1, cs, sn )
               vl( k, i+1 ) = zero
            END IF
20       CONTINUE
      END IF
!
      IF( wantvr ) THEN
!
!        Undo balancing of right eigenvectors
!        (Workspace: need N)
!
         CALL dgebak( 'B', 'R', n, ilo, ihi, work( ibal ), n, vr, &
     &                ldvr, &
     &                ierr )
!
!        Normalize right eigenvectors and make largest component real
!
         DO 40 i = 1, n
            IF( wi( i ).EQ.zero ) THEN
               scl = one / dnrm2( n, vr( 1, i ), 1 )
               CALL dscal( n, scl, vr( 1, i ), 1 )
            ELSE IF( wi( i ).GT.zero ) THEN
               scl = one / dlapy2( dnrm2( n, vr( 1, i ), 1 ), &
     &               dnrm2( n, vr( 1, i+1 ), 1 ) )
               CALL dscal( n, scl, vr( 1, i ), 1 )
               CALL dscal( n, scl, vr( 1, i+1 ), 1 )
               DO 30 k = 1, n
                  work( iwrk+k-1 ) = vr( k, i )**2 + vr( k, i+1 )**2
30             CONTINUE
               k = idamax( n, work( iwrk ), 1 )
               CALL dlartg( vr( k, i ), vr( k, i+1 ), cs, sn, r )
               CALL drot( n, vr( 1, i ), 1, vr( 1, i+1 ), 1, cs, sn )
               vr( k, i+1 ) = zero
            END IF
40       CONTINUE
      END IF
!
!     Undo scaling if necessary
!
50    CONTINUE
      IF( scalea ) THEN
         CALL dlascl( 'G', 0, 0, cscale, anrm, n-info, 1, &
     &                wr( info+1 ), &
     &                max( n-info, 1 ), ierr )
         CALL dlascl( 'G', 0, 0, cscale, anrm, n-info, 1, &
     &                wi( info+1 ), &
     &                max( n-info, 1 ), ierr )
         IF( info.GT.0 ) THEN
            CALL dlascl( 'G', 0, 0, cscale, anrm, ilo-1, 1, wr, n, &
     &                   ierr )
            CALL dlascl( 'G', 0, 0, cscale, anrm, ilo-1, 1, wi, n, &
     &                   ierr )
         END IF
      END IF
!
      work( 1 ) = maxwrk
      RETURN
!
!     End of DGEEV
!

      END

!> \brief \b XERBLA
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download XERBLA + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/xerbla.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/xerbla.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/xerbla.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       SUBROUTINE XERBLA( SRNAME, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER*(*)      SRNAME
!       INTEGER            INFO
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> XERBLA  is an error handler for the LAPACK routines.
!> It is called by an LAPACK routine if an input parameter has an
!> invalid value.  A message is printed and execution stops.
!>
!> Installers may consider modifying the STOP statement in order to
!> call system-specific exception-handling facilities.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SRNAME
!> \verbatim
!>          SRNAME is CHARACTER*(*)
!>          The name of the routine which called XERBLA.
!> \endverbatim
!>
!> \param[in] INFO
!> \verbatim
!>          INFO is INTEGER
!>          The position of the invalid parameter in the parameter list
!>          of the calling routine.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup xerbla
!
!  =====================================================================

      SUBROUTINE xerbla( SRNAME, INFO )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER*(*)      SRNAME
      INTEGER            INFO
!     ..
!
! =====================================================================
!
!     .. Intrinsic Functions ..
      INTRINSIC          len_trim
!     ..
!     .. Executable Statements ..
!
      WRITE( *, fmt = 9999 )srname( 1:len_trim( srname ) ), info
!
      stop
!
9999  FORMAT( ' ** On entry to ', a, ' parameter number ', i2, ' had ', &
     &      'an illegal value' )
!
!     End of XERBLA
!

      END
!> \brief \b LSAME
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!      LOGICAL FUNCTION LSAME( CA, CB )
!
!     .. Scalar Arguments ..
!      CHARACTER          CA, CB
!     ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> LSAME returns .TRUE. if CA is the same letter as CB regardless of
!> case.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] CA
!> \param[in] CB
!> \verbatim
!>          CA and CB specify the single characters to be compared.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup lsame
!
!  =====================================================================

      LOGICAL FUNCTION lsame( CA, CB )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          ca, cb
!     ..
!
! =====================================================================
!
!     .. Intrinsic Functions ..
      INTRINSIC          ichar
!     ..
!     .. Local Scalars ..
      INTEGER            inta, intb, zcode
!     ..
!     .. Executable Statements ..
!
!     Test if the characters are equal
!
      lsame = ca.EQ.cb
      IF( lsame ) &
     &   RETURN
!
!     Now test for equivalence if both characters are alphabetic.
!
      zcode = ichar( 'Z' )
!
!     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
!     machines, on which ICHAR returns a value with bit 8 set.
!     ICHAR('A') on Prime machines returns 193 which is the same as
!     ICHAR('A') on an EBCDIC machine.
!
      inta = ichar( ca )
      intb = ichar( cb )
!
      IF( zcode.EQ.90 .OR. zcode.EQ.122 ) THEN
!
!        ASCII is assumed - ZCODE is the ASCII code of either lower or
!        upper case 'Z'.
!
         IF( inta.GE.97 .AND. inta.LE.122 ) inta = inta - 32
         IF( intb.GE.97 .AND. intb.LE.122 ) intb = intb - 32
!
      ELSE IF( zcode.EQ.233 .OR. zcode.EQ.169 ) THEN
!
!        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
!        upper case 'Z'.
!
         IF( inta.GE.129 .AND. inta.LE.137 .OR. &
     &       inta.GE.145 .AND. inta.LE.153 .OR. &
     &       inta.GE.162 .AND. inta.LE.169 ) inta = inta + 64
         IF( intb.GE.129 .AND. intb.LE.137 .OR. &
     &       intb.GE.145 .AND. intb.LE.153 .OR. &
     &       intb.GE.162 .AND. intb.LE.169 ) intb = intb + 64
!
      ELSE IF( zcode.EQ.218 .OR. zcode.EQ.250 ) THEN
!
!        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
!        plus 128 of either lower or upper case 'Z'.
!
         IF( inta.GE.225 .AND. inta.LE.250 ) inta = inta - 32
         IF( intb.GE.225 .AND. intb.LE.250 ) intb = intb - 32
      END IF
      lsame = inta.EQ.intb
!
!     RETURN
!
!     End of LSAME
!

      END
!> \brief \b ILAENV
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download ILAENV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilaenv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ilaenv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilaenv.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
!
!       .. Scalar Arguments ..
!       CHARACTER*( * )    NAME, OPTS
!       INTEGER            ISPEC, N1, N2, N3, N4
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ILAENV is called from the LAPACK routines to choose problem-dependent
!> parameters for the local environment.  See ISPEC for a description of
!> the parameters.
!>
!> ILAENV returns an INTEGER
!> if ILAENV >= 0: ILAENV returns the value of the parameter specified by ISPEC
!> if ILAENV < 0:  if ILAENV = -k, the k-th argument had an illegal value.
!>
!> This version provides a set of parameters which should give good,
!> but not optimal, performance on many of the currently available
!> computers.  Users are encouraged to modify this subroutine to set
!> the tuning parameters for their particular machine using the option
!> and problem size information in the arguments.
!>
!> This routine will not function correctly if it is converted to all
!> lower case.  Converting it to all upper case is allowed.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ISPEC
!> \verbatim
!>          ISPEC is INTEGER
!>          Specifies the parameter to be returned as the value of
!>          ILAENV.
!>          = 1: the optimal blocksize; if this value is 1, an unblocked
!>               algorithm will give the best performance.
!>          = 2: the minimum block size for which the block routine
!>               should be used; if the usable block size is less than
!>               this value, an unblocked routine should be used.
!>          = 3: the crossover point (in a block routine, for N less
!>               than this value, an unblocked routine should be used)
!>          = 4: the number of shifts, used in the nonsymmetric
!>               eigenvalue routines (DEPRECATED)
!>          = 5: the minimum column dimension for blocking to be used;
!>               rectangular blocks must have dimension at least k by m,
!>               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
!>          = 6: the crossover point for the SVD (when reducing an m by n
!>               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
!>               this value, a QR factorization is used first to reduce
!>               the matrix to a triangular form.)
!>          = 7: the number of processors
!>          = 8: the crossover point for the multishift QR method
!>               for nonsymmetric eigenvalue problems (DEPRECATED)
!>          = 9: maximum size of the subproblems at the bottom of the
!>               computation tree in the divide-and-conquer algorithm
!>               (used by xGELSD and xGESDD)
!>          =10: ieee infinity and NaN arithmetic can be trusted not to trap
!>          =11: infinity arithmetic can be trusted not to trap
!>          12 <= ISPEC <= 17:
!>               xHSEQR or related subroutines,
!>               see IPARMQ for detailed explanation
!> \endverbatim
!>
!> \param[in] NAME
!> \verbatim
!>          NAME is CHARACTER*(*)
!>          The name of the calling subroutine, in either upper case or
!>          lower case.
!> \endverbatim
!>
!> \param[in] OPTS
!> \verbatim
!>          OPTS is CHARACTER*(*)
!>          The character options to the subroutine NAME, concatenated
!>          into a single character string.  For example, UPLO = 'U',
!>          TRANS = 'T', and DIAG = 'N' for a triangular routine would
!>          be specified as OPTS = 'UTN'.
!> \endverbatim
!>
!> \param[in] N1
!> \verbatim
!>          N1 is INTEGER
!> \endverbatim
!>
!> \param[in] N2
!> \verbatim
!>          N2 is INTEGER
!> \endverbatim
!>
!> \param[in] N3
!> \verbatim
!>          N3 is INTEGER
!> \endverbatim
!>
!> \param[in] N4
!> \verbatim
!>          N4 is INTEGER
!>          Problem dimensions for the subroutine NAME; these may not all
!>          be required.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup ilaenv
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The following conventions have been used when calling ILAENV from the
!>  LAPACK routines:
!>  1)  OPTS is a concatenation of all of the character options to
!>      subroutine NAME, in the same order that they appear in the
!>      argument list for NAME, even if they are not used in determining
!>      the value of the parameter specified by ISPEC.
!>  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
!>      that they appear in the argument list for NAME.  N1 is used
!>      first, N2 second, and so on, and unused problem dimensions are
!>      passed a value of -1.
!>  3)  The parameter value returned by ILAENV is checked for validity in
!>      the calling subroutine.  For example, ILAENV is used to retrieve
!>      the optimal blocksize for STRTRI as follows:
!>
!>      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
!>      IF( NB.LE.1 ) NB = MAX( 1, N )
!> \endverbatim
!>
!  =====================================================================

      INTEGER FUNCTION ilaenv( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER*( * )    name, opts
      INTEGER            ispec, n1, n2, n3, n4
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER            i, ic, iz, nb, nbmin, nx
      LOGICAL            cname, sname, twostage
      CHARACTER          c1*1, c2*2, c4*2, c3*3, subnam*16
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          char, ichar, int, min, real
!     ..
!     .. External Functions ..
!      INTEGER            ieeeck, iparmq, iparam2stage
!      EXTERNAL           ieeeck, iparmq, iparam2stage
!     ..
!     .. Executable Statements ..
!
      GO TO ( 10, 10, 10, 80, 90, 100, 110, 120, &
     &        130, 140, 150, 160, 160, 160, 160, 160, 160)ispec
!
!     Invalid value for ISPEC
!
      ilaenv = -1
      RETURN
!
10    CONTINUE
!
!     Convert NAME to upper case if the first character is lower case.
!
      ilaenv = 1
      subnam = name
      ic = ichar( subnam( 1: 1 ) )
      iz = ichar( 'Z' )
      IF( iz.EQ.90 .OR. iz.EQ.122 ) THEN
!
!        ASCII character set
!
         IF( ic.GE.97 .AND. ic.LE.122 ) THEN
            subnam( 1: 1 ) = char( ic-32 )
            DO 20 i = 2, 6
               ic = ichar( subnam( i: i ) )
               IF( ic.GE.97 .AND. ic.LE.122 ) &
     &            subnam( i: i ) = char( ic-32 )
20          CONTINUE
         END IF
!
      ELSE IF( iz.EQ.233 .OR. iz.EQ.169 ) THEN
!
!        EBCDIC character set
!
         IF( ( ic.GE.129 .AND. ic.LE.137 ) .OR. &
     &       ( ic.GE.145 .AND. ic.LE.153 ) .OR. &
     &       ( ic.GE.162 .AND. ic.LE.169 ) ) THEN
            subnam( 1: 1 ) = char( ic+64 )
            DO 30 i = 2, 6
               ic = ichar( subnam( i: i ) )
               IF( ( ic.GE.129 .AND. ic.LE.137 ) .OR. &
     &             ( ic.GE.145 .AND. ic.LE.153 ) .OR. &
     &             ( ic.GE.162 .AND. ic.LE.169 ) )subnam( i: &
     &             i ) = char( ic+64 )
30          CONTINUE
         END IF
!
      ELSE IF( iz.EQ.218 .OR. iz.EQ.250 ) THEN
!
!        Prime machines:  ASCII+128
!
         IF( ic.GE.225 .AND. ic.LE.250 ) THEN
            subnam( 1: 1 ) = char( ic-32 )
            DO 40 i = 2, 6
               ic = ichar( subnam( i: i ) )
               IF( ic.GE.225 .AND. ic.LE.250 ) &
     &            subnam( i: i ) = char( ic-32 )
40          CONTINUE
         END IF
      END IF
!
      c1 = subnam( 1: 1 )
      sname = c1.EQ.'S' .OR. c1.EQ.'D'
      cname = c1.EQ.'C' .OR. c1.EQ.'Z'
      IF( .NOT.( cname .OR. sname ) ) &
     &   RETURN
      c2 = subnam( 2: 3 )
      c3 = subnam( 4: 6 )
      c4 = c3( 2: 3 )
      twostage = len( subnam ).GE.11 &
     &           .AND. subnam( 11: 11 ).EQ.'2'
!
      GO TO ( 50, 60, 70 )ispec
!
50    CONTINUE
!
!     ISPEC = 1:  block size
!
!     In these examples, separate code is provided for setting NB for
!     real and complex.  We assume that NB will take the same value in
!     single or double precision.
!
      nb = 1
!
      IF( subnam(2:6).EQ.'LAORH' ) THEN
!
!        This is for *LAORHR_GETRFNP routine
!
         IF( sname ) THEN
             nb = 32
         ELSE
             nb = 32
         END IF
      ELSE IF( c2.EQ.'GE' ) THEN
         IF( c3.EQ.'TRF' ) THEN
            IF( sname ) THEN
               nb = 64
            ELSE
               nb = 64
            END IF
         ELSE IF( c3.EQ.'QRF' .OR. c3.EQ.'RQF' .OR. c3.EQ.'LQF' .OR. &
     &            c3.EQ.'QLF' ) THEN
            IF( sname ) THEN
               nb = 32
            ELSE
               nb = 32
            END IF
         ELSE IF( c3.EQ.'QR ') THEN
            IF( n3 .EQ. 1) THEN
               IF( sname ) THEN
!     M*N
                  IF ((n1*n2.LE.131072).OR.(n1.LE.8192)) THEN
                     nb = n1
                  ELSE
                     nb = 32768/n2
                  END IF
               ELSE
                  IF ((n1*n2.LE.131072).OR.(n1.LE.8192)) THEN
                     nb = n1
                  ELSE
                     nb = 32768/n2
                  END IF
               END IF
            ELSE
               IF( sname ) THEN
                  nb = 1
               ELSE
                  nb = 1
               END IF
            END IF
         ELSE IF( c3.EQ.'LQ ') THEN
            IF( n3 .EQ. 2) THEN
               IF( sname ) THEN
!     M*N
                  IF ((n1*n2.LE.131072).OR.(n1.LE.8192)) THEN
                     nb = n1
                  ELSE
                     nb = 32768/n2
                  END IF
               ELSE
                  IF ((n1*n2.LE.131072).OR.(n1.LE.8192)) THEN
                     nb = n1
                  ELSE
                     nb = 32768/n2
                  END IF
               END IF
            ELSE
               IF( sname ) THEN
                  nb = 1
               ELSE
                  nb = 1
               END IF
            END IF
         ELSE IF( c3.EQ.'HRD' ) THEN
            IF( sname ) THEN
               nb = 32
            ELSE
               nb = 32
            END IF
         ELSE IF( c3.EQ.'BRD' ) THEN
            IF( sname ) THEN
               nb = 32
            ELSE
               nb = 32
            END IF
         ELSE IF( c3.EQ.'TRI' ) THEN
            IF( sname ) THEN
               nb = 64
            ELSE
               nb = 64
            END IF
         ELSE IF( subnam( 4: 7 ).EQ.'QP3RK' ) THEN
            IF( sname ) THEN
               nb = 32
            ELSE
               nb = 32
            END IF
         END IF
      ELSE IF( c2.EQ.'PO' ) THEN
         IF( c3.EQ.'TRF' ) THEN
            IF( sname ) THEN
               nb = 64
            ELSE
               nb = 64
            END IF
         END IF
      ELSE IF( c2.EQ.'SY' ) THEN
         IF( c3.EQ.'TRF' ) THEN
            IF( sname ) THEN
               IF( twostage ) THEN
                  nb = 192
               ELSE
                  nb = 64
               END IF
            ELSE
               IF( twostage ) THEN
                  nb = 192
               ELSE
                  nb = 64
               END IF
            END IF
         ELSE IF( sname .AND. c3.EQ.'TRD' ) THEN
            nb = 32
         ELSE IF( sname .AND. c3.EQ.'GST' ) THEN
            nb = 64
         END IF
      ELSE IF( cname .AND. c2.EQ.'HE' ) THEN
         IF( c3.EQ.'TRF' ) THEN
            IF( twostage ) THEN
               nb = 192
            ELSE
               nb = 64
            END IF
         ELSE IF( c3.EQ.'TRD' ) THEN
            nb = 32
         ELSE IF( c3.EQ.'GST' ) THEN
            nb = 64
         END IF
      ELSE IF( sname .AND. c2.EQ.'OR' ) THEN
         IF( c3( 1: 1 ).EQ.'G' ) THEN
            IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ. &
     &          'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' ) &
     &           THEN
               nb = 32
            END IF
         ELSE IF( c3( 1: 1 ).EQ.'M' ) THEN
            IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ. &
     &          'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' ) &
     &           THEN
               nb = 32
            END IF
         END IF
      ELSE IF( cname .AND. c2.EQ.'UN' ) THEN
         IF( c3( 1: 1 ).EQ.'G' ) THEN
            IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ. &
     &          'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' ) &
     &           THEN
               nb = 32
            END IF
         ELSE IF( c3( 1: 1 ).EQ.'M' ) THEN
            IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ. &
     &          'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' ) &
     &           THEN
               nb = 32
            END IF
         END IF
      ELSE IF( c2.EQ.'GB' ) THEN
         IF( c3.EQ.'TRF' ) THEN
            IF( sname ) THEN
               IF( n4.LE.64 ) THEN
                  nb = 1
               ELSE
                  nb = 32
               END IF
            ELSE
               IF( n4.LE.64 ) THEN
                  nb = 1
               ELSE
                  nb = 32
               END IF
            END IF
         END IF
      ELSE IF( c2.EQ.'PB' ) THEN
         IF( c3.EQ.'TRF' ) THEN
            IF( sname ) THEN
               IF( n2.LE.64 ) THEN
                  nb = 1
               ELSE
                  nb = 32
               END IF
            ELSE
               IF( n2.LE.64 ) THEN
                  nb = 1
               ELSE
                  nb = 32
               END IF
            END IF
         END IF
      ELSE IF( c2.EQ.'TR' ) THEN
         IF( c3.EQ.'TRI' ) THEN
            IF( sname ) THEN
               nb = 64
            ELSE
               nb = 64
            END IF
         ELSE IF ( c3.EQ.'EVC' ) THEN
            IF( sname ) THEN
               nb = 64
            ELSE
               nb = 64
            END IF
         ELSE IF( c3.EQ.'SYL' ) THEN
!           The upper bound is to prevent overly aggressive scaling.
            IF( sname ) THEN
               nb = min( max( 48, int( ( min( n1, n2 ) * 16 ) / 100) ), &
     &                   240 )
            ELSE
               nb = min( max( 24, int( ( min( n1, n2 ) * 8 ) / 100) ), &
     &                   80 )
            END IF
         END IF
      ELSE IF( c2.EQ.'LA' ) THEN
         IF( c3.EQ.'UUM' ) THEN
            IF( sname ) THEN
               nb = 64
            ELSE
               nb = 64
            END IF
         ELSE IF( c3.EQ.'TRS' ) THEN
            IF( sname ) THEN
               nb = 32
            ELSE
               nb = 32
            END IF
         END IF
      ELSE IF( sname .AND. c2.EQ.'ST' ) THEN
         IF( c3.EQ.'EBZ' ) THEN
            nb = 1
         END IF
      ELSE IF( c2.EQ.'GG' ) THEN
         nb = 32
         IF( c3.EQ.'HD3' ) THEN
            IF( sname ) THEN
               nb = 32
            ELSE
               nb = 32
            END IF
         END IF
      END IF
      ilaenv = nb
      RETURN
!
60    CONTINUE
!
!     ISPEC = 2:  minimum block size
!
      nbmin = 2
      IF( c2.EQ.'GE' ) THEN
         IF( c3.EQ.'QRF' .OR. c3.EQ.'RQF' .OR. c3.EQ.'LQF' .OR. c3.EQ. &
     &       'QLF' ) THEN
            IF( sname ) THEN
               nbmin = 2
            ELSE
               nbmin = 2
            END IF
         ELSE IF( c3.EQ.'HRD' ) THEN
            IF( sname ) THEN
               nbmin = 2
            ELSE
               nbmin = 2
            END IF
         ELSE IF( c3.EQ.'BRD' ) THEN
            IF( sname ) THEN
               nbmin = 2
            ELSE
               nbmin = 2
            END IF
         ELSE IF( c3.EQ.'TRI' ) THEN
            IF( sname ) THEN
               nbmin = 2
            ELSE
               nbmin = 2
            END IF
         ELSE IF( subnam( 4: 7 ).EQ.'QP3RK' ) THEN
            IF( sname ) THEN
               nbmin = 2
            ELSE
               nbmin = 2
            END IF
         END IF

      ELSE IF( c2.EQ.'SY' ) THEN
         IF( c3.EQ.'TRF' ) THEN
            IF( sname ) THEN
               nbmin = 8
            ELSE
               nbmin = 8
            END IF
         ELSE IF( sname .AND. c3.EQ.'TRD' ) THEN
            nbmin = 2
         END IF
      ELSE IF( cname .AND. c2.EQ.'HE' ) THEN
         IF( c3.EQ.'TRD' ) THEN
            nbmin = 2
         END IF
      ELSE IF( sname .AND. c2.EQ.'OR' ) THEN
         IF( c3( 1: 1 ).EQ.'G' ) THEN
            IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ. &
     &          'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' ) &
     &           THEN
               nbmin = 2
            END IF
         ELSE IF( c3( 1: 1 ).EQ.'M' ) THEN
            IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ. &
     &          'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' ) &
     &           THEN
               nbmin = 2
            END IF
         END IF
      ELSE IF( cname .AND. c2.EQ.'UN' ) THEN
         IF( c3( 1: 1 ).EQ.'G' ) THEN
            IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ. &
     &          'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' ) &
     &           THEN
               nbmin = 2
            END IF
         ELSE IF( c3( 1: 1 ).EQ.'M' ) THEN
            IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ. &
     &          'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' ) &
     &           THEN
               nbmin = 2
            END IF
         END IF
      ELSE IF( c2.EQ.'GG' ) THEN
         nbmin = 2
         IF( c3.EQ.'HD3' ) THEN
            nbmin = 2
         END IF
      END IF
      ilaenv = nbmin
      RETURN
!
70    CONTINUE
!
!     ISPEC = 3:  crossover point
!
      nx = 0
      IF( c2.EQ.'GE' ) THEN
         IF( c3.EQ.'QRF' .OR. c3.EQ.'RQF' .OR. c3.EQ.'LQF' .OR. c3.EQ. &
     &       'QLF' ) THEN
            IF( sname ) THEN
               nx = 128
            ELSE
               nx = 128
            END IF
         ELSE IF( c3.EQ.'HRD' ) THEN
            IF( sname ) THEN
               nx = 128
            ELSE
               nx = 128
            END IF
         ELSE IF( c3.EQ.'BRD' ) THEN
            IF( sname ) THEN
               nx = 128
            ELSE
               nx = 128
            END IF
         ELSE IF( subnam( 4: 7 ).EQ.'QP3RK' ) THEN
            IF( sname ) THEN
               nx = 128
            ELSE
               nx = 128
            END IF
         END IF
      ELSE IF( c2.EQ.'SY' ) THEN
         IF( sname .AND. c3.EQ.'TRD' ) THEN
            nx = 32
         END IF
      ELSE IF( cname .AND. c2.EQ.'HE' ) THEN
         IF( c3.EQ.'TRD' ) THEN
            nx = 32
         END IF
      ELSE IF( sname .AND. c2.EQ.'OR' ) THEN
         IF( c3( 1: 1 ).EQ.'G' ) THEN
            IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ. &
     &          'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' ) &
     &           THEN
               nx = 128
            END IF
         END IF
      ELSE IF( cname .AND. c2.EQ.'UN' ) THEN
         IF( c3( 1: 1 ).EQ.'G' ) THEN
            IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ. &
     &          'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' ) &
     &           THEN
               nx = 128
            END IF
         END IF
      ELSE IF( c2.EQ.'GG' ) THEN
         nx = 128
         IF( c3.EQ.'HD3' ) THEN
            nx = 128
         END IF
      END IF
      ilaenv = nx
      RETURN
!
80    CONTINUE
!
!     ISPEC = 4:  number of shifts (used by xHSEQR)
!
      ilaenv = 6
      RETURN
!
90    CONTINUE
!
!     ISPEC = 5:  minimum column dimension (not used)
!
      ilaenv = 2
      RETURN
!
100   CONTINUE
!
!     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
!
      ilaenv = int( real( min( n1, n2 ) )*1.6e0 )
      RETURN
!
110   CONTINUE
!
!     ISPEC = 7:  number of processors (not used)
!
      ilaenv = 1
      RETURN
!
120   CONTINUE
!
!     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
!
      ilaenv = 50
      RETURN
!
130   CONTINUE
!
!     ISPEC = 9:  maximum size of the subproblems at the bottom of the
!                 computation tree in the divide-and-conquer algorithm
!                 (used by xGELSD and xGESDD)
!
      ilaenv = 25
      RETURN
!
140   CONTINUE
!
!     ISPEC = 10: ieee and infinity NaN arithmetic can be trusted not to trap
!
!     ILAENV = 0
      ilaenv = 1
      IF( ilaenv.EQ.1 ) THEN
         ilaenv = ieeeck( 1, 0.0, 1.0 )
      END IF
      RETURN
!
150   CONTINUE
!
!     ISPEC = 11: ieee infinity arithmetic can be trusted not to trap
!
!     ILAENV = 0
      ilaenv = 1
      IF( ilaenv.EQ.1 ) THEN
         ilaenv = ieeeck( 0, 0.0, 1.0 )
      END IF
      RETURN
!
160   CONTINUE
!
!     12 <= ISPEC <= 17: xHSEQR or related subroutines.
!
      ilaenv = iparmq( ispec, name, opts, n1, n2, n3, n4 )
      RETURN
!
!     End of ILAENV
!

      END
!> \brief \b IPARMQ
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download IPARMQ + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/iparmq.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/iparmq.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/iparmq.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       INTEGER FUNCTION IPARMQ( ISPEC, NAME, OPTS, N, ILO, IHI, LWORK )
!
!       .. Scalar Arguments ..
!       INTEGER            IHI, ILO, ISPEC, LWORK, N
!       CHARACTER          NAME*( * ), OPTS*( * )
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>      This program sets problem and machine dependent parameters
!>      useful for xHSEQR and related subroutines for eigenvalue
!>      problems. It is called whenever
!>      IPARMQ is called with 12 <= ISPEC <= 16
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ISPEC
!> \verbatim
!>          ISPEC is INTEGER
!>              ISPEC specifies which tunable parameter IPARMQ should
!>              return.
!>
!>              ISPEC=12: (INMIN)  Matrices of order nmin or less
!>                        are sent directly to xLAHQR, the implicit
!>                        double shift QR algorithm.  NMIN must be
!>                        at least 11.
!>
!>              ISPEC=13: (INWIN)  Size of the deflation window.
!>                        This is best set greater than or equal to
!>                        the number of simultaneous shifts NS.
!>                        Larger matrices benefit from larger deflation
!>                        windows.
!>
!>              ISPEC=14: (INIBL) Determines when to stop nibbling and
!>                        invest in an (expensive) multi-shift QR sweep.
!>                        If the aggressive early deflation subroutine
!>                        finds LD converged eigenvalues from an order
!>                        NW deflation window and LD > (NW*NIBBLE)/100,
!>                        then the next QR sweep is skipped and early
!>                        deflation is applied immediately to the
!>                        remaining active diagonal block.  Setting
!>                        IPARMQ(ISPEC=14) = 0 causes TTQRE to skip a
!>                        multi-shift QR sweep whenever early deflation
!>                        finds a converged eigenvalue.  Setting
!>                        IPARMQ(ISPEC=14) greater than or equal to 100
!>                        prevents TTQRE from skipping a multi-shift
!>                        QR sweep.
!>
!>              ISPEC=15: (NSHFTS) The number of simultaneous shifts in
!>                        a multi-shift QR iteration.
!>
!>              ISPEC=16: (IACC22) IPARMQ is set to 0, 1 or 2 with the
!>                        following meanings.
!>                        0:  During the multi-shift QR/QZ sweep,
!>                            blocked eigenvalue reordering, blocked
!>                            Hessenberg-triangular reduction,
!>                            reflections and/or rotations are not
!>                            accumulated when updating the
!>                            far-from-diagonal matrix entries.
!>                        1:  During the multi-shift QR/QZ sweep,
!>                            blocked eigenvalue reordering, blocked
!>                            Hessenberg-triangular reduction,
!>                            reflections and/or rotations are
!>                            accumulated, and matrix-matrix
!>                            multiplication is used to update the
!>                            far-from-diagonal matrix entries.
!>                        2:  During the multi-shift QR/QZ sweep,
!>                            blocked eigenvalue reordering, blocked
!>                            Hessenberg-triangular reduction,
!>                            reflections and/or rotations are
!>                            accumulated, and 2-by-2 block structure
!>                            is exploited during matrix-matrix
!>                            multiplies.
!>                        (If xTRMM is slower than xGEMM, then
!>                        IPARMQ(ISPEC=16)=1 may be more efficient than
!>                        IPARMQ(ISPEC=16)=2 despite the greater level of
!>                        arithmetic work implied by the latter choice.)
!>
!>              ISPEC=17: (ICOST) An estimate of the relative cost of flops
!>                        within the near-the-diagonal shift chase compared
!>                        to flops within the BLAS calls of a QZ sweep.
!> \endverbatim
!>
!> \param[in] NAME
!> \verbatim
!>          NAME is CHARACTER string
!>               Name of the calling subroutine
!> \endverbatim
!>
!> \param[in] OPTS
!> \verbatim
!>          OPTS is CHARACTER string
!>               This is a concatenation of the string arguments to
!>               TTQRE.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>               N is the order of the Hessenberg matrix H.
!> \endverbatim
!>
!> \param[in] ILO
!> \verbatim
!>          ILO is INTEGER
!> \endverbatim
!>
!> \param[in] IHI
!> \verbatim
!>          IHI is INTEGER
!>               It is assumed that H is already upper triangular
!>               in rows and columns 1:ILO-1 and IHI+1:N.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>               The amount of workspace available.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup iparmq
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>       Little is known about how best to choose these parameters.
!>       It is possible to use different values of the parameters
!>       for each of CHSEQR, DHSEQR, SHSEQR and ZHSEQR.
!>
!>       It is probably best to choose different parameters for
!>       different matrices and different parameters at different
!>       times during the iteration, but this has not been
!>       implemented --- yet.
!>
!>
!>       The best choices of most of the parameters depend
!>       in an ill-understood way on the relative execution
!>       rate of xLAQR3 and xLAQR5 and on the nature of each
!>       particular eigenvalue problem.  Experiment may be the
!>       only practical way to determine which choices are most
!>       effective.
!>
!>       Following is a list of default values supplied by IPARMQ.
!>       These defaults may be adjusted in order to attain better
!>       performance in any particular computational environment.
!>
!>       IPARMQ(ISPEC=12) The xLAHQR vs xLAQR0 crossover point.
!>                        Default: 75. (Must be at least 11.)
!>
!>       IPARMQ(ISPEC=13) Recommended deflation window size.
!>                        This depends on ILO, IHI and NS, the
!>                        number of simultaneous shifts returned
!>                        by IPARMQ(ISPEC=15).  The default for
!>                        (IHI-ILO+1) <= 500 is NS.  The default
!>                        for (IHI-ILO+1) > 500 is 3*NS/2.
!>
!>       IPARMQ(ISPEC=14) Nibble crossover point.  Default: 14.
!>
!>       IPARMQ(ISPEC=15) Number of simultaneous shifts, NS.
!>                        a multi-shift QR iteration.
!>
!>                        If IHI-ILO+1 is ...
!>
!>                        greater than      ...but less    ... the
!>                        or equal to ...      than        default is
!>
!>                                0               30       NS =   2+
!>                               30               60       NS =   4+
!>                               60              150       NS =  10
!>                              150              590       NS =  **
!>                              590             3000       NS =  64
!>                             3000             6000       NS = 128
!>                             6000             infinity   NS = 256
!>
!>                    (+)  By default matrices of this order are
!>                         passed to the implicit double shift routine
!>                         xLAHQR.  See IPARMQ(ISPEC=12) above.   These
!>                         values of NS are used only in case of a rare
!>                         xLAHQR failure.
!>
!>                    (**) The asterisks (**) indicate an ad-hoc
!>                         function increasing from 10 to 64.
!>
!>       IPARMQ(ISPEC=16) Select structured matrix multiply.
!>                        (See ISPEC=16 above for details.)
!>                        Default: 3.
!>
!>       IPARMQ(ISPEC=17) Relative cost heuristic for blocksize selection.
!>                        Expressed as a percentage.
!>                        Default: 10.
!> \endverbatim
!>
!  =====================================================================

      INTEGER FUNCTION iparmq( ISPEC, NAME, OPTS, N, ILO, IHI, &
     &                         LWORK )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            ihi, ilo, ispec, lwork, n
      CHARACTER          name*( * ), opts*( * )
!
!  ================================================================
!     .. Parameters ..
      INTEGER            inmin, inwin, inibl, ishfts, iacc22, icost
      parameter( inmin = 12, inwin = 13, inibl = 14, &
     &                   ishfts = 15, iacc22 = 16, icost = 17 )
      INTEGER            nmin, k22min, kacmin, nibble, knwswp, rcost
      parameter( nmin = 75, k22min = 14, kacmin = 14, &
     &                   nibble = 14, knwswp = 500, rcost = 10 )
      REAL               two
      parameter( two = 2.0 )
!     ..
!     .. Local Scalars ..
      INTEGER            nh, ns
      INTEGER            i, ic, iz
      CHARACTER          subnam*6
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          log, max, mod, nint, real
!     ..
!     .. Executable Statements ..
      IF( ( ispec.EQ.ishfts ) .OR. ( ispec.EQ.inwin ) .OR. &
     &    ( ispec.EQ.iacc22 ) ) THEN
!
!        ==== Set the number simultaneous shifts ====
!
         nh = ihi - ilo + 1
         ns = 2
         IF( nh.GE.30 ) &
     &      ns = 4
         IF( nh.GE.60 ) &
     &      ns = 10
         IF( nh.GE.150 ) &
     &      ns = max( 10, nh / nint( log( real( nh ) ) / log( two ) ) )
         IF( nh.GE.590 ) &
     &      ns = 64
         IF( nh.GE.3000 ) &
     &      ns = 128
         IF( nh.GE.6000 ) &
     &      ns = 256
         ns = max( 2, ns-mod( ns, 2 ) )
      END IF
!
      IF( ispec.EQ.inmin ) THEN
!
!
!        ===== Matrices of order smaller than NMIN get sent
!        .     to xLAHQR, the classic double shift algorithm.
!        .     This must be at least 11. ====
!
         iparmq = nmin
!
      ELSE IF( ispec.EQ.inibl ) THEN
!
!        ==== INIBL: skip a multi-shift qr iteration and
!        .    whenever aggressive early deflation finds
!        .    at least (NIBBLE*(window size)/100) deflations. ====
!
         iparmq = nibble
!
      ELSE IF( ispec.EQ.ishfts ) THEN
!
!        ==== NSHFTS: The number of simultaneous shifts =====
!
         iparmq = ns
!
      ELSE IF( ispec.EQ.inwin ) THEN
!
!        ==== NW: deflation window size.  ====
!
         IF( nh.LE.knwswp ) THEN
            iparmq = ns
         ELSE
            iparmq = 3*ns / 2
         END IF
!
      ELSE IF( ispec.EQ.iacc22 ) THEN
!
!        ==== IACC22: Whether to accumulate reflections
!        .     before updating the far-from-diagonal elements
!        .     and whether to use 2-by-2 block structure while
!        .     doing it.  A small amount of work could be saved
!        .     by making this choice dependent also upon the
!        .     NH=IHI-ILO+1.
!
!
!        Convert NAME to upper case if the first character is lower case.
!
         iparmq = 0
         subnam = name
         ic = ichar( subnam( 1: 1 ) )
         iz = ichar( 'Z' )
         IF( iz.EQ.90 .OR. iz.EQ.122 ) THEN
!
!           ASCII character set
!
            IF( ic.GE.97 .AND. ic.LE.122 ) THEN
               subnam( 1: 1 ) = char( ic-32 )
               DO i = 2, 6
                  ic = ichar( subnam( i: i ) )
                  IF( ic.GE.97 .AND. ic.LE.122 ) &
     &               subnam( i: i ) = char( ic-32 )
               END DO
            END IF
!
         ELSE IF( iz.EQ.233 .OR. iz.EQ.169 ) THEN
!
!           EBCDIC character set
!
            IF( ( ic.GE.129 .AND. ic.LE.137 ) .OR. &
     &          ( ic.GE.145 .AND. ic.LE.153 ) .OR. &
     &          ( ic.GE.162 .AND. ic.LE.169 ) ) THEN
               subnam( 1: 1 ) = char( ic+64 )
               DO i = 2, 6
                  ic = ichar( subnam( i: i ) )
                  IF( ( ic.GE.129 .AND. ic.LE.137 ) .OR. &
     &                ( ic.GE.145 .AND. ic.LE.153 ) .OR. &
     &                ( ic.GE.162 .AND. ic.LE.169 ) )subnam( i: &
     &                i ) = char( ic+64 )
               END DO
            END IF
!
         ELSE IF( iz.EQ.218 .OR. iz.EQ.250 ) THEN
!
!           Prime machines:  ASCII+128
!
            IF( ic.GE.225 .AND. ic.LE.250 ) THEN
               subnam( 1: 1 ) = char( ic-32 )
               DO i = 2, 6
                  ic = ichar( subnam( i: i ) )
                  IF( ic.GE.225 .AND. ic.LE.250 ) &
     &               subnam( i: i ) = char( ic-32 )
               END DO
            END IF
         END IF
!
         IF( subnam( 2:6 ).EQ.'GGHRD' .OR. &
     &       subnam( 2:6 ).EQ.'GGHD3' ) THEN
            iparmq = 1
            IF( nh.GE.k22min ) &
     &         iparmq = 2
         ELSE IF ( subnam( 4:6 ).EQ.'EXC' ) THEN
            IF( nh.GE.kacmin ) &
     &         iparmq = 1
            IF( nh.GE.k22min ) &
     &         iparmq = 2
         ELSE IF ( subnam( 2:6 ).EQ.'HSEQR' .OR. &
     &             subnam( 2:5 ).EQ.'LAQR' ) THEN
            IF( ns.GE.kacmin ) &
     &         iparmq = 1
            IF( ns.GE.k22min ) &
     &         iparmq = 2
         END IF
!
      ELSE IF( ispec.EQ.icost ) THEN
!
!        === Relative cost of near-the-diagonal chase vs
!            BLAS updates ===
!
         iparmq = rcost
      ELSE
!        ===== invalid value of ispec =====
         iparmq = -1
!
      END IF
!
!     ==== End of IPARMQ ====
!

      END
!> \brief \b IEEECK
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download IEEECK + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ieeeck.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ieeeck.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ieeeck.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       INTEGER          FUNCTION IEEECK( ISPEC, ZERO, ONE )
!
!       .. Scalar Arguments ..
!       INTEGER            ISPEC
!       REAL               ONE, ZERO
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> IEEECK is called from the ILAENV to verify that Infinity and
!> possibly NaN arithmetic is safe (i.e. will not trap).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ISPEC
!> \verbatim
!>          ISPEC is INTEGER
!>          Specifies whether to test just for infinity arithmetic
!>          or whether to test for infinity and NaN arithmetic.
!>          = 0: Verify infinity arithmetic only.
!>          = 1: Verify infinity and NaN arithmetic.
!> \endverbatim
!>
!> \param[in] ZERO
!> \verbatim
!>          ZERO is REAL
!>          Must contain the value 0.0
!>          This is passed to prevent the compiler from optimizing
!>          away this code.
!> \endverbatim
!>
!> \param[in] ONE
!> \verbatim
!>          ONE is REAL
!>          Must contain the value 1.0
!>          This is passed to prevent the compiler from optimizing
!>          away this code.
!>
!>  RETURN VALUE:  INTEGER
!>          = 0:  Arithmetic failed to produce the correct answers
!>          = 1:  Arithmetic produced the correct answers
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup ieeeck
!
!  =====================================================================

      INTEGER          FUNCTION ieeeck( ISPEC, ZERO, ONE )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            ispec
      REAL               one, zero
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      REAL               nan1, nan2, nan3, nan4, nan5, nan6, neginf, &
     &                   negzro, newzro, posinf
!     ..
!     .. Executable Statements ..
      ieeeck = 1
!
      posinf = one / zero
      IF( posinf.LE.one ) THEN
         ieeeck = 0
         RETURN
      END IF
!
      neginf = -one / zero
      IF( neginf.GE.zero ) THEN
         ieeeck = 0
         RETURN
      END IF
!
      negzro = one / ( neginf+one )
      IF( negzro.NE.zero ) THEN
         ieeeck = 0
         RETURN
      END IF
!
      neginf = one / negzro
      IF( neginf.GE.zero ) THEN
         ieeeck = 0
         RETURN
      END IF
!
      newzro = negzro + zero
      IF( newzro.NE.zero ) THEN
         ieeeck = 0
         RETURN
      END IF
!
      posinf = one / newzro
      IF( posinf.LE.one ) THEN
         ieeeck = 0
         RETURN
      END IF
!
      neginf = neginf*posinf
      IF( neginf.GE.zero ) THEN
         ieeeck = 0
         RETURN
      END IF
!
      posinf = posinf*posinf
      IF( posinf.LE.one ) THEN
         ieeeck = 0
         RETURN
      END IF
!
!
!
!
!     Return if we were only asked to check infinity arithmetic
!
      IF( ispec.EQ.0 ) &
     &   RETURN
!
      nan1 = posinf + neginf
!
      nan2 = posinf / neginf
!
      nan3 = posinf / posinf
!
      nan4 = posinf*zero
!
      nan5 = neginf*negzro
!
      nan6 = nan5*zero
!
      IF( nan1.EQ.nan1 ) THEN
         ieeeck = 0
         RETURN
      END IF
!
      IF( nan2.EQ.nan2 ) THEN
         ieeeck = 0
         RETURN
      END IF
!
      IF( nan3.EQ.nan3 ) THEN
         ieeeck = 0
         RETURN
      END IF
!
      IF( nan4.EQ.nan4 ) THEN
         ieeeck = 0
         RETURN
      END IF
!
      IF( nan5.EQ.nan5 ) THEN
         ieeeck = 0
         RETURN
      END IF
!
      IF( nan6.EQ.nan6 ) THEN
         ieeeck = 0
         RETURN
      END IF
!
      RETURN

      END
!> \brief \b IDAMAX
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       INTEGER FUNCTION IDAMAX(N,DX,INCX)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION DX(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    IDAMAX finds the index of the first element having maximum absolute value.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         number of elements in input vector(s)
!> \endverbatim
!>
!> \param[in] DX
!> \verbatim
!>          DX is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>         storage spacing between elements of DX
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup iamax
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, linpack, 3/11/78.
!>     modified 3/93 to return if incx .le. 0.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================

      INTEGER FUNCTION idamax(N,DX,INCX)
!
!  -- Reference BLAS level1 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER incx,n
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION dx(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      DOUBLE PRECISION dmax
      INTEGER i,ix
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC dabs
!     ..
      idamax = 0
      IF (n.LT.1 .OR. incx.LE.0) RETURN
      idamax = 1
      IF (n.EQ.1) RETURN
      IF (incx.EQ.1) THEN
!
!        code for increment equal to 1
!
         dmax = dabs(dx(1))
         DO i = 2,n
            IF (dabs(dx(i)).GT.dmax) THEN
               idamax = i
               dmax = dabs(dx(i))
            END IF
         END DO
      ELSE
!
!        code for increment not equal to 1
!
         ix = 1
         dmax = dabs(dx(1))
         ix = ix + incx
         DO i = 2,n
            IF (dabs(dx(ix)).GT.dmax) THEN
               idamax = i
               dmax = dabs(dx(ix))
            END IF
            ix = ix + incx
         END DO
      END IF
      RETURN
!
!     End of IDAMAX
!

      END
!> \brief \b DTREVC3
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download DTREVC3 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtrevc3.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtrevc3.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtrevc3.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       SUBROUTINE DTREVC3( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL,
!                           VR, LDVR, MM, M, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          HOWMNY, SIDE
!       INTEGER            INFO, LDT, LDVL, LDVR, LWORK, M, MM, N
!       ..
!       .. Array Arguments ..
!       LOGICAL            SELECT( * )
!       DOUBLE PRECISION   T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ),
!      $                   WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DTREVC3 computes some or all of the right and/or left eigenvectors of
!> a real upper quasi-triangular matrix T.
!> Matrices of this type are produced by the Schur factorization of
!> a real general matrix:  A = Q*T*Q**T, as computed by DHSEQR.
!>
!> The right eigenvector x and the left eigenvector y of T corresponding
!> to an eigenvalue w are defined by:
!>
!>    T*x = w*x,     (y**T)*T = w*(y**T)
!>
!> where y**T denotes the transpose of the vector y.
!> The eigenvalues are not input to this routine, but are read directly
!> from the diagonal blocks of T.
!>
!> This routine returns the matrices X and/or Y of right and left
!> eigenvectors of T, or the products Q*X and/or Q*Y, where Q is an
!> input matrix. If Q is the orthogonal factor that reduces a matrix
!> A to Schur form T, then Q*X and Q*Y are the matrices of right and
!> left eigenvectors of A.
!>
!> This uses a Level 3 BLAS version of the back transformation.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'R':  compute right eigenvectors only;
!>          = 'L':  compute left eigenvectors only;
!>          = 'B':  compute both right and left eigenvectors.
!> \endverbatim
!>
!> \param[in] HOWMNY
!> \verbatim
!>          HOWMNY is CHARACTER*1
!>          = 'A':  compute all right and/or left eigenvectors;
!>          = 'B':  compute all right and/or left eigenvectors,
!>                  backtransformed by the matrices in VR and/or VL;
!>          = 'S':  compute selected right and/or left eigenvectors,
!>                  as indicated by the logical array SELECT.
!> \endverbatim
!>
!> \param[in,out] SELECT
!> \verbatim
!>          SELECT is LOGICAL array, dimension (N)
!>          If HOWMNY = 'S', SELECT specifies the eigenvectors to be
!>          computed.
!>          If w(j) is a real eigenvalue, the corresponding real
!>          eigenvector is computed if SELECT(j) is .TRUE..
!>          If w(j) and w(j+1) are the real and imaginary parts of a
!>          complex eigenvalue, the corresponding complex eigenvector is
!>          computed if either SELECT(j) or SELECT(j+1) is .TRUE., and
!>          on exit SELECT(j) is set to .TRUE. and SELECT(j+1) is set to
!>          .FALSE..
!>          Not referenced if HOWMNY = 'A' or 'B'.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix T. N >= 0.
!> \endverbatim
!>
!> \param[in] T
!> \verbatim
!>          T is DOUBLE PRECISION array, dimension (LDT,N)
!>          The upper quasi-triangular matrix T in Schur canonical form.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of the array T. LDT >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] VL
!> \verbatim
!>          VL is DOUBLE PRECISION array, dimension (LDVL,MM)
!>          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must
!>          contain an N-by-N matrix Q (usually the orthogonal matrix Q
!>          of Schur vectors returned by DHSEQR).
!>          On exit, if SIDE = 'L' or 'B', VL contains:
!>          if HOWMNY = 'A', the matrix Y of left eigenvectors of T;
!>          if HOWMNY = 'B', the matrix Q*Y;
!>          if HOWMNY = 'S', the left eigenvectors of T specified by
!>                           SELECT, stored consecutively in the columns
!>                           of VL, in the same order as their
!>                           eigenvalues.
!>          A complex eigenvector corresponding to a complex eigenvalue
!>          is stored in two consecutive columns, the first holding the
!>          real part, and the second the imaginary part.
!>          Not referenced if SIDE = 'R'.
!> \endverbatim
!>
!> \param[in] LDVL
!> \verbatim
!>          LDVL is INTEGER
!>          The leading dimension of the array VL.
!>          LDVL >= 1, and if SIDE = 'L' or 'B', LDVL >= N.
!> \endverbatim
!>
!> \param[in,out] VR
!> \verbatim
!>          VR is DOUBLE PRECISION array, dimension (LDVR,MM)
!>          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must
!>          contain an N-by-N matrix Q (usually the orthogonal matrix Q
!>          of Schur vectors returned by DHSEQR).
!>          On exit, if SIDE = 'R' or 'B', VR contains:
!>          if HOWMNY = 'A', the matrix X of right eigenvectors of T;
!>          if HOWMNY = 'B', the matrix Q*X;
!>          if HOWMNY = 'S', the right eigenvectors of T specified by
!>                           SELECT, stored consecutively in the columns
!>                           of VR, in the same order as their
!>                           eigenvalues.
!>          A complex eigenvector corresponding to a complex eigenvalue
!>          is stored in two consecutive columns, the first holding the
!>          real part and the second the imaginary part.
!>          Not referenced if SIDE = 'L'.
!> \endverbatim
!>
!> \param[in] LDVR
!> \verbatim
!>          LDVR is INTEGER
!>          The leading dimension of the array VR.
!>          LDVR >= 1, and if SIDE = 'R' or 'B', LDVR >= N.
!> \endverbatim
!>
!> \param[in] MM
!> \verbatim
!>          MM is INTEGER
!>          The number of columns in the arrays VL and/or VR. MM >= M.
!> \endverbatim
!>
!> \param[out] M
!> \verbatim
!>          M is INTEGER
!>          The number of columns in the arrays VL and/or VR actually
!>          used to store the eigenvectors.
!>          If HOWMNY = 'A' or 'B', M is set to N.
!>          Each selected real eigenvector occupies one column and each
!>          selected complex eigenvector occupies two columns.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of array WORK. LWORK >= max(1,3*N).
!>          For optimum performance, LWORK >= N + 2*N*NB, where NB is
!>          the optimal blocksize.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup trevc3
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The algorithm used in this program is basically backward (forward)
!>  substitution, with scaling to make the the code robust against
!>  possible overflow.
!>
!>  Each eigenvector is normalized so that the element of largest
!>  magnitude has magnitude 1; here the magnitude of a complex number
!>  (x,y) is taken to be |x| + |y|.
!> \endverbatim
!>
!  =====================================================================

      SUBROUTINE dtrevc3( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, &
     &                    VR, LDVR, MM, M, WORK, LWORK, INFO )
      IMPLICIT NONE
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          HOWMNY, SIDE
      INTEGER            INFO, LDT, LDVL, LDVR, LWORK, M, MM, N
!     ..
!     .. Array Arguments ..
      LOGICAL            SELECT( * )
      DOUBLE PRECISION   T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ), &
     &                   work( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      parameter( zero = 0.0d+0, one = 1.0d+0 )
      INTEGER            NBMIN, NBMAX
      parameter( nbmin = 8, nbmax = 128 )
!     ..
!     .. Local Scalars ..
      LOGICAL            ALLV, BOTHV, LEFTV, LQUERY, OVER, PAIR, &
     &                   rightv, somev
      INTEGER            I, IERR, II, IP, IS, J, J1, J2, JNXT, K, KI, &
     &                   iv, maxwrk, nb, ki2
      DOUBLE PRECISION   BETA, BIGNUM, EMAX, OVFL, REC, REMAX, SCALE, &
     &                   smin, smlnum, ulp, unfl, vcrit, vmax, wi, wr, &
     &                   xnorm
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      INTEGER            IDAMAX, ILAENV
!      DOUBLE PRECISION   DDOT, DLAMCH
!      EXTERNAL           lsame, idamax, ilaenv, ddot, dlamch
!     ..
!     .. External Subroutines ..
!      EXTERNAL           daxpy, dcopy, dgemv, dlaln2, dscal, &
!     &                   xerbla, &
!     &                   dgemm, dlaset, dlacpy
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          abs, max, sqrt
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION   X( 2, 2 )
      INTEGER            ISCOMPLEX( NBMAX )
!     ..
!     .. Executable Statements ..
!
!     Decode and test the input parameters
!
      bothv  = lsame( side, 'B' )
      rightv = lsame( side, 'R' ) .OR. bothv
      leftv  = lsame( side, 'L' ) .OR. bothv
!
      allv  = lsame( howmny, 'A' )
      over  = lsame( howmny, 'B' )
      somev = lsame( howmny, 'S' )
!
      info = 0
      nb = ilaenv( 1, 'DTREVC', side // howmny, n, -1, -1, -1 )
      maxwrk = max( 1, n + 2*n*nb )
      work(1) = maxwrk
      lquery = ( lwork.EQ.-1 )
      IF( .NOT.rightv .AND. .NOT.leftv ) THEN
         info = -1
      ELSE IF( .NOT.allv .AND. .NOT.over .AND. .NOT.somev ) THEN
         info = -2
      ELSE IF( n.LT.0 ) THEN
         info = -4
      ELSE IF( ldt.LT.max( 1, n ) ) THEN
         info = -6
      ELSE IF( ldvl.LT.1 .OR. ( leftv .AND. ldvl.LT.n ) ) THEN
         info = -8
      ELSE IF( ldvr.LT.1 .OR. ( rightv .AND. ldvr.LT.n ) ) THEN
         info = -10
      ELSE IF( lwork.LT.max( 1, 3*n ) .AND. .NOT.lquery ) THEN
         info = -14
      ELSE
!
!        Set M to the number of columns required to store the selected
!        eigenvectors, standardize the array SELECT if necessary, and
!        test MM.
!
         IF( somev ) THEN
            m = 0
            pair = .false.
            DO 10 j = 1, n
               IF( pair ) THEN
                  pair = .false.
                  SELECT( j ) = .false.
               ELSE
                  IF( j.LT.n ) THEN
                     IF( t( j+1, j ).EQ.zero ) THEN
                        IF( SELECT( j ) ) &
     &                     m = m + 1
                     ELSE
                        pair = .true.
                        IF( SELECT( j ) .OR. SELECT( j+1 ) ) THEN
                           SELECT( j ) = .true.
                           m = m + 2
                        END IF
                     END IF
                  ELSE
                     IF( SELECT( n ) ) &
     &                  m = m + 1
                  END IF
               END IF
10          CONTINUE
         ELSE
            m = n
         END IF
!
         IF( mm.LT.m ) THEN
            info = -11
         END IF
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DTREVC3', -info )
         RETURN
      ELSE IF( lquery ) THEN
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( n.EQ.0 ) &
     &   RETURN
!
!     Use blocked version of back-transformation if sufficient workspace.
!     Zero-out the workspace to avoid potential NaN propagation.
!
      IF( over .AND. lwork .GE. n + 2*n*nbmin ) THEN
         nb = (lwork - n) / (2*n)
         nb = min( nb, nbmax )
         CALL dlaset( 'F', n, 1+2*nb, zero, zero, work, n )
      ELSE
         nb = 1
      END IF
!
!     Set the constants to control overflow.
!
      unfl = dlamch( 'Safe minimum' )
      ovfl = one / unfl
      ulp = dlamch( 'Precision' )
      smlnum = unfl*( n / ulp )
      bignum = ( one-ulp ) / smlnum
!
!     Compute 1-norm of each column of strictly upper triangular
!     part of T to control overflow in triangular solver.
!
      work( 1 ) = zero
      DO 30 j = 2, n
         work( j ) = zero
         DO 20 i = 1, j - 1
            work( j ) = work( j ) + abs( t( i, j ) )
20       CONTINUE
30    CONTINUE
!
!     Index IP is used to specify the real or complex eigenvalue:
!       IP = 0, real eigenvalue,
!            1, first  of conjugate complex pair: (wr,wi)
!           -1, second of conjugate complex pair: (wr,wi)
!       ISCOMPLEX array stores IP for each column in current block.
!
      IF( rightv ) THEN
!
!        ============================================================
!        Compute right eigenvectors.
!
!        IV is index of column in current block.
!        For complex right vector, uses IV-1 for real part and IV for complex part.
!        Non-blocked version always uses IV=2;
!        blocked     version starts with IV=NB, goes down to 1 or 2.
!        (Note the "0-th" column is used for 1-norms computed above.)
         iv = 2
         IF( nb.GT.2 ) THEN
            iv = nb
         END IF

         ip = 0
         is = m
         DO 140 ki = n, 1, -1
            IF( ip.EQ.-1 ) THEN
!              previous iteration (ki+1) was second of conjugate pair,
!              so this ki is first of conjugate pair; skip to end of loop
               ip = 1
               GO TO 140
            ELSE IF( ki.EQ.1 ) THEN
!              last column, so this ki must be real eigenvalue
               ip = 0
            ELSE IF( t( ki, ki-1 ).EQ.zero ) THEN
!              zero on sub-diagonal, so this ki is real eigenvalue
               ip = 0
            ELSE
!              non-zero on sub-diagonal, so this ki is second of conjugate pair
               ip = -1
            END IF

            IF( somev ) THEN
               IF( ip.EQ.0 ) THEN
                  IF( .NOT.SELECT( ki ) ) &
     &               GO TO 140
               ELSE
                  IF( .NOT.SELECT( ki-1 ) ) &
     &               GO TO 140
               END IF
            END IF
!
!           Compute the KI-th eigenvalue (WR,WI).
!
            wr = t( ki, ki )
            wi = zero
            IF( ip.NE.0 ) &
     &         wi = sqrt( abs( t( ki, ki-1 ) ) )* &
     &              sqrt( abs( t( ki-1, ki ) ) )
            smin = max( ulp*( abs( wr )+abs( wi ) ), smlnum )
!
            IF( ip.EQ.0 ) THEN
!
!              --------------------------------------------------------
!              Real right eigenvector
!
               work( ki + iv*n ) = one
!
!              Form right-hand side.
!
               DO 50 k = 1, ki - 1
                  work( k + iv*n ) = -t( k, ki )
50             CONTINUE
!
!              Solve upper quasi-triangular system:
!              [ T(1:KI-1,1:KI-1) - WR ]*X = SCALE*WORK.
!
               jnxt = ki - 1
               DO 60 j = ki - 1, 1, -1
                  IF( j.GT.jnxt ) &
     &               GO TO 60
                  j1 = j
                  j2 = j
                  jnxt = j - 1
                  IF( j.GT.1 ) THEN
                     IF( t( j, j-1 ).NE.zero ) THEN
                        j1   = j - 1
                        jnxt = j - 2
                     END IF
                  END IF
!
                  IF( j1.EQ.j2 ) THEN
!
!                    1-by-1 diagonal block
!
                     CALL dlaln2( .false., 1, 1, smin, one, t( j, &
     &                            j ), &
     &                            ldt, one, one, work( j+iv*n ), n, wr, &
     &                            zero, x, 2, scale, xnorm, ierr )
!
!                    Scale X(1,1) to avoid overflow when updating
!                    the right-hand side.
!
                     IF( xnorm.GT.one ) THEN
                        IF( work( j ).GT.bignum / xnorm ) THEN
                           x( 1, 1 ) = x( 1, 1 ) / xnorm
                           scale = scale / xnorm
                        END IF
                     END IF
!
!                    Scale if necessary
!
                     IF( scale.NE.one ) &
     &                  CALL dscal( ki, scale, work( 1+iv*n ), 1 )
                     work( j+iv*n ) = x( 1, 1 )
!
!                    Update right-hand side
!
                     CALL daxpy( j-1, -x( 1, 1 ), t( 1, j ), 1, &
     &                           work( 1+iv*n ), 1 )
!
                  ELSE
!
!                    2-by-2 diagonal block
!
                     CALL dlaln2( .false., 2, 1, smin, one, &
     &                            t( j-1, j-1 ), ldt, one, one, &
     &                            work( j-1+iv*n ), n, wr, zero, x, 2, &
     &                            scale, xnorm, ierr )
!
!                    Scale X(1,1) and X(2,1) to avoid overflow when
!                    updating the right-hand side.
!
                     IF( xnorm.GT.one ) THEN
                        beta = max( work( j-1 ), work( j ) )
                        IF( beta.GT.bignum / xnorm ) THEN
                           x( 1, 1 ) = x( 1, 1 ) / xnorm
                           x( 2, 1 ) = x( 2, 1 ) / xnorm
                           scale = scale / xnorm
                        END IF
                     END IF
!
!                    Scale if necessary
!
                     IF( scale.NE.one ) &
     &                  CALL dscal( ki, scale, work( 1+iv*n ), 1 )
                     work( j-1+iv*n ) = x( 1, 1 )
                     work( j  +iv*n ) = x( 2, 1 )
!
!                    Update right-hand side
!
                     CALL daxpy( j-2, -x( 1, 1 ), t( 1, j-1 ), 1, &
     &                           work( 1+iv*n ), 1 )
                     CALL daxpy( j-2, -x( 2, 1 ), t( 1, j ), 1, &
     &                           work( 1+iv*n ), 1 )
                  END IF
60             CONTINUE
!
!              Copy the vector x or Q*x to VR and normalize.
!
               IF( .NOT.over ) THEN
!                 ------------------------------
!                 no back-transform: copy x to VR and normalize.
                  CALL dcopy( ki, work( 1 + iv*n ), 1, vr( 1, is ), &
     &                        1 )
!
                  ii = idamax( ki, vr( 1, is ), 1 )
                  remax = one / abs( vr( ii, is ) )
                  CALL dscal( ki, remax, vr( 1, is ), 1 )
!
                  DO 70 k = ki + 1, n
                     vr( k, is ) = zero
70                CONTINUE
!
               ELSE IF( nb.EQ.1 ) THEN
!                 ------------------------------
!                 version 1: back-transform each vector with GEMV, Q*x.
                  IF( ki.GT.1 ) &
     &               CALL dgemv( 'N', n, ki-1, one, vr, ldvr, &
     &                           work( 1 + iv*n ), 1, work( ki + iv*n ), &
     &                           vr( 1, ki ), 1 )
!
                  ii = idamax( n, vr( 1, ki ), 1 )
                  remax = one / abs( vr( ii, ki ) )
                  CALL dscal( n, remax, vr( 1, ki ), 1 )
!
               ELSE
!                 ------------------------------
!                 version 2: back-transform block of vectors with GEMM
!                 zero out below vector
                  DO k = ki + 1, n
                     work( k + iv*n ) = zero
                  END DO
                  iscomplex( iv ) = ip
!                 back-transform and normalization is done below
               END IF
            ELSE
!
!              --------------------------------------------------------
!              Complex right eigenvector.
!
!              Initial solve
!              [ ( T(KI-1,KI-1) T(KI-1,KI) ) - (WR + I*WI) ]*X = 0.
!              [ ( T(KI,  KI-1) T(KI,  KI) )               ]
!
               IF( abs( t( ki-1, ki ) ).GE.abs( t( ki, ki-1 ) ) ) THEN
                  work( ki-1 + (iv-1)*n ) = one
                  work( ki   + (iv  )*n ) = wi / t( ki-1, ki )
               ELSE
                  work( ki-1 + (iv-1)*n ) = -wi / t( ki, ki-1 )
                  work( ki   + (iv  )*n ) = one
               END IF
               work( ki   + (iv-1)*n ) = zero
               work( ki-1 + (iv  )*n ) = zero
!
!              Form right-hand side.
!
               DO 80 k = 1, ki - 2
                  work( k+(iv-1)*n ) = -work( ki-1+(iv-1)*n )*t(k,ki-1)
                  work( k+(iv  )*n ) = -work( ki  +(iv  )*n )*t(k,ki  )
80             CONTINUE
!
!              Solve upper quasi-triangular system:
!              [ T(1:KI-2,1:KI-2) - (WR+i*WI) ]*X = SCALE*(WORK+i*WORK2)
!
               jnxt = ki - 2
               DO 90 j = ki - 2, 1, -1
                  IF( j.GT.jnxt ) &
     &               GO TO 90
                  j1 = j
                  j2 = j
                  jnxt = j - 1
                  IF( j.GT.1 ) THEN
                     IF( t( j, j-1 ).NE.zero ) THEN
                        j1   = j - 1
                        jnxt = j - 2
                     END IF
                  END IF
!
                  IF( j1.EQ.j2 ) THEN
!
!                    1-by-1 diagonal block
!
                     CALL dlaln2( .false., 1, 2, smin, one, t( j, &
     &                            j ), &
     &                            ldt, one, one, work( j+(iv-1)*n ), n, &
     &                            wr, wi, x, 2, scale, xnorm, ierr )
!
!                    Scale X(1,1) and X(1,2) to avoid overflow when
!                    updating the right-hand side.
!
                     IF( xnorm.GT.one ) THEN
                        IF( work( j ).GT.bignum / xnorm ) THEN
                           x( 1, 1 ) = x( 1, 1 ) / xnorm
                           x( 1, 2 ) = x( 1, 2 ) / xnorm
                           scale = scale / xnorm
                        END IF
                     END IF
!
!                    Scale if necessary
!
                     IF( scale.NE.one ) THEN
                        CALL dscal( ki, scale, work( 1+(iv-1)*n ), &
     &                              1 )
                        CALL dscal( ki, scale, work( 1+(iv  )*n ), &
     &                              1 )
                     END IF
                     work( j+(iv-1)*n ) = x( 1, 1 )
                     work( j+(iv  )*n ) = x( 1, 2 )
!
!                    Update the right-hand side
!
                     CALL daxpy( j-1, -x( 1, 1 ), t( 1, j ), 1, &
     &                           work( 1+(iv-1)*n ), 1 )
                     CALL daxpy( j-1, -x( 1, 2 ), t( 1, j ), 1, &
     &                           work( 1+(iv  )*n ), 1 )
!
                  ELSE
!
!                    2-by-2 diagonal block
!
                     CALL dlaln2( .false., 2, 2, smin, one, &
     &                            t( j-1, j-1 ), ldt, one, one, &
     &                            work( j-1+(iv-1)*n ), n, wr, wi, x, 2, &
     &                            scale, xnorm, ierr )
!
!                    Scale X to avoid overflow when updating
!                    the right-hand side.
!
                     IF( xnorm.GT.one ) THEN
                        beta = max( work( j-1 ), work( j ) )
                        IF( beta.GT.bignum / xnorm ) THEN
                           rec = one / xnorm
                           x( 1, 1 ) = x( 1, 1 )*rec
                           x( 1, 2 ) = x( 1, 2 )*rec
                           x( 2, 1 ) = x( 2, 1 )*rec
                           x( 2, 2 ) = x( 2, 2 )*rec
                           scale = scale*rec
                        END IF
                     END IF
!
!                    Scale if necessary
!
                     IF( scale.NE.one ) THEN
                        CALL dscal( ki, scale, work( 1+(iv-1)*n ), &
     &                              1 )
                        CALL dscal( ki, scale, work( 1+(iv  )*n ), &
     &                              1 )
                     END IF
                     work( j-1+(iv-1)*n ) = x( 1, 1 )
                     work( j  +(iv-1)*n ) = x( 2, 1 )
                     work( j-1+(iv  )*n ) = x( 1, 2 )
                     work( j  +(iv  )*n ) = x( 2, 2 )
!
!                    Update the right-hand side
!
                     CALL daxpy( j-2, -x( 1, 1 ), t( 1, j-1 ), 1, &
     &                           work( 1+(iv-1)*n   ), 1 )
                     CALL daxpy( j-2, -x( 2, 1 ), t( 1, j ), 1, &
     &                           work( 1+(iv-1)*n   ), 1 )
                     CALL daxpy( j-2, -x( 1, 2 ), t( 1, j-1 ), 1, &
     &                           work( 1+(iv  )*n ), 1 )
                     CALL daxpy( j-2, -x( 2, 2 ), t( 1, j ), 1, &
     &                           work( 1+(iv  )*n ), 1 )
                  END IF
90             CONTINUE
!
!              Copy the vector x or Q*x to VR and normalize.
!
               IF( .NOT.over ) THEN
!                 ------------------------------
!                 no back-transform: copy x to VR and normalize.
                  CALL dcopy( ki, work( 1+(iv-1)*n ), 1, vr(1,is-1), &
     &                        1 )
                  CALL dcopy( ki, work( 1+(iv  )*n ), 1, vr(1,is  ), &
     &                        1 )
!
                  emax = zero
                  DO 100 k = 1, ki
                     emax = max( emax, abs( vr( k, is-1 ) )+ &
     &                                 abs( vr( k, is   ) ) )
100               CONTINUE
                  remax = one / emax
                  CALL dscal( ki, remax, vr( 1, is-1 ), 1 )
                  CALL dscal( ki, remax, vr( 1, is   ), 1 )
!
                  DO 110 k = ki + 1, n
                     vr( k, is-1 ) = zero
                     vr( k, is   ) = zero
110               CONTINUE
!
               ELSE IF( nb.EQ.1 ) THEN
!                 ------------------------------
!                 version 1: back-transform each vector with GEMV, Q*x.
                  IF( ki.GT.2 ) THEN
                     CALL dgemv( 'N', n, ki-2, one, vr, ldvr, &
     &                           work( 1    + (iv-1)*n ), 1, &
     &                           work( ki-1 + (iv-1)*n ), vr(1,ki-1), 1)
                     CALL dgemv( 'N', n, ki-2, one, vr, ldvr, &
     &                           work( 1  + (iv)*n ), 1, &
     &                           work( ki + (iv)*n ), vr( 1, ki ), 1 )
                  ELSE
                     CALL dscal( n, work(ki-1+(iv-1)*n), vr(1,ki-1), &
     &                           1)
                     CALL dscal( n, work(ki  +(iv  )*n), vr(1,ki  ), &
     &                           1)
                  END IF
!
                  emax = zero
                  DO 120 k = 1, n
                     emax = max( emax, abs( vr( k, ki-1 ) )+ &
     &                                 abs( vr( k, ki   ) ) )
120               CONTINUE
                  remax = one / emax
                  CALL dscal( n, remax, vr( 1, ki-1 ), 1 )
                  CALL dscal( n, remax, vr( 1, ki   ), 1 )
!
               ELSE
!                 ------------------------------
!                 version 2: back-transform block of vectors with GEMM
!                 zero out below vector
                  DO k = ki + 1, n
                     work( k + (iv-1)*n ) = zero
                     work( k + (iv  )*n ) = zero
                  END DO
                  iscomplex( iv-1 ) = -ip
                  iscomplex( iv   ) =  ip
                  iv = iv - 1
!                 back-transform and normalization is done below
               END IF
            END IF

            IF( nb.GT.1 ) THEN
!              --------------------------------------------------------
!              Blocked version of back-transform
!              For complex case, KI2 includes both vectors (KI-1 and KI)
               IF( ip.EQ.0 ) THEN
                  ki2 = ki
               ELSE
                  ki2 = ki - 1
               END IF

!              Columns IV:NB of work are valid vectors.
!              When the number of vectors stored reaches NB-1 or NB,
!              or if this was last vector, do the GEMM
               IF( (iv.LE.2) .OR. (ki2.EQ.1) ) THEN
                  CALL dgemm( 'N', 'N', n, nb-iv+1, ki2+nb-iv, one, &
     &                        vr, ldvr, &
     &                        work( 1 + (iv)*n    ), n, &
     &                        zero, &
     &                        work( 1 + (nb+iv)*n ), n )
!                 normalize vectors
                  DO k = iv, nb
                     IF( iscomplex(k).EQ.0 ) THEN
!                       real eigenvector
                        ii = idamax( n, work( 1 + (nb+k)*n ), 1 )
                        remax = one / abs( work( ii + (nb+k)*n ) )
                     ELSE IF( iscomplex(k).EQ.1 ) THEN
!                       first eigenvector of conjugate pair
                        emax = zero
                        DO ii = 1, n
                           emax = max( emax, &
     &                                 abs( work( ii + (nb+k  )*n ) )+ &
     &                                 abs( work( ii + (nb+k+1)*n ) ) )
                        END DO
                        remax = one / emax
!                    else if ISCOMPLEX(K).EQ.-1
!                       second eigenvector of conjugate pair
!                       reuse same REMAX as previous K
                     END IF
                     CALL dscal( n, remax, work( 1 + (nb+k)*n ), 1 )
                  END DO
                  CALL dlacpy( 'F', n, nb-iv+1, &
     &                         work( 1 + (nb+iv)*n ), n, &
     &                         vr( 1, ki2 ), ldvr )
                  iv = nb
               ELSE
                  iv = iv - 1
               END IF
            END IF ! blocked back-transform
!
            is = is - 1
            IF( ip.NE.0 ) &
     &         is = is - 1
140      CONTINUE
      END IF

      IF( leftv ) THEN
!
!        ============================================================
!        Compute left eigenvectors.
!
!        IV is index of column in current block.
!        For complex left vector, uses IV for real part and IV+1 for complex part.
!        Non-blocked version always uses IV=1;
!        blocked     version starts with IV=1, goes up to NB-1 or NB.
!        (Note the "0-th" column is used for 1-norms computed above.)
         iv = 1
         ip = 0
         is = 1
         DO 260 ki = 1, n
            IF( ip.EQ.1 ) THEN
!              previous iteration (ki-1) was first of conjugate pair,
!              so this ki is second of conjugate pair; skip to end of loop
               ip = -1
               GO TO 260
            ELSE IF( ki.EQ.n ) THEN
!              last column, so this ki must be real eigenvalue
               ip = 0
            ELSE IF( t( ki+1, ki ).EQ.zero ) THEN
!              zero on sub-diagonal, so this ki is real eigenvalue
               ip = 0
            ELSE
!              non-zero on sub-diagonal, so this ki is first of conjugate pair
               ip = 1
            END IF
!
            IF( somev ) THEN
               IF( .NOT.SELECT( ki ) ) &
     &            GO TO 260
            END IF
!
!           Compute the KI-th eigenvalue (WR,WI).
!
            wr = t( ki, ki )
            wi = zero
            IF( ip.NE.0 ) &
     &         wi = sqrt( abs( t( ki, ki+1 ) ) )* &
     &              sqrt( abs( t( ki+1, ki ) ) )
            smin = max( ulp*( abs( wr )+abs( wi ) ), smlnum )
!
            IF( ip.EQ.0 ) THEN
!
!              --------------------------------------------------------
!              Real left eigenvector
!
               work( ki + iv*n ) = one
!
!              Form right-hand side.
!
               DO 160 k = ki + 1, n
                  work( k + iv*n ) = -t( ki, k )
160            CONTINUE
!
!              Solve transposed quasi-triangular system:
!              [ T(KI+1:N,KI+1:N) - WR ]**T * X = SCALE*WORK
!
               vmax = one
               vcrit = bignum
!
               jnxt = ki + 1
               DO 170 j = ki + 1, n
                  IF( j.LT.jnxt ) &
     &               GO TO 170
                  j1 = j
                  j2 = j
                  jnxt = j + 1
                  IF( j.LT.n ) THEN
                     IF( t( j+1, j ).NE.zero ) THEN
                        j2 = j + 1
                        jnxt = j + 2
                     END IF
                  END IF
!
                  IF( j1.EQ.j2 ) THEN
!
!                    1-by-1 diagonal block
!
!                    Scale if necessary to avoid overflow when forming
!                    the right-hand side.
!
                     IF( work( j ).GT.vcrit ) THEN
                        rec = one / vmax
                        CALL dscal( n-ki+1, rec, work( ki+iv*n ), 1 )
                        vmax = one
                        vcrit = bignum
                     END IF
!
                     work( j+iv*n ) = work( j+iv*n ) - &
     &                                ddot( j-ki-1, t( ki+1, j ), 1, &
     &                                      work( ki+1+iv*n ), 1 )
!
!                    Solve [ T(J,J) - WR ]**T * X = WORK
!
                     CALL dlaln2( .false., 1, 1, smin, one, t( j, &
     &                            j ), &
     &                            ldt, one, one, work( j+iv*n ), n, wr, &
     &                            zero, x, 2, scale, xnorm, ierr )
!
!                    Scale if necessary
!
                     IF( scale.NE.one ) &
     &                  CALL dscal( n-ki+1, scale, work( ki+iv*n ), &
     &                              1 )
                     work( j+iv*n ) = x( 1, 1 )
                     vmax = max( abs( work( j+iv*n ) ), vmax )
                     vcrit = bignum / vmax
!
                  ELSE
!
!                    2-by-2 diagonal block
!
!                    Scale if necessary to avoid overflow when forming
!                    the right-hand side.
!
                     beta = max( work( j ), work( j+1 ) )
                     IF( beta.GT.vcrit ) THEN
                        rec = one / vmax
                        CALL dscal( n-ki+1, rec, work( ki+iv*n ), 1 )
                        vmax = one
                        vcrit = bignum
                     END IF
!
                     work( j+iv*n ) = work( j+iv*n ) - &
     &                                ddot( j-ki-1, t( ki+1, j ), 1, &
     &                                      work( ki+1+iv*n ), 1 )
!
                     work( j+1+iv*n ) = work( j+1+iv*n ) - &
     &                                  ddot( j-ki-1, t( ki+1, j+1 ), &
     &                                        1, &
     &                                        work( ki+1+iv*n ), 1 )
!
!                    Solve
!                    [ T(J,J)-WR   T(J,J+1)      ]**T * X = SCALE*( WORK1 )
!                    [ T(J+1,J)    T(J+1,J+1)-WR ]                ( WORK2 )
!
                     CALL dlaln2( .true., 2, 1, smin, one, t( j, j ), &
     &                            ldt, one, one, work( j+iv*n ), n, wr, &
     &                            zero, x, 2, scale, xnorm, ierr )
!
!                    Scale if necessary
!
                     IF( scale.NE.one ) &
     &                  CALL dscal( n-ki+1, scale, work( ki+iv*n ), &
     &                              1 )
                     work( j  +iv*n ) = x( 1, 1 )
                     work( j+1+iv*n ) = x( 2, 1 )
!
                     vmax = max( abs( work( j  +iv*n ) ), &
     &                           abs( work( j+1+iv*n ) ), vmax )
                     vcrit = bignum / vmax
!
                  END IF
170            CONTINUE
!
!              Copy the vector x or Q*x to VL and normalize.
!
               IF( .NOT.over ) THEN
!                 ------------------------------
!                 no back-transform: copy x to VL and normalize.
                  CALL dcopy( n-ki+1, work( ki + iv*n ), 1, &
     &                                vl( ki, is ), 1 )
!
                  ii = idamax( n-ki+1, vl( ki, is ), 1 ) + ki - 1
                  remax = one / abs( vl( ii, is ) )
                  CALL dscal( n-ki+1, remax, vl( ki, is ), 1 )
!
                  DO 180 k = 1, ki - 1
                     vl( k, is ) = zero
180               CONTINUE
!
               ELSE IF( nb.EQ.1 ) THEN
!                 ------------------------------
!                 version 1: back-transform each vector with GEMV, Q*x.
                  IF( ki.LT.n ) &
     &               CALL dgemv( 'N', n, n-ki, one, &
     &                           vl( 1, ki+1 ), ldvl, &
     &                           work( ki+1 + iv*n ), 1, &
     &                           work( ki   + iv*n ), vl( 1, ki ), 1 )
!
                  ii = idamax( n, vl( 1, ki ), 1 )
                  remax = one / abs( vl( ii, ki ) )
                  CALL dscal( n, remax, vl( 1, ki ), 1 )
!
               ELSE
!                 ------------------------------
!                 version 2: back-transform block of vectors with GEMM
!                 zero out above vector
!                 could go from KI-NV+1 to KI-1
                  DO k = 1, ki - 1
                     work( k + iv*n ) = zero
                  END DO
                  iscomplex( iv ) = ip
!                 back-transform and normalization is done below
               END IF
            ELSE
!
!              --------------------------------------------------------
!              Complex left eigenvector.
!
!              Initial solve:
!              [ ( T(KI,KI)    T(KI,KI+1)  )**T - (WR - I* WI) ]*X = 0.
!              [ ( T(KI+1,KI) T(KI+1,KI+1) )                   ]
!
               IF( abs( t( ki, ki+1 ) ).GE.abs( t( ki+1, ki ) ) ) THEN
                  work( ki   + (iv  )*n ) = wi / t( ki, ki+1 )
                  work( ki+1 + (iv+1)*n ) = one
               ELSE
                  work( ki   + (iv  )*n ) = one
                  work( ki+1 + (iv+1)*n ) = -wi / t( ki+1, ki )
               END IF
               work( ki+1 + (iv  )*n ) = zero
               work( ki   + (iv+1)*n ) = zero
!
!              Form right-hand side.
!
               DO 190 k = ki + 2, n
                  work( k+(iv  )*n ) = -work( ki  +(iv  )*n )*t(ki,  k)
                  work( k+(iv+1)*n ) = -work( ki+1+(iv+1)*n )*t(ki+1,k)
190            CONTINUE
!
!              Solve transposed quasi-triangular system:
!              [ T(KI+2:N,KI+2:N)**T - (WR-i*WI) ]*X = WORK1+i*WORK2
!
               vmax = one
               vcrit = bignum
!
               jnxt = ki + 2
               DO 200 j = ki + 2, n
                  IF( j.LT.jnxt ) &
     &               GO TO 200
                  j1 = j
                  j2 = j
                  jnxt = j + 1
                  IF( j.LT.n ) THEN
                     IF( t( j+1, j ).NE.zero ) THEN
                        j2 = j + 1
                        jnxt = j + 2
                     END IF
                  END IF
!
                  IF( j1.EQ.j2 ) THEN
!
!                    1-by-1 diagonal block
!
!                    Scale if necessary to avoid overflow when
!                    forming the right-hand side elements.
!
                     IF( work( j ).GT.vcrit ) THEN
                        rec = one / vmax
                        CALL dscal( n-ki+1, rec, work(ki+(iv  )*n), &
     &                              1 )
                        CALL dscal( n-ki+1, rec, work(ki+(iv+1)*n), &
     &                              1 )
                        vmax = one
                        vcrit = bignum
                     END IF
!
                     work( j+(iv  )*n ) = work( j+(iv)*n ) - &
     &                                  ddot( j-ki-2, t( ki+2, j ), &
     &                                        1, &
     &                                        work( ki+2+(iv)*n ), 1 )
                     work( j+(iv+1)*n ) = work( j+(iv+1)*n ) - &
     &                                  ddot( j-ki-2, t( ki+2, j ), &
     &                                        1, &
     &                                        work( ki+2+(iv+1)*n ), 1 )
!
!                    Solve [ T(J,J)-(WR-i*WI) ]*(X11+i*X12)= WK+I*WK2
!
                     CALL dlaln2( .false., 1, 2, smin, one, t( j, &
     &                            j ), &
     &                            ldt, one, one, work( j+iv*n ), n, wr, &
     &                            -wi, x, 2, scale, xnorm, ierr )
!
!                    Scale if necessary
!
                     IF( scale.NE.one ) THEN
                        CALL dscal( n-ki+1, scale, work(ki+(iv  )*n), &
     &                              1)
                        CALL dscal( n-ki+1, scale, work(ki+(iv+1)*n), &
     &                              1)
                     END IF
                     work( j+(iv  )*n ) = x( 1, 1 )
                     work( j+(iv+1)*n ) = x( 1, 2 )
                     vmax = max( abs( work( j+(iv  )*n ) ), &
     &                           abs( work( j+(iv+1)*n ) ), vmax )
                     vcrit = bignum / vmax
!
                  ELSE
!
!                    2-by-2 diagonal block
!
!                    Scale if necessary to avoid overflow when forming
!                    the right-hand side elements.
!
                     beta = max( work( j ), work( j+1 ) )
                     IF( beta.GT.vcrit ) THEN
                        rec = one / vmax
                        CALL dscal( n-ki+1, rec, work(ki+(iv  )*n), &
     &                              1 )
                        CALL dscal( n-ki+1, rec, work(ki+(iv+1)*n), &
     &                              1 )
                        vmax = one
                        vcrit = bignum
                     END IF
!
                     work( j  +(iv  )*n ) = work( j+(iv)*n ) - &
     &                                ddot( j-ki-2, t( ki+2, j ), 1, &
     &                                      work( ki+2+(iv)*n ), 1 )
!
                     work( j  +(iv+1)*n ) = work( j+(iv+1)*n ) - &
     &                                ddot( j-ki-2, t( ki+2, j ), 1, &
     &                                      work( ki+2+(iv+1)*n ), 1 )
!
                     work( j+1+(iv  )*n ) = work( j+1+(iv)*n ) - &
     &                                ddot( j-ki-2, t( ki+2, j+1 ), &
     &                                      1, &
     &                                      work( ki+2+(iv)*n ), 1 )
!
                     work( j+1+(iv+1)*n ) = work( j+1+(iv+1)*n ) - &
     &                                ddot( j-ki-2, t( ki+2, j+1 ), &
     &                                      1, &
     &                                      work( ki+2+(iv+1)*n ), 1 )
!
!                    Solve 2-by-2 complex linear equation
!                    [ (T(j,j)   T(j,j+1)  )**T - (wr-i*wi)*I ]*X = SCALE*B
!                    [ (T(j+1,j) T(j+1,j+1))                  ]
!
                     CALL dlaln2( .true., 2, 2, smin, one, t( j, j ), &
     &                            ldt, one, one, work( j+iv*n ), n, wr, &
     &                            -wi, x, 2, scale, xnorm, ierr )
!
!                    Scale if necessary
!
                     IF( scale.NE.one ) THEN
                        CALL dscal( n-ki+1, scale, work(ki+(iv  )*n), &
     &                              1)
                        CALL dscal( n-ki+1, scale, work(ki+(iv+1)*n), &
     &                              1)
                     END IF
                     work( j  +(iv  )*n ) = x( 1, 1 )
                     work( j  +(iv+1)*n ) = x( 1, 2 )
                     work( j+1+(iv  )*n ) = x( 2, 1 )
                     work( j+1+(iv+1)*n ) = x( 2, 2 )
                     vmax = max( abs( x( 1, 1 ) ), abs( x( 1, 2 ) ), &
     &                           abs( x( 2, 1 ) ), abs( x( 2, 2 ) ), &
     &                           vmax )
                     vcrit = bignum / vmax
!
                  END IF
200            CONTINUE
!
!              Copy the vector x or Q*x to VL and normalize.
!
               IF( .NOT.over ) THEN
!                 ------------------------------
!                 no back-transform: copy x to VL and normalize.
                  CALL dcopy( n-ki+1, work( ki + (iv  )*n ), 1, &
     &                        vl( ki, is   ), 1 )
                  CALL dcopy( n-ki+1, work( ki + (iv+1)*n ), 1, &
     &                        vl( ki, is+1 ), 1 )
!
                  emax = zero
                  DO 220 k = ki, n
                     emax = max( emax, abs( vl( k, is   ) )+ &
     &                                 abs( vl( k, is+1 ) ) )
220               CONTINUE
                  remax = one / emax
                  CALL dscal( n-ki+1, remax, vl( ki, is   ), 1 )
                  CALL dscal( n-ki+1, remax, vl( ki, is+1 ), 1 )
!
                  DO 230 k = 1, ki - 1
                     vl( k, is   ) = zero
                     vl( k, is+1 ) = zero
230               CONTINUE
!
               ELSE IF( nb.EQ.1 ) THEN
!                 ------------------------------
!                 version 1: back-transform each vector with GEMV, Q*x.
                  IF( ki.LT.n-1 ) THEN
                     CALL dgemv( 'N', n, n-ki-1, one, &
     &                           vl( 1, ki+2 ), ldvl, &
     &                           work( ki+2 + (iv)*n ), 1, &
     &                           work( ki   + (iv)*n ), &
     &                           vl( 1, ki ), 1 )
                     CALL dgemv( 'N', n, n-ki-1, one, &
     &                           vl( 1, ki+2 ), ldvl, &
     &                           work( ki+2 + (iv+1)*n ), 1, &
     &                           work( ki+1 + (iv+1)*n ), &
     &                           vl( 1, ki+1 ), 1 )
                  ELSE
                     CALL dscal( n, work(ki+  (iv  )*n), vl(1, ki  ), &
     &                           1)
                     CALL dscal( n, work(ki+1+(iv+1)*n), vl(1, ki+1), &
     &                           1)
                  END IF
!
                  emax = zero
                  DO 240 k = 1, n
                     emax = max( emax, abs( vl( k, ki   ) )+ &
     &                                 abs( vl( k, ki+1 ) ) )
240               CONTINUE
                  remax = one / emax
                  CALL dscal( n, remax, vl( 1, ki   ), 1 )
                  CALL dscal( n, remax, vl( 1, ki+1 ), 1 )
!
               ELSE
!                 ------------------------------
!                 version 2: back-transform block of vectors with GEMM
!                 zero out above vector
!                 could go from KI-NV+1 to KI-1
                  DO k = 1, ki - 1
                     work( k + (iv  )*n ) = zero
                     work( k + (iv+1)*n ) = zero
                  END DO
                  iscomplex( iv   ) =  ip
                  iscomplex( iv+1 ) = -ip
                  iv = iv + 1
!                 back-transform and normalization is done below
               END IF
            END IF

            IF( nb.GT.1 ) THEN
!              --------------------------------------------------------
!              Blocked version of back-transform
!              For complex case, KI2 includes both vectors (KI and KI+1)
               IF( ip.EQ.0 ) THEN
                  ki2 = ki
               ELSE
                  ki2 = ki + 1
               END IF

!              Columns 1:IV of work are valid vectors.
!              When the number of vectors stored reaches NB-1 or NB,
!              or if this was last vector, do the GEMM
               IF( (iv.GE.nb-1) .OR. (ki2.EQ.n) ) THEN
                  CALL dgemm( 'N', 'N', n, iv, n-ki2+iv, one, &
     &                        vl( 1, ki2-iv+1 ), ldvl, &
     &                        work( ki2-iv+1 + (1)*n ), n, &
     &                        zero, &
     &                        work( 1 + (nb+1)*n ), n )
!                 normalize vectors
                  DO k = 1, iv
                     IF( iscomplex(k).EQ.0) THEN
!                       real eigenvector
                        ii = idamax( n, work( 1 + (nb+k)*n ), 1 )
                        remax = one / abs( work( ii + (nb+k)*n ) )
                     ELSE IF( iscomplex(k).EQ.1) THEN
!                       first eigenvector of conjugate pair
                        emax = zero
                        DO ii = 1, n
                           emax = max( emax, &
     &                                 abs( work( ii + (nb+k  )*n ) )+ &
     &                                 abs( work( ii + (nb+k+1)*n ) ) )
                        END DO
                        remax = one / emax
!                    else if ISCOMPLEX(K).EQ.-1
!                       second eigenvector of conjugate pair
!                       reuse same REMAX as previous K
                     END IF
                     CALL dscal( n, remax, work( 1 + (nb+k)*n ), 1 )
                  END DO
                  CALL dlacpy( 'F', n, iv, &
     &                         work( 1 + (nb+1)*n ), n, &
     &                         vl( 1, ki2-iv+1 ), ldvl )
                  iv = 1
               ELSE
                  iv = iv + 1
               END IF
            END IF ! blocked back-transform
!
            is = is + 1
            IF( ip.NE.0 ) &
     &         is = is + 1
260      CONTINUE
      END IF
!
      RETURN
!
!     End of DTREVC3
!

      END
!> \brief \b DSCAL
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSCAL(N,DA,DX,INCX)
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION DA
!       INTEGER INCX,N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION DX(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    DSCAL scales a vector by a constant.
!>    uses unrolled loops for increment equal to 1.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         number of elements in input vector(s)
!> \endverbatim
!>
!> \param[in] DA
!> \verbatim
!>          DA is DOUBLE PRECISION
!>           On entry, DA specifies the scalar alpha.
!> \endverbatim
!>
!> \param[in,out] DX
!> \verbatim
!>          DX is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>         storage spacing between elements of DX
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup scal
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, linpack, 3/11/78.
!>     modified 3/93 to return if incx .le. 0.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================

      SUBROUTINE dscal(N,DA,DX,INCX)
!
!  -- Reference BLAS level1 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION DA
      INTEGER INCX,N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER I,M,MP1,NINCX
!     .. Parameters ..
      DOUBLE PRECISION ONE
      parameter(one=1.0d+0)
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC mod
!     ..
      IF (n.LE.0 .OR. incx.LE.0 .OR. da.EQ.one) RETURN
      IF (incx.EQ.1) THEN
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
         m = mod(n,5)
         IF (m.NE.0) THEN
            DO i = 1,m
               dx(i) = da*dx(i)
            END DO
            IF (n.LT.5) RETURN
         END IF
         mp1 = m + 1
         DO i = mp1,n,5
            dx(i) = da*dx(i)
            dx(i+1) = da*dx(i+1)
            dx(i+2) = da*dx(i+2)
            dx(i+3) = da*dx(i+3)
            dx(i+4) = da*dx(i+4)
         END DO
      ELSE
!
!        code for increment not equal to 1
!
         nincx = n*incx
         DO i = 1,nincx,incx
            dx(i) = da*dx(i)
         END DO
      END IF
      RETURN
!
!     End of DSCAL
!

      END
!> \brief \b DROT
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DROT(N,DX,INCX,DY,INCY,C,S)
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION C,S
!       INTEGER INCX,INCY,N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION DX(*),DY(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    DROT applies a plane rotation.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         number of elements in input vector(s)
!> \endverbatim
!>
!> \param[in,out] DX
!> \verbatim
!>          DX is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>         storage spacing between elements of DX
!> \endverbatim
!>
!> \param[in,out] DY
!> \verbatim
!>          DY is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>         storage spacing between elements of DY
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] S
!> \verbatim
!>          S is DOUBLE PRECISION
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup rot
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, linpack, 3/11/78.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================

      SUBROUTINE drot(N,DX,INCX,DY,INCY,C,S)
!
!  -- Reference BLAS level1 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION C,S
      INTEGER INCX,INCY,N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER I,IX,IY
!     ..
      IF (n.LE.0) RETURN
      IF (incx.EQ.1 .AND. incy.EQ.1) THEN
!
!       code for both increments equal to 1
!
         DO i = 1,n
            dtemp = c*dx(i) + s*dy(i)
            dy(i) = c*dy(i) - s*dx(i)
            dx(i) = dtemp
         END DO
      ELSE
!
!       code for unequal increments or equal increments not equal
!         to 1
!
         ix = 1
         iy = 1
         IF (incx.LT.0) ix = (-n+1)*incx + 1
         IF (incy.LT.0) iy = (-n+1)*incy + 1
         DO i = 1,n
            dtemp = c*dx(ix) + s*dy(iy)
            dy(iy) = c*dy(iy) - s*dx(ix)
            dx(ix) = dtemp
            ix = ix + incx
            iy = iy + incy
         END DO
      END IF
      RETURN
!
!     End of DROT
!

      END
!> \brief \b DORGHR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download DORGHR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorghr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorghr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorghr.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       SUBROUTINE DORGHR( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            IHI, ILO, INFO, LDA, LWORK, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DORGHR generates a real orthogonal matrix Q which is defined as the
!> product of IHI-ILO elementary reflectors of order N, as returned by
!> DGEHRD:
!>
!> Q = H(ilo) H(ilo+1) . . . H(ihi-1).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix Q. N >= 0.
!> \endverbatim
!>
!> \param[in] ILO
!> \verbatim
!>          ILO is INTEGER
!> \endverbatim
!>
!> \param[in] IHI
!> \verbatim
!>          IHI is INTEGER
!>
!>          ILO and IHI must have the same values as in the previous call
!>          of DGEHRD. Q is equal to the unit matrix except in the
!>          submatrix Q(ilo+1:ihi,ilo+1:ihi).
!>          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the vectors which define the elementary reflectors,
!>          as returned by DGEHRD.
!>          On exit, the N-by-N orthogonal matrix Q.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION array, dimension (N-1)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by DGEHRD.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK. LWORK >= IHI-ILO.
!>          For optimum performance LWORK >= (IHI-ILO)*NB, where NB is
!>          the optimal blocksize.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup unghr
!
!  =====================================================================

      SUBROUTINE dorghr( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, &
     &                   INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            IHI, ILO, INFO, LDA, LWORK, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      parameter( zero = 0.0d+0, one = 1.0d+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IINFO, J, LWKOPT, NB, NH
!     ..
!     .. External Subroutines ..
!      EXTERNAL           dorgqr, xerbla
!     ..
!     .. External Functions ..
!      INTEGER            ILAENV
!      EXTERNAL           ilaenv
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          max, min
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      info = 0
      nh = ihi - ilo
      lquery = ( lwork.EQ.-1 )
      IF( n.LT.0 ) THEN
         info = -1
      ELSE IF( ilo.LT.1 .OR. ilo.GT.max( 1, n ) ) THEN
         info = -2
      ELSE IF( ihi.LT.min( ilo, n ) .OR. ihi.GT.n ) THEN
         info = -3
      ELSE IF( lda.LT.max( 1, n ) ) THEN
         info = -5
      ELSE IF( lwork.LT.max( 1, nh ) .AND. .NOT.lquery ) THEN
         info = -8
      END IF
!
      IF( info.EQ.0 ) THEN
         nb = ilaenv( 1, 'DORGQR', ' ', nh, nh, nh, -1 )
         lwkopt = max( 1, nh )*nb
         work( 1 ) = lwkopt
      END IF
!
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DORGHR', -info )
         RETURN
      ELSE IF( lquery ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( n.EQ.0 ) THEN
         work( 1 ) = 1
         RETURN
      END IF
!
!     Shift the vectors which define the elementary reflectors one
!     column to the right, and set the first ilo and the last n-ihi
!     rows and columns to those of the unit matrix
!
      DO 40 j = ihi, ilo + 1, -1
         DO 10 i = 1, j - 1
            a( i, j ) = zero
10       CONTINUE
         DO 20 i = j + 1, ihi
            a( i, j ) = a( i, j-1 )
20       CONTINUE
         DO 30 i = ihi + 1, n
            a( i, j ) = zero
30       CONTINUE
40    CONTINUE
      DO 60 j = 1, ilo
         DO 50 i = 1, n
            a( i, j ) = zero
50       CONTINUE
         a( j, j ) = one
60    CONTINUE
      DO 80 j = ihi + 1, n
         DO 70 i = 1, n
            a( i, j ) = zero
70       CONTINUE
         a( j, j ) = one
80    CONTINUE
!
      IF( nh.GT.0 ) THEN
!
!        Generate Q(ilo+1:ihi,ilo+1:ihi)
!
         CALL dorgqr( nh, nh, nh, a( ilo+1, ilo+1 ), lda, tau( ilo ), &
     &                work, lwork, iinfo )
      END IF
      work( 1 ) = lwkopt
      RETURN
!
!     End of DORGHR
!

      END
!> \brief \b DORGQR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download DORGQR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorgqr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorgqr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorgqr.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       SUBROUTINE DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, K, LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DORGQR generates an M-by-N real matrix Q with orthonormal columns,
!> which is defined as the first N columns of a product of K elementary
!> reflectors of order M
!>
!>       Q  =  H(1) H(2) . . . H(k)
!>
!> as returned by DGEQRF.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix Q. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix Q. M >= N >= 0.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of elementary reflectors whose product defines the
!>          matrix Q. N >= K >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the i-th column must contain the vector which
!>          defines the elementary reflector H(i), for i = 1,2,...,k, as
!>          returned by DGEQRF in the first k columns of its array
!>          argument A.
!>          On exit, the M-by-N matrix Q.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The first dimension of the array A. LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by DGEQRF.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK. LWORK >= max(1,N).
!>          For optimum performance LWORK >= N*NB, where NB is the
!>          optimal blocksize.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument has an illegal value
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup ungqr
!
!  =====================================================================

      SUBROUTINE dorgqr( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LWORK, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      parameter( zero = 0.0d+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, J, KI, KK, L, LDWORK, &
     &                   LWKOPT, NB, NBMIN, NX
!     ..
!     .. External Subroutines ..
!      EXTERNAL           dlarfb, dlarft, dorg2r, xerbla
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          max, min
!     ..
!     .. External Functions ..
!      INTEGER            ILAENV
!      EXTERNAL           ilaenv
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      info = 0
      nb = ilaenv( 1, 'DORGQR', ' ', m, n, k, -1 )
      lwkopt = max( 1, n )*nb
      work( 1 ) = lwkopt
      lquery = ( lwork.EQ.-1 )
      IF( m.LT.0 ) THEN
         info = -1
      ELSE IF( n.LT.0 .OR. n.GT.m ) THEN
         info = -2
      ELSE IF( k.LT.0 .OR. k.GT.n ) THEN
         info = -3
      ELSE IF( lda.LT.max( 1, m ) ) THEN
         info = -5
      ELSE IF( lwork.LT.max( 1, n ) .AND. .NOT.lquery ) THEN
         info = -8
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DORGQR', -info )
         RETURN
      ELSE IF( lquery ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( n.LE.0 ) THEN
         work( 1 ) = 1
         RETURN
      END IF
!
      nbmin = 2
      nx = 0
      iws = n
      IF( nb.GT.1 .AND. nb.LT.k ) THEN
!
!        Determine when to cross over from blocked to unblocked code.
!
         nx = max( 0, ilaenv( 3, 'DORGQR', ' ', m, n, k, -1 ) )
         IF( nx.LT.k ) THEN
!
!           Determine if workspace is large enough for blocked code.
!
            ldwork = n
            iws = ldwork*nb
            IF( lwork.LT.iws ) THEN
!
!              Not enough workspace to use optimal NB:  reduce NB and
!              determine the minimum value of NB.
!
               nb = lwork / ldwork
               nbmin = max( 2, ilaenv( 2, 'DORGQR', ' ', m, n, k, &
     &                      -1 ) )
            END IF
         END IF
      END IF
!
      IF( nb.GE.nbmin .AND. nb.LT.k .AND. nx.LT.k ) THEN
!
!        Use blocked code after the last block.
!        The first kk columns are handled by the block method.
!
         ki = ( ( k-nx-1 ) / nb )*nb
         kk = min( k, ki+nb )
!
!        Set A(1:kk,kk+1:n) to zero.
!
         DO 20 j = kk + 1, n
            DO 10 i = 1, kk
               a( i, j ) = zero
10          CONTINUE
20       CONTINUE
      ELSE
         kk = 0
      END IF
!
!     Use unblocked code for the last or only block.
!
      IF( kk.LT.n ) &
     &   CALL dorg2r( m-kk, n-kk, k-kk, a( kk+1, kk+1 ), lda, &
     &                tau( kk+1 ), work, iinfo )
!
      IF( kk.GT.0 ) THEN
!
!        Use blocked code
!
         DO 50 i = ki + 1, 1, -nb
            ib = min( nb, k-i+1 )
            IF( i+ib.LE.n ) THEN
!
!              Form the triangular factor of the block reflector
!              H = H(i) H(i+1) . . . H(i+ib-1)
!
               CALL dlarft( 'Forward', 'Columnwise', m-i+1, ib, &
     &                      a( i, i ), lda, tau( i ), work, ldwork )
!
!              Apply H to A(i:m,i+ib:n) from the left
!
               CALL dlarfb( 'Left', 'No transpose', 'Forward', &
     &                      'Columnwise', m-i+1, n-i-ib+1, ib, &
     &                      a( i, i ), lda, work, ldwork, a( i, i+ib ), &
     &                      lda, work( ib+1 ), ldwork )
            END IF
!
!           Apply H to rows i:m of current block
!
            CALL dorg2r( m-i+1, ib, ib, a( i, i ), lda, tau( i ), &
     &                   work, &
     &                   iinfo )
!
!           Set rows 1:i-1 of current block to zero
!
            DO 40 j = i, i + ib - 1
               DO 30 l = 1, i - 1
                  a( l, j ) = zero
30             CONTINUE
40          CONTINUE
50       CONTINUE
      END IF
!
      work( 1 ) = iws
      RETURN
!
!     End of DORGQR
!

      END
!> \brief \b DORGR2 generates all or part of the orthogonal matrix Q from an RQ factorization determined by sgerqf (unblocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download DORGR2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorgr2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorgr2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorgr2.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       SUBROUTINE DORGR2( M, N, K, A, LDA, TAU, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, K, LDA, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DORGR2 generates an m by n real matrix Q with orthonormal rows,
!> which is defined as the last m rows of a product of k elementary
!> reflectors of order n
!>
!>       Q  =  H(1) H(2) . . . H(k)
!>
!> as returned by DGERQF.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix Q. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix Q. N >= M.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of elementary reflectors whose product defines the
!>          matrix Q. M >= K >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the (m-k+i)-th row must contain the vector which
!>          defines the elementary reflector H(i), for i = 1,2,...,k, as
!>          returned by DGERQF in the last k rows of its array argument
!>          A.
!>          On exit, the m by n matrix Q.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The first dimension of the array A. LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by DGERQF.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (M)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument has an illegal value
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup ungr2
!
!  =====================================================================

      SUBROUTINE dorgr2( M, N, K, A, LDA, TAU, WORK, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      parameter( one = 1.0d+0, zero = 0.0d+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, II, J, L
!     ..
!     .. External Subroutines ..
!      EXTERNAL           dlarf, dscal, xerbla
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          max
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      info = 0
      IF( m.LT.0 ) THEN
         info = -1
      ELSE IF( n.LT.m ) THEN
         info = -2
      ELSE IF( k.LT.0 .OR. k.GT.m ) THEN
         info = -3
      ELSE IF( lda.LT.max( 1, m ) ) THEN
         info = -5
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DORGR2', -info )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( m.LE.0 ) &
     &   RETURN
!
      IF( k.LT.m ) THEN
!
!        Initialise rows 1:m-k to rows of the unit matrix
!
         DO 20 j = 1, n
            DO 10 l = 1, m - k
               a( l, j ) = zero
10          CONTINUE
            IF( j.GT.n-m .AND. j.LE.n-k ) &
     &         a( m-n+j, j ) = one
20       CONTINUE
      END IF
!
      DO 40 i = 1, k
         ii = m - k + i
!
!        Apply H(i) to A(1:m-k+i,1:n-k+i) from the right
!
         !A( II, N-M+II ) = ONE
         CALL dlarf1l( 'Right', ii-1, n-m+ii, a( ii, 1 ), lda, &
     &               tau( i ), &
     &               a, lda, work )
         CALL dscal( n-m+ii-1, -tau( i ), a( ii, 1 ), lda )
         a( ii, n-m+ii ) = one - tau( i )
!
!        Set A(m-k+i,n-k+i+1:n) to zero
!
         DO 30 l = n - m + ii + 1, n
            a( ii, l ) = zero
30       CONTINUE
40    CONTINUE
      RETURN
!
!     End of DORGR2
!

      END
!> \brief \b DLASET initializes the off-diagonal elements and the diagonal elements of a matrix to given values.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download DLASET + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaset.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaset.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaset.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            LDA, M, N
!       DOUBLE PRECISION   ALPHA, BETA
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLASET initializes an m-by-n matrix A to BETA on the diagonal and
!> ALPHA on the offdiagonals.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies the part of the matrix A to be set.
!>          = 'U':      Upper triangular part is set; the strictly lower
!>                      triangular part of A is not changed.
!>          = 'L':      Lower triangular part is set; the strictly upper
!>                      triangular part of A is not changed.
!>          Otherwise:  All of the matrix A is set.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is DOUBLE PRECISION
!>          The constant to which the offdiagonal elements are to be set.
!> \endverbatim
!>
!> \param[in] BETA
!> \verbatim
!>          BETA is DOUBLE PRECISION
!>          The constant to which the diagonal elements are to be set.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On exit, the leading m-by-n submatrix of A is set as follows:
!>
!>          if UPLO = 'U', A(i,j) = ALPHA, 1<=i<=j-1, 1<=j<=n,
!>          if UPLO = 'L', A(i,j) = ALPHA, j+1<=i<=m, 1<=j<=n,
!>          otherwise,     A(i,j) = ALPHA, 1<=i<=m, 1<=j<=n, i.ne.j,
!>
!>          and, for all UPLO, A(i,i) = BETA, 1<=i<=min(m,n).
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup laset
!
!  =====================================================================

      SUBROUTINE dlaset( UPLO, M, N, ALPHA, BETA, A, LDA )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, M, N
      DOUBLE PRECISION   ALPHA, BETA
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
!     ..
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, J
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      EXTERNAL           lsame
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          min
!     ..
!     .. Executable Statements ..
!
      IF( lsame( uplo, 'U' ) ) THEN
!
!        Set the strictly upper triangular or trapezoidal part of the
!        array to ALPHA.
!
         DO 20 j = 2, n
            DO 10 i = 1, min( j-1, m )
               a( i, j ) = alpha
10          CONTINUE
20       CONTINUE
!
      ELSE IF( lsame( uplo, 'L' ) ) THEN
!
!        Set the strictly lower triangular or trapezoidal part of the
!        array to ALPHA.
!
         DO 40 j = 1, min( m, n )
            DO 30 i = j + 1, m
               a( i, j ) = alpha
30          CONTINUE
40       CONTINUE
!
      ELSE
!
!        Set the leading m-by-n submatrix to ALPHA.
!
         DO 60 j = 1, n
            DO 50 i = 1, m
               a( i, j ) = alpha
50          CONTINUE
60       CONTINUE
      END IF
!
!     Set the first min(M,N) diagonal elements to BETA.
!
      DO 70 i = 1, min( m, n )
         a( i, i ) = beta
70    CONTINUE
!
      RETURN
!
!     End of DLASET
!

      END
!> \brief \b DLASCL multiplies a general rectangular matrix by a real scalar defined as cto/cfrom.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download DLASCL + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlascl.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlascl.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlascl.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          TYPE
!       INTEGER            INFO, KL, KU, LDA, M, N
!       DOUBLE PRECISION   CFROM, CTO
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLASCL multiplies the M by N real matrix A by the real scalar
!> CTO/CFROM.  This is done without over/underflow as long as the final
!> result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that
!> A may be full, upper triangular, lower triangular, upper Hessenberg,
!> or banded.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TYPE
!> \verbatim
!>          TYPE is CHARACTER*1
!>          TYPE indices the storage type of the input matrix.
!>          = 'G':  A is a full matrix.
!>          = 'L':  A is a lower triangular matrix.
!>          = 'U':  A is an upper triangular matrix.
!>          = 'H':  A is an upper Hessenberg matrix.
!>          = 'B':  A is a symmetric band matrix with lower bandwidth KL
!>                  and upper bandwidth KU and with the only the lower
!>                  half stored.
!>          = 'Q':  A is a symmetric band matrix with lower bandwidth KL
!>                  and upper bandwidth KU and with the only the upper
!>                  half stored.
!>          = 'Z':  A is a band matrix with lower bandwidth KL and upper
!>                  bandwidth KU. See DGBTRF for storage details.
!> \endverbatim
!>
!> \param[in] KL
!> \verbatim
!>          KL is INTEGER
!>          The lower bandwidth of A.  Referenced only if TYPE = 'B',
!>          'Q' or 'Z'.
!> \endverbatim
!>
!> \param[in] KU
!> \verbatim
!>          KU is INTEGER
!>          The upper bandwidth of A.  Referenced only if TYPE = 'B',
!>          'Q' or 'Z'.
!> \endverbatim
!>
!> \param[in] CFROM
!> \verbatim
!>          CFROM is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] CTO
!> \verbatim
!>          CTO is DOUBLE PRECISION
!>
!>          The matrix A is multiplied by CTO/CFROM. A(I,J) is computed
!>          without over/underflow if the final result CTO*A(I,J)/CFROM
!>          can be represented without over/underflow.  CFROM must be
!>          nonzero.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          The matrix to be multiplied by CTO/CFROM.  See TYPE for the
!>          storage type.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.
!>          If TYPE = 'G', 'L', 'U', 'H', LDA >= max(1,M);
!>             TYPE = 'B', LDA >= KL+1;
!>             TYPE = 'Q', LDA >= KU+1;
!>             TYPE = 'Z', LDA >= 2*KL+KU+1.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          0  - successful exit
!>          <0 - if INFO = -i, the i-th argument had an illegal value.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup lascl
!
!  =====================================================================

      SUBROUTINE dlascl( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, &
     &                   INFO )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          TYPE
      INTEGER            INFO, KL, KU, LDA, M, N
      DOUBLE PRECISION   CFROM, CTO
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      parameter( zero = 0.0d0, one = 1.0d0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            DONE
      INTEGER            I, ITYPE, J, K1, K2, K3, K4
      DOUBLE PRECISION   BIGNUM, CFROM1, CFROMC, CTO1, CTOC, MUL, SMLNUM
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME, DISNAN
!      DOUBLE PRECISION   DLAMCH
!      EXTERNAL           lsame, dlamch, disnan
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          abs, max, min
!     ..
!     .. External Subroutines ..
!      EXTERNAL           xerbla
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      info = 0
!
      IF( lsame( TYPE, 'G' ) ) then
         itype = 0
      ELSE IF( lsame( TYPE, 'L' ) ) then
         itype = 1
      ELSE IF( lsame( TYPE, 'U' ) ) then
         itype = 2
      ELSE IF( lsame( TYPE, 'H' ) ) then
         itype = 3
      ELSE IF( lsame( TYPE, 'B' ) ) then
         itype = 4
      ELSE IF( lsame( TYPE, 'Q' ) ) then
         itype = 5
      ELSE IF( lsame( TYPE, 'Z' ) ) then
         itype = 6
      ELSE
         itype = -1
      END IF
!
      IF( itype.EQ.-1 ) THEN
         info = -1
      ELSE IF( cfrom.EQ.zero .OR. disnan(cfrom) ) THEN
         info = -4
      ELSE IF( disnan(cto) ) THEN
         info = -5
      ELSE IF( m.LT.0 ) THEN
         info = -6
      ELSE IF( n.LT.0 .OR. ( itype.EQ.4 .AND. n.NE.m ) .OR. &
     &         ( itype.EQ.5 .AND. n.NE.m ) ) THEN
         info = -7
      ELSE IF( itype.LE.3 .AND. lda.LT.max( 1, m ) ) THEN
         info = -9
      ELSE IF( itype.GE.4 ) THEN
         IF( kl.LT.0 .OR. kl.GT.max( m-1, 0 ) ) THEN
            info = -2
         ELSE IF( ku.LT.0 .OR. ku.GT.max( n-1, 0 ) .OR. &
     &            ( ( itype.EQ.4 .OR. itype.EQ.5 ) .AND. kl.NE.ku ) ) &
     &             THEN
            info = -3
         ELSE IF( ( itype.EQ.4 .AND. lda.LT.kl+1 ) .OR. &
     &            ( itype.EQ.5 .AND. lda.LT.ku+1 ) .OR. &
     &            ( itype.EQ.6 .AND. lda.LT.2*kl+ku+1 ) ) THEN
            info = -9
         END IF
      END IF
!
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DLASCL', -info )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( n.EQ.0 .OR. m.EQ.0 ) &
     &   RETURN
!
!     Get machine parameters
!
      smlnum = dlamch( 'S' )
      bignum = one / smlnum
!
      cfromc = cfrom
      ctoc = cto
!
10    CONTINUE
      cfrom1 = cfromc*smlnum
      IF( cfrom1.EQ.cfromc ) THEN
!        CFROMC is an inf.  Multiply by a correctly signed zero for
!        finite CTOC, or a NaN if CTOC is infinite.
         mul = ctoc / cfromc
         done = .true.
         cto1 = ctoc
      ELSE
         cto1 = ctoc / bignum
         IF( cto1.EQ.ctoc ) THEN
!           CTOC is either 0 or an inf.  In both cases, CTOC itself
!           serves as the correct multiplication factor.
            mul = ctoc
            done = .true.
            cfromc = one
         ELSE IF( abs( cfrom1 ).GT.abs( ctoc ) .AND. ctoc.NE.zero ) THEN
            mul = smlnum
            done = .false.
            cfromc = cfrom1
         ELSE IF( abs( cto1 ).GT.abs( cfromc ) ) THEN
            mul = bignum
            done = .false.
            ctoc = cto1
         ELSE
            mul = ctoc / cfromc
            done = .true.
            IF (mul .EQ. one) &
     &         RETURN
         END IF
      END IF
!
      IF( itype.EQ.0 ) THEN
!
!        Full matrix
!
         DO 30 j = 1, n
            DO 20 i = 1, m
               a( i, j ) = a( i, j )*mul
20          CONTINUE
30       CONTINUE
!
      ELSE IF( itype.EQ.1 ) THEN
!
!        Lower triangular matrix
!
         DO 50 j = 1, n
            DO 40 i = j, m
               a( i, j ) = a( i, j )*mul
40          CONTINUE
50       CONTINUE
!
      ELSE IF( itype.EQ.2 ) THEN
!
!        Upper triangular matrix
!
         DO 70 j = 1, n
            DO 60 i = 1, min( j, m )
               a( i, j ) = a( i, j )*mul
60          CONTINUE
70       CONTINUE
!
      ELSE IF( itype.EQ.3 ) THEN
!
!        Upper Hessenberg matrix
!
         DO 90 j = 1, n
            DO 80 i = 1, min( j+1, m )
               a( i, j ) = a( i, j )*mul
80          CONTINUE
90       CONTINUE
!
      ELSE IF( itype.EQ.4 ) THEN
!
!        Lower half of a symmetric band matrix
!
         k3 = kl + 1
         k4 = n + 1
         DO 110 j = 1, n
            DO 100 i = 1, min( k3, k4-j )
               a( i, j ) = a( i, j )*mul
100         CONTINUE
110      CONTINUE
!
      ELSE IF( itype.EQ.5 ) THEN
!
!        Upper half of a symmetric band matrix
!
         k1 = ku + 2
         k3 = ku + 1
         DO 130 j = 1, n
            DO 120 i = max( k1-j, 1 ), k3
               a( i, j ) = a( i, j )*mul
120         CONTINUE
130      CONTINUE
!
      ELSE IF( itype.EQ.6 ) THEN
!
!        Band matrix
!
         k1 = kl + ku + 2
         k2 = kl + 1
         k3 = 2*kl + ku + 1
         k4 = kl + ku + 1 + m
         DO 150 j = 1, n
            DO 140 i = max( k1-j, k2 ), min( k3, k4-j )
               a( i, j ) = a( i, j )*mul
140         CONTINUE
150      CONTINUE
!
      END IF
!
      IF( .NOT.done ) &
     &   GO TO 10
!
      RETURN
!
!     End of DLASCL
!

      END
!> \brief \b DORG2R generates all or part of the orthogonal matrix Q from a QR factorization determined by sgeqrf (unblocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download DORG2R + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorg2r.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorg2r.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorg2r.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       SUBROUTINE DORG2R( M, N, K, A, LDA, TAU, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, K, LDA, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DORG2R generates an m by n real matrix Q with orthonormal columns,
!> which is defined as the first n columns of a product of k elementary
!> reflectors of order m
!>
!>       Q  =  H(1) H(2) . . . H(k)
!>
!> as returned by DGEQRF.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix Q. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix Q. M >= N >= 0.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of elementary reflectors whose product defines the
!>          matrix Q. N >= K >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the i-th column must contain the vector which
!>          defines the elementary reflector H(i), for i = 1,2,...,k, as
!>          returned by DGEQRF in the first k columns of its array
!>          argument A.
!>          On exit, the m-by-n matrix Q.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The first dimension of the array A. LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by DGEQRF.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument has an illegal value
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup ung2r
!
!  =====================================================================

      SUBROUTINE dorg2r( M, N, K, A, LDA, TAU, WORK, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      parameter( one = 1.0d+0, zero = 0.0d+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J, L
!     ..
!     .. External Subroutines ..
!      EXTERNAL           dlarf1f, dscal, xerbla
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          max
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      info = 0
      IF( m.LT.0 ) THEN
         info = -1
      ELSE IF( n.LT.0 .OR. n.GT.m ) THEN
         info = -2
      ELSE IF( k.LT.0 .OR. k.GT.n ) THEN
         info = -3
      ELSE IF( lda.LT.max( 1, m ) ) THEN
         info = -5
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DORG2R', -info )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( n.LE.0 ) &
     &   RETURN
!
!     Initialise columns k+1:n to columns of the unit matrix
!
      DO 20 j = k + 1, n
         DO 10 l = 1, m
            a( l, j ) = zero
10       CONTINUE
         a( j, j ) = one
20    CONTINUE
!
      DO 40 i = k, 1, -1
!
!        Apply H(i) to A(i:m,i:n) from the left
!
         IF( i.LT.n ) THEN
            CALL dlarf1f( 'Left', m-i+1, n-i, a( i, i ), 1, tau( i ), &
     &                  a( i, i+1 ), lda, work )
         END IF
         IF( i.LT.m ) &
     &      CALL dscal( m-i, -tau( i ), a( i+1, i ), 1 )
         a( i, i ) = one - tau( i )
!
!        Set A(1:i-1,i) to zero
!
         DO 30 l = 1, i - 1
            a( l, i ) = zero
30       CONTINUE
40    CONTINUE
      RETURN
!
!     End of DORG2R
!

      END
!> \brief \b DHSEQR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download DHSEQR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dhseqr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dhseqr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dhseqr.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       SUBROUTINE DHSEQR( JOB, COMPZ, N, ILO, IHI, H, LDH, WR, WI, Z,
!                          LDZ, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            IHI, ILO, INFO, LDH, LDZ, LWORK, N
!       CHARACTER          COMPZ, JOB
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   H( LDH, * ), WI( * ), WORK( * ), WR( * ),
!      $                   Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    DHSEQR computes the eigenvalues of a Hessenberg matrix H
!>    and, optionally, the matrices T and Z from the Schur decomposition
!>    H = Z T Z**T, where T is an upper quasi-triangular matrix (the
!>    Schur form), and Z is the orthogonal matrix of Schur vectors.
!>
!>    Optionally Z may be postmultiplied into an input orthogonal
!>    matrix Q so that this routine can give the Schur factorization
!>    of a matrix A which has been reduced to the Hessenberg form H
!>    by the orthogonal matrix Q:  A = Q*H*Q**T = (QZ)*T*(QZ)**T.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOB
!> \verbatim
!>          JOB is CHARACTER*1
!>           = 'E':  compute eigenvalues only;
!>           = 'S':  compute eigenvalues and the Schur form T.
!> \endverbatim
!>
!> \param[in] COMPZ
!> \verbatim
!>          COMPZ is CHARACTER*1
!>           = 'N':  no Schur vectors are computed;
!>           = 'I':  Z is initialized to the unit matrix and the matrix Z
!>                   of Schur vectors of H is returned;
!>           = 'V':  Z must contain an orthogonal matrix Q on entry, and
!>                   the product Q*Z is returned.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           The order of the matrix H.  N >= 0.
!> \endverbatim
!>
!> \param[in] ILO
!> \verbatim
!>          ILO is INTEGER
!> \endverbatim
!>
!> \param[in] IHI
!> \verbatim
!>          IHI is INTEGER
!>
!>           It is assumed that H is already upper triangular in rows
!>           and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally
!>           set by a previous call to DGEBAL, and then passed to ZGEHRD
!>           when the matrix output by DGEBAL is reduced to Hessenberg
!>           form. Otherwise ILO and IHI should be set to 1 and N
!>           respectively.  If N > 0, then 1 <= ILO <= IHI <= N.
!>           If N = 0, then ILO = 1 and IHI = 0.
!> \endverbatim
!>
!> \param[in,out] H
!> \verbatim
!>          H is DOUBLE PRECISION array, dimension (LDH,N)
!>           On entry, the upper Hessenberg matrix H.
!>           On exit, if INFO = 0 and JOB = 'S', then H contains the
!>           upper quasi-triangular matrix T from the Schur decomposition
!>           (the Schur form); 2-by-2 diagonal blocks (corresponding to
!>           complex conjugate pairs of eigenvalues) are returned in
!>           standard form, with H(i,i) = H(i+1,i+1) and
!>           H(i+1,i)*H(i,i+1) < 0. If INFO = 0 and JOB = 'E', the
!>           contents of H are unspecified on exit.  (The output value of
!>           H when INFO > 0 is given under the description of INFO
!>           below.)
!>
!>           Unlike earlier versions of DHSEQR, this subroutine may
!>           explicitly H(i,j) = 0 for i > j and j = 1, 2, ... ILO-1
!>           or j = IHI+1, IHI+2, ... N.
!> \endverbatim
!>
!> \param[in] LDH
!> \verbatim
!>          LDH is INTEGER
!>           The leading dimension of the array H. LDH >= max(1,N).
!> \endverbatim
!>
!> \param[out] WR
!> \verbatim
!>          WR is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] WI
!> \verbatim
!>          WI is DOUBLE PRECISION array, dimension (N)
!>
!>           The real and imaginary parts, respectively, of the computed
!>           eigenvalues. If two eigenvalues are computed as a complex
!>           conjugate pair, they are stored in consecutive elements of
!>           WR and WI, say the i-th and (i+1)th, with WI(i) > 0 and
!>           WI(i+1) < 0. If JOB = 'S', the eigenvalues are stored in
!>           the same order as on the diagonal of the Schur form returned
!>           in H, with WR(i) = H(i,i) and, if H(i:i+1,i:i+1) is a 2-by-2
!>           diagonal block, WI(i) = sqrt(-H(i+1,i)*H(i,i+1)) and
!>           WI(i+1) = -WI(i).
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is DOUBLE PRECISION array, dimension (LDZ,N)
!>           If COMPZ = 'N', Z is not referenced.
!>           If COMPZ = 'I', on entry Z need not be set and on exit,
!>           if INFO = 0, Z contains the orthogonal matrix Z of the Schur
!>           vectors of H.  If COMPZ = 'V', on entry Z must contain an
!>           N-by-N matrix Q, which is assumed to be equal to the unit
!>           matrix except for the submatrix Z(ILO:IHI,ILO:IHI). On exit,
!>           if INFO = 0, Z contains Q*Z.
!>           Normally Q is the orthogonal matrix generated by DORGHR
!>           after the call to DGEHRD which formed the Hessenberg matrix
!>           H. (The output value of Z when INFO > 0 is given under
!>           the description of INFO below.)
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>           The leading dimension of the array Z.  if COMPZ = 'I' or
!>           COMPZ = 'V', then LDZ >= MAX(1,N).  Otherwise, LDZ >= 1.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (LWORK)
!>           On exit, if INFO = 0, WORK(1) returns an estimate of
!>           the optimal value for LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>           The dimension of the array WORK.  LWORK >= max(1,N)
!>           is sufficient and delivers very good and sometimes
!>           optimal performance.  However, LWORK as large as 11*N
!>           may be required for optimal performance.  A workspace
!>           query is recommended to determine the optimal workspace
!>           size.
!>
!>           If LWORK = -1, then DHSEQR does a workspace query.
!>           In this case, DHSEQR checks the input parameters and
!>           estimates the optimal workspace size for the given
!>           values of N, ILO and IHI.  The estimate is returned
!>           in WORK(1).  No error message related to LWORK is
!>           issued by XERBLA.  Neither H nor Z are accessed.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>             = 0:  successful exit
!>             < 0:  if INFO = -i, the i-th argument had an illegal
!>                    value
!>             > 0:  if INFO = i, DHSEQR failed to compute all of
!>                the eigenvalues.  Elements 1:ilo-1 and i+1:n of WR
!>                and WI contain those eigenvalues which have been
!>                successfully computed.  (Failures are rare.)
!>
!>                If INFO > 0 and JOB = 'E', then on exit, the
!>                remaining unconverged eigenvalues are the eigen-
!>                values of the upper Hessenberg matrix rows and
!>                columns ILO through INFO of the final, output
!>                value of H.
!>
!>                If INFO > 0 and JOB   = 'S', then on exit
!>
!>           (*)  (initial value of H)*U  = U*(final value of H)
!>
!>                where U is an orthogonal matrix.  The final
!>                value of H is upper Hessenberg and quasi-triangular
!>                in rows and columns INFO+1 through IHI.
!>
!>                If INFO > 0 and COMPZ = 'V', then on exit
!>
!>                  (final value of Z)  =  (initial value of Z)*U
!>
!>                where U is the orthogonal matrix in (*) (regard-
!>                less of the value of JOB.)
!>
!>                If INFO > 0 and COMPZ = 'I', then on exit
!>                      (final value of Z)  = U
!>                where U is the orthogonal matrix in (*) (regard-
!>                less of the value of JOB.)
!>
!>                If INFO > 0 and COMPZ = 'N', then Z is not
!>                accessed.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup hseqr
!
!> \par Contributors:
!  ==================
!>
!>       Karen Braman and Ralph Byers, Department of Mathematics,
!>       University of Kansas, USA
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>             Default values supplied by
!>             ILAENV(ISPEC,'DHSEQR',JOB(:1)//COMPZ(:1),N,ILO,IHI,LWORK).
!>             It is suggested that these defaults be adjusted in order
!>             to attain best performance in each particular
!>             computational environment.
!>
!>            ISPEC=12: The DLAHQR vs DLAQR0 crossover point.
!>                      Default: 75. (Must be at least 11.)
!>
!>            ISPEC=13: Recommended deflation window size.
!>                      This depends on ILO, IHI and NS.  NS is the
!>                      number of simultaneous shifts returned
!>                      by ILAENV(ISPEC=15).  (See ISPEC=15 below.)
!>                      The default for (IHI-ILO+1) <= 500 is NS.
!>                      The default for (IHI-ILO+1) >  500 is 3*NS/2.
!>
!>            ISPEC=14: Nibble crossover point. (See IPARMQ for
!>                      details.)  Default: 14% of deflation window
!>                      size.
!>
!>            ISPEC=15: Number of simultaneous shifts in a multishift
!>                      QR iteration.
!>
!>                      If IHI-ILO+1 is ...
!>
!>                      greater than      ...but less    ... the
!>                      or equal to ...      than        default is
!>
!>                           1               30          NS =   2(+)
!>                          30               60          NS =   4(+)
!>                          60              150          NS =  10(+)
!>                         150              590          NS =  **
!>                         590             3000          NS =  64
!>                        3000             6000          NS = 128
!>                        6000             infinity      NS = 256
!>
!>                  (+)  By default some or all matrices of this order
!>                       are passed to the implicit double shift routine
!>                       DLAHQR and this parameter is ignored.  See
!>                       ISPEC=12 above and comments in IPARMQ for
!>                       details.
!>
!>                 (**)  The asterisks (**) indicate an ad-hoc
!>                       function of N increasing from 10 to 64.
!>
!>            ISPEC=16: Select structured matrix multiply.
!>                      If the number of simultaneous shifts (specified
!>                      by ISPEC=15) is less than 14, then the default
!>                      for ISPEC=16 is 0.  Otherwise the default for
!>                      ISPEC=16 is 2.
!> \endverbatim
!
!> \par References:
!  ================
!>
!>       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
!>       Algorithm Part I: Maintaining Well Focused Shifts, and Level 3
!>       Performance, SIAM Journal of Matrix Analysis, volume 23, pages
!>       929--947, 2002.
!> \n
!>       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
!>       Algorithm Part II: Aggressive Early Deflation, SIAM Journal
!>       of Matrix Analysis, volume 23, pages 948--973, 2002.
!
!  =====================================================================

      SUBROUTINE dhseqr( JOB, COMPZ, N, ILO, IHI, H, LDH, WR, WI, Z, &
     &                   LDZ, WORK, LWORK, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            IHI, ILO, INFO, LDH, LDZ, LWORK, N
      CHARACTER          COMPZ, JOB
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   H( LDH, * ), WI( * ), WORK( * ), WR( * ), &
     &                   z( ldz, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
!
!     ==== Matrices of order NTINY or smaller must be processed by
!     .    DLAHQR because of insufficient subdiagonal scratch space.
!     .    (This is a hard limit.) ====
      INTEGER            NTINY
      parameter( ntiny = 15 )
!
!     ==== NL allocates some local workspace to help small matrices
!     .    through a rare DLAHQR failure.  NL > NTINY = 15 is
!     .    required and NL <= NMIN = ILAENV(ISPEC=12,...) is recom-
!     .    mended.  (The default value of NMIN is 75.)  Using NL = 49
!     .    allows up to six simultaneous shifts and a 16-by-16
!     .    deflation window.  ====
      INTEGER            NL
      parameter( nl = 49 )
      DOUBLE PRECISION   ZERO, ONE
      parameter( zero = 0.0d0, one = 1.0d0 )
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION   HL( NL, NL ), WORKL( NL )
!     ..
!     .. Local Scalars ..
      INTEGER            I, KBOT, NMIN
      LOGICAL            INITZ, LQUERY, WANTT, WANTZ
!     ..
!     .. External Functions ..
!      INTEGER            ILAENV
!      LOGICAL            LSAME
!      EXTERNAL           ilaenv, lsame
!     ..
!     .. External Subroutines ..
!      EXTERNAL           dlacpy, dlahqr, dlaqr0, dlaset, &
!     &                   xerbla
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          dble, max, min
!     ..
!     .. Executable Statements ..
!
!     ==== Decode and check the input parameters. ====
!
      wantt = lsame( job, 'S' )
      initz = lsame( compz, 'I' )
      wantz = initz .OR. lsame( compz, 'V' )
      work( 1 ) = dble( max( 1, n ) )
      lquery = lwork.EQ.-1
!
      info = 0
      IF( .NOT.lsame( job, 'E' ) .AND. .NOT.wantt ) THEN
         info = -1
      ELSE IF( .NOT.lsame( compz, 'N' ) .AND. .NOT.wantz ) THEN
         info = -2
      ELSE IF( n.LT.0 ) THEN
         info = -3
      ELSE IF( ilo.LT.1 .OR. ilo.GT.max( 1, n ) ) THEN
         info = -4
      ELSE IF( ihi.LT.min( ilo, n ) .OR. ihi.GT.n ) THEN
         info = -5
      ELSE IF( ldh.LT.max( 1, n ) ) THEN
         info = -7
      ELSE IF( ldz.LT.1 .OR. ( wantz .AND. ldz.LT.max( 1, n ) ) ) THEN
         info = -11
      ELSE IF( lwork.LT.max( 1, n ) .AND. .NOT.lquery ) THEN
         info = -13
      END IF
!
      IF( info.NE.0 ) THEN
!
!        ==== Quick return in case of invalid argument. ====
!
         CALL xerbla( 'DHSEQR', -info )
         RETURN
!
      ELSE IF( n.EQ.0 ) THEN
!
!        ==== Quick return in case N = 0; nothing to do. ====
!
         RETURN
!
      ELSE IF( lquery ) THEN
!
!        ==== Quick return in case of a workspace query ====
!
         CALL dlaqr0( wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, ilo, &
     &                ihi, z, ldz, work, lwork, info )
!        ==== Ensure reported workspace size is backward-compatible with
!        .    previous LAPACK versions. ====
         work( 1 ) = max( dble( max( 1, n ) ), work( 1 ) )
         RETURN
!
      ELSE
!
!        ==== copy eigenvalues isolated by DGEBAL ====
!
         DO 10 i = 1, ilo - 1
            wr( i ) = h( i, i )
            wi( i ) = zero
10       CONTINUE
         DO 20 i = ihi + 1, n
            wr( i ) = h( i, i )
            wi( i ) = zero
20       CONTINUE
!
!        ==== Initialize Z, if requested ====
!
         IF( initz ) &
     &      CALL dlaset( 'A', n, n, zero, one, z, ldz )
!
!        ==== Quick return if possible ====
!
         IF( ilo.EQ.ihi ) THEN
            wr( ilo ) = h( ilo, ilo )
            wi( ilo ) = zero
            RETURN
         END IF
!
!        ==== DLAHQR/DLAQR0 crossover point ====
!
         nmin = ilaenv( 12, 'DHSEQR', job( : 1 ) // compz( : 1 ), n, &
     &          ilo, ihi, lwork )
         nmin = max( ntiny, nmin )
!
!        ==== DLAQR0 for big matrices; DLAHQR for small ones ====
!
         IF( n.GT.nmin ) THEN
            CALL dlaqr0( wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, &
     &                   ilo, &
     &                   ihi, z, ldz, work, lwork, info )
         ELSE
!
!           ==== Small matrix ====
!
            CALL dlahqr( wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, &
     &                   ilo, &
     &                   ihi, z, ldz, info )
!
            IF( info.GT.0 ) THEN
!
!              ==== A rare DLAHQR failure!  DLAQR0 sometimes succeeds
!              .    when DLAHQR fails. ====
!
               kbot = info
!
               IF( n.GE.nl ) THEN
!
!                 ==== Larger matrices have enough subdiagonal scratch
!                 .    space to call DLAQR0 directly. ====
!
                  CALL dlaqr0( wantt, wantz, n, ilo, kbot, h, ldh, &
     &                         wr, &
     &                         wi, ilo, ihi, z, ldz, work, lwork, info )
!
               ELSE
!
!                 ==== Tiny matrices don't have enough subdiagonal
!                 .    scratch space to benefit from DLAQR0.  Hence,
!                 .    tiny matrices must be copied into a larger
!                 .    array before calling DLAQR0. ====
!
                  CALL dlacpy( 'A', n, n, h, ldh, hl, nl )
                  hl( n+1, n ) = zero
                  CALL dlaset( 'A', nl, nl-n, zero, zero, hl( 1, &
     &                         n+1 ), &
     &                         nl )
                  CALL dlaqr0( wantt, wantz, nl, ilo, kbot, hl, nl, &
     &                         wr, &
     &                         wi, ilo, ihi, z, ldz, workl, nl, info )
                  IF( wantt .OR. info.NE.0 ) &
     &               CALL dlacpy( 'A', n, n, hl, nl, h, ldh )
               END IF
            END IF
         END IF
!
!        ==== Clear out the trash, if necessary. ====
!
         IF( ( wantt .OR. info.NE.0 ) .AND. n.GT.2 ) &
     &      CALL dlaset( 'L', n-2, n-2, zero, zero, h( 3, 1 ), ldh )
!
!        ==== Ensure reported workspace size is backward-compatible with
!        .    previous LAPACK versions. ====
!
         work( 1 ) = max( dble( max( 1, n ) ), work( 1 ) )
      END IF
!
!     ==== End of DHSEQR ====
!

      END
!> \brief \b DGEBAK
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download DGEBAK + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgebak.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgebak.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgebak.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGEBAK( JOB, SIDE, N, ILO, IHI, SCALE, M, V, LDV,
!                          INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOB, SIDE
!       INTEGER            IHI, ILO, INFO, LDV, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   SCALE( * ), V( LDV, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGEBAK forms the right or left eigenvectors of a real general matrix
!> by backward transformation on the computed eigenvectors of the
!> balanced matrix output by DGEBAL.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOB
!> \verbatim
!>          JOB is CHARACTER*1
!>          Specifies the type of backward transformation required:
!>          = 'N': do nothing, return immediately;
!>          = 'P': do backward transformation for permutation only;
!>          = 'S': do backward transformation for scaling only;
!>          = 'B': do backward transformations for both permutation and
!>                 scaling.
!>          JOB must be the same as the argument JOB supplied to DGEBAL.
!> \endverbatim
!>
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'R':  V contains right eigenvectors;
!>          = 'L':  V contains left eigenvectors.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of rows of the matrix V.  N >= 0.
!> \endverbatim
!>
!> \param[in] ILO
!> \verbatim
!>          ILO is INTEGER
!> \endverbatim
!>
!> \param[in] IHI
!> \verbatim
!>          IHI is INTEGER
!>          The integers ILO and IHI determined by DGEBAL.
!>          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
!> \endverbatim
!>
!> \param[in] SCALE
!> \verbatim
!>          SCALE is DOUBLE PRECISION array, dimension (N)
!>          Details of the permutation and scaling factors, as returned
!>          by DGEBAL.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of columns of the matrix V.  M >= 0.
!> \endverbatim
!>
!> \param[in,out] V
!> \verbatim
!>          V is DOUBLE PRECISION array, dimension (LDV,M)
!>          On entry, the matrix of right or left eigenvectors to be
!>          transformed, as returned by DHSEIN or DTREVC.
!>          On exit, V is overwritten by the transformed eigenvectors.
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is INTEGER
!>          The leading dimension of the array V. LDV >= max(1,N).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup gebak
!
!  =====================================================================

      SUBROUTINE dgebak( JOB, SIDE, N, ILO, IHI, SCALE, M, V, LDV, &
     &                   INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          JOB, SIDE
      INTEGER            IHI, ILO, INFO, LDV, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   SCALE( * ), V( LDV, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      parameter( one = 1.0d+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LEFTV, RIGHTV
      INTEGER            I, II, K
      DOUBLE PRECISION   S
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      EXTERNAL           lsame
!     ..
!     .. External Subroutines ..
!      EXTERNAL           dscal, dswap, xerbla
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          max, min
!     ..
!     .. Executable Statements ..
!
!     Decode and Test the input parameters
!
      rightv = lsame( side, 'R' )
      leftv = lsame( side, 'L' )
!
      info = 0
      IF( .NOT.lsame( job, 'N' ) .AND. &
     &    .NOT.lsame( job, 'P' ) .AND. &
     &    .NOT.lsame( job, 'S' ) .AND. &
     &                .NOT.lsame( job, 'B' ) ) THEN
         info = -1
      ELSE IF( .NOT.rightv .AND. .NOT.leftv ) THEN
         info = -2
      ELSE IF( n.LT.0 ) THEN
         info = -3
      ELSE IF( ilo.LT.1 .OR. ilo.GT.max( 1, n ) ) THEN
         info = -4
      ELSE IF( ihi.LT.min( ilo, n ) .OR. ihi.GT.n ) THEN
         info = -5
      ELSE IF( m.LT.0 ) THEN
         info = -7
      ELSE IF( ldv.LT.max( 1, n ) ) THEN
         info = -9
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DGEBAK', -info )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( n.EQ.0 ) &
     &   RETURN
      IF( m.EQ.0 ) &
     &   RETURN
      IF( lsame( job, 'N' ) ) &
     &   RETURN
!
      IF( ilo.EQ.ihi ) &
     &   GO TO 30
!
!     Backward balance
!
      IF( lsame( job, 'S' ) .OR. lsame( job, 'B' ) ) THEN
!
         IF( rightv ) THEN
            DO 10 i = ilo, ihi
               s = scale( i )
               CALL dscal( m, s, v( i, 1 ), ldv )
10          CONTINUE
         END IF
!
         IF( leftv ) THEN
            DO 20 i = ilo, ihi
               s = one / scale( i )
               CALL dscal( m, s, v( i, 1 ), ldv )
20          CONTINUE
         END IF
!
      END IF
!
!     Backward permutation
!
!     For  I = ILO-1 step -1 until 1,
!              IHI+1 step 1 until N do --
!
30    CONTINUE
      IF( lsame( job, 'P' ) .OR. lsame( job, 'B' ) ) THEN
         IF( rightv ) THEN
            DO 40 ii = 1, n
               i = ii
               IF( i.GE.ilo .AND. i.LE.ihi ) &
     &            GO TO 40
               IF( i.LT.ilo ) &
     &            i = ilo - ii
               k = int( scale( i ) )
               IF( k.EQ.i ) &
     &            GO TO 40
               CALL dswap( m, v( i, 1 ), ldv, v( k, 1 ), ldv )
40          CONTINUE
         END IF
!
         IF( leftv ) THEN
            DO 50 ii = 1, n
               i = ii
               IF( i.GE.ilo .AND. i.LE.ihi ) &
     &            GO TO 50
               IF( i.LT.ilo ) &
     &            i = ilo - ii
               k = int( scale( i ) )
               IF( k.EQ.i ) &
     &            GO TO 50
               CALL dswap( m, v( i, 1 ), ldv, v( k, 1 ), ldv )
50          CONTINUE
         END IF
      END IF
!
      RETURN
!
!     End of DGEBAK
!

      END
!> \brief \b DGEHRD
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download DGEHRD + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgehrd.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgehrd.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgehrd.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGEHRD( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            IHI, ILO, INFO, LDA, LWORK, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION  A( LDA, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGEHRD reduces a real general matrix A to upper Hessenberg form H by
!> an orthogonal similarity transformation:  Q**T * A * Q = H .
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] ILO
!> \verbatim
!>          ILO is INTEGER
!> \endverbatim
!>
!> \param[in] IHI
!> \verbatim
!>          IHI is INTEGER
!>
!>          It is assumed that A is already upper triangular in rows
!>          and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally
!>          set by a previous call to DGEBAL; otherwise they should be
!>          set to 1 and N respectively. See Further Details.
!>          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the N-by-N general matrix to be reduced.
!>          On exit, the upper triangle and the first subdiagonal of A
!>          are overwritten with the upper Hessenberg matrix H, and the
!>          elements below the first subdiagonal, with the array TAU,
!>          represent the orthogonal matrix Q as a product of elementary
!>          reflectors. See Further Details.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION array, dimension (N-1)
!>          The scalar factors of the elementary reflectors (see Further
!>          Details). Elements 1:ILO-1 and IHI:N-1 of TAU are set to
!>          zero.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The length of the array WORK.  LWORK >= max(1,N).
!>          For good performance, LWORK should generally be larger.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup gehrd
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The matrix Q is represented as a product of (ihi-ilo) elementary
!>  reflectors
!>
!>     Q = H(ilo) H(ilo+1) . . . H(ihi-1).
!>
!>  Each H(i) has the form
!>
!>     H(i) = I - tau * v * v**T
!>
!>  where tau is a real scalar, and v is a real vector with
!>  v(1:i) = 0, v(i+1) = 1 and v(ihi+1:n) = 0; v(i+2:ihi) is stored on
!>  exit in A(i+2:ihi,i), and tau in TAU(i).
!>
!>  The contents of A are illustrated by the following example, with
!>  n = 7, ilo = 2 and ihi = 6:
!>
!>  on entry,                        on exit,
!>
!>  ( a   a   a   a   a   a   a )    (  a   a   h   h   h   h   a )
!>  (     a   a   a   a   a   a )    (      a   h   h   h   h   a )
!>  (     a   a   a   a   a   a )    (      h   h   h   h   h   h )
!>  (     a   a   a   a   a   a )    (      v2  h   h   h   h   h )
!>  (     a   a   a   a   a   a )    (      v2  v3  h   h   h   h )
!>  (     a   a   a   a   a   a )    (      v2  v3  v4  h   h   h )
!>  (                         a )    (                          a )
!>
!>  where a denotes an element of the original matrix A, h denotes a
!>  modified element of the upper Hessenberg matrix H, and vi denotes an
!>  element of the vector defining H(i).
!>
!>  This file is a slight modification of LAPACK-3.0's DGEHRD
!>  subroutine incorporating improvements proposed by Quintana-Orti and
!>  Van de Geijn (2006). (See DLAHR2.)
!> \endverbatim
!>
!  =====================================================================

      SUBROUTINE dgehrd( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, &
     &                   INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            IHI, ILO, INFO, LDA, LWORK, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            NBMAX, LDT, TSIZE
      parameter( nbmax = 64, ldt = nbmax+1, &
     &                     tsize = ldt*nbmax )
      DOUBLE PRECISION   ZERO, ONE
      parameter( zero = 0.0d+0, &
     &                     one = 1.0d+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWT, J, LDWORK, LWKOPT, NB, &
     &                   nbmin, nh, nx
      DOUBLE PRECISION   EI
!     ..
!     .. External Subroutines ..
!      EXTERNAL           daxpy, dgehd2, dgemm, dlahr2, dlarfb, &
!     &                   dtrmm, &
!     &                   xerbla
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          max, min
!     ..
!     .. External Functions ..
!      INTEGER            ILAENV
!      EXTERNAL           ilaenv
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      info = 0
      lquery = ( lwork.EQ.-1 )
      IF( n.LT.0 ) THEN
         info = -1
      ELSE IF( ilo.LT.1 .OR. ilo.GT.max( 1, n ) ) THEN
         info = -2
      ELSE IF( ihi.LT.min( ilo, n ) .OR. ihi.GT.n ) THEN
         info = -3
      ELSE IF( lda.LT.max( 1, n ) ) THEN
         info = -5
      ELSE IF( lwork.LT.max( 1, n ) .AND. .NOT.lquery ) THEN
         info = -8
      END IF
!
      nh = ihi - ilo + 1
      IF( info.EQ.0 ) THEN
!
!        Compute the workspace requirements
!
         IF( nh.LE.1 ) THEN
            lwkopt = 1
         ELSE
            nb = min( nbmax, ilaenv( 1, 'DGEHRD', ' ', n, ilo, ihi, &
     &                              -1 ) )
            lwkopt = n*nb + tsize
         ENDIF
         work( 1 ) = lwkopt
      END IF
!
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DGEHRD', -info )
         RETURN
      ELSE IF( lquery ) THEN
         RETURN
      END IF
!
!     Set elements 1:ILO-1 and IHI:N-1 of TAU to zero
!
      DO 10 i = 1, ilo - 1
         tau( i ) = zero
10    CONTINUE
      DO 20 i = max( 1, ihi ), n - 1
         tau( i ) = zero
20    CONTINUE
!
!     Quick return if possible
!
      IF( nh.LE.1 ) THEN
         work( 1 ) = 1
         RETURN
      END IF
!
!     Determine the block size
!
      nb = min( nbmax, ilaenv( 1, 'DGEHRD', ' ', n, ilo, ihi, -1 ) )
      nbmin = 2
      IF( nb.GT.1 .AND. nb.LT.nh ) THEN
!
!        Determine when to cross over from blocked to unblocked code
!        (last block is always handled by unblocked code)
!
         nx = max( nb, ilaenv( 3, 'DGEHRD', ' ', n, ilo, ihi, -1 ) )
         IF( nx.LT.nh ) THEN
!
!           Determine if workspace is large enough for blocked code
!
            IF( lwork.LT.lwkopt ) THEN
!
!              Not enough workspace to use optimal NB:  determine the
!              minimum value of NB, and reduce NB or force use of
!              unblocked code
!
               nbmin = max( 2, ilaenv( 2, 'DGEHRD', ' ', n, ilo, ihi, &
     &                 -1 ) )
               IF( lwork.GE.(n*nbmin + tsize) ) THEN
                  nb = (lwork-tsize) / n
               ELSE
                  nb = 1
               END IF
            END IF
         END IF
      END IF
      ldwork = n
!
      IF( nb.LT.nbmin .OR. nb.GE.nh ) THEN
!
!        Use unblocked code below
!
         i = ilo
!
      ELSE
!
!        Use blocked code
!
         iwt = 1 + n*nb
         DO 40 i = ilo, ihi - 1 - nx, nb
            ib = min( nb, ihi-i )
!
!           Reduce columns i:i+ib-1 to Hessenberg form, returning the
!           matrices V and T of the block reflector H = I - V*T*V**T
!           which performs the reduction, and also the matrix Y = A*V*T
!
            CALL dlahr2( ihi, i, ib, a( 1, i ), lda, tau( i ), &
     &                   work( iwt ), ldt, work, ldwork )
!
!           Apply the block reflector H to A(1:ihi,i+ib:ihi) from the
!           right, computing  A := A - Y * V**T. V(i+ib,ib-1) must be set
!           to 1
!
            ei = a( i+ib, i+ib-1 )
            a( i+ib, i+ib-1 ) = one
            CALL dgemm( 'No transpose', 'Transpose', &
     &                  ihi, ihi-i-ib+1, &
     &                  ib, -one, work, ldwork, a( i+ib, i ), lda, one, &
     &                  a( 1, i+ib ), lda )
            a( i+ib, i+ib-1 ) = ei
!
!           Apply the block reflector H to A(1:i,i+1:i+ib-1) from the
!           right
!
            CALL dtrmm( 'Right', 'Lower', 'Transpose', &
     &                  'Unit', i, ib-1, &
     &                  one, a( i+1, i ), lda, work, ldwork )
            DO 30 j = 0, ib-2
               CALL daxpy( i, -one, work( ldwork*j+1 ), 1, &
     &                     a( 1, i+j+1 ), 1 )
30          CONTINUE
!
!           Apply the block reflector H to A(i+1:ihi,i+ib:n) from the
!           left
!
            CALL dlarfb( 'Left', 'Transpose', 'Forward', &
     &                   'Columnwise', &
     &                   ihi-i, n-i-ib+1, ib, a( i+1, i ), lda, &
     &                   work( iwt ), ldt, a( i+1, i+ib ), lda, &
     &                   work, ldwork )
40       CONTINUE
      END IF
!
!     Use unblocked code to reduce the rest of the matrix
!
      CALL dgehd2( n, i, ihi, a, lda, tau, work, iinfo )
!
      work( 1 ) = lwkopt
!
      RETURN
!
!     End of DGEHRD
!

      END
!> \brief \b DLANGE returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute value of any element of a general rectangular matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download DLANGE + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlange.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlange.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlange.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION DLANGE( NORM, M, N, A, LDA, WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          NORM
!       INTEGER            LDA, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLANGE  returns the value of the one norm,  or the Frobenius norm, or
!> the  infinity norm,  or the  element of  largest absolute value  of a
!> real matrix A.
!> \endverbatim
!>
!> \return DLANGE
!> \verbatim
!>
!>    DLANGE = ( max(abs(A(i,j))), NORM = 'M' or 'm'
!>             (
!>             ( norm1(A),         NORM = '1', 'O' or 'o'
!>             (
!>             ( normI(A),         NORM = 'I' or 'i'
!>             (
!>             ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
!>
!> where  norm1  denotes the  one norm of a matrix (maximum column sum),
!> normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
!> normF  denotes the  Frobenius norm of a matrix (square root of sum of
!> squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] NORM
!> \verbatim
!>          NORM is CHARACTER*1
!>          Specifies the value to be returned in DLANGE as described
!>          above.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.  When M = 0,
!>          DLANGE is set to zero.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.  When N = 0,
!>          DLANGE is set to zero.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          The m by n matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(M,1).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)),
!>          where LWORK >= M when NORM = 'I'; otherwise, WORK is not
!>          referenced.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup lange
!
!  =====================================================================

      DOUBLE PRECISION FUNCTION dlange( NORM, M, N, A, LDA, WORK )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          norm
      INTEGER            lda, m, n
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   a( lda, * ), work( * )
!     ..
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   one, zero
      parameter( one = 1.0d+0, zero = 0.0d+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            i, j
      DOUBLE PRECISION   scale, sum, VALUE, temp
!     ..
!     .. External Subroutines ..
!      EXTERNAL           dlassq
!     ..
!     .. External Functions ..
!      LOGICAL            lsame, disnan
!      EXTERNAL           lsame, disnan
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          abs, min, sqrt
!     ..
!     .. Executable Statements ..
!
      IF( min( m, n ).EQ.0 ) THEN
         VALUE = zero
      ELSE IF( lsame( norm, 'M' ) ) THEN
!
!        Find max(abs(A(i,j))).
!
         VALUE = zero
         DO 20 j = 1, n
            DO 10 i = 1, m
               temp = abs( a( i, j ) )
               IF( VALUE.LT.temp .OR. disnan( temp ) ) VALUE = temp
10          CONTINUE
20       CONTINUE
      ELSE IF( ( lsame( norm, 'O' ) ) .OR. ( norm.EQ.'1' ) ) THEN
!
!        Find norm1(A).
!
         VALUE = zero
         DO 40 j = 1, n
            sum = zero
            DO 30 i = 1, m
               sum = sum + abs( a( i, j ) )
30          CONTINUE
            IF( VALUE.LT.sum .OR. disnan( sum ) ) VALUE = sum
40       CONTINUE
      ELSE IF( lsame( norm, 'I' ) ) THEN
!
!        Find normI(A).
!
         DO 50 i = 1, m
            work( i ) = zero
50       CONTINUE
         DO 70 j = 1, n
            DO 60 i = 1, m
               work( i ) = work( i ) + abs( a( i, j ) )
60          CONTINUE
70       CONTINUE
         VALUE = zero
         DO 80 i = 1, m
            temp = work( i )
            IF( VALUE.LT.temp .OR. disnan( temp ) ) VALUE = temp
80       CONTINUE
      ELSE IF( ( lsame( norm, 'F' ) ) .OR. &
     &         ( lsame( norm, 'E' ) ) ) THEN
!
!        Find normF(A).
!
         scale = zero
         sum = one
         DO 90 j = 1, n
            CALL dlassq( m, a( 1, j ), 1, scale, sum )
90       CONTINUE
         VALUE = scale*sqrt( sum )
      END IF
!
      dlange = VALUE
      RETURN
!
!     End of DLANGE
!

      END
!> \brief \b DDOT
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,INCY,N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION DX(*),DY(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    DDOT forms the dot product of two vectors.
!>    uses unrolled loops for increments equal to one.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         number of elements in input vector(s)
!> \endverbatim
!>
!> \param[in] DX
!> \verbatim
!>          DX is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>         storage spacing between elements of DX
!> \endverbatim
!>
!> \param[in] DY
!> \verbatim
!>          DY is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>         storage spacing between elements of DY
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup dot
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, linpack, 3/11/78.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================

      DOUBLE PRECISION FUNCTION ddot(N,DX,INCX,DY,INCY)
!
!  -- Reference BLAS level1 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER incx,incy,n
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION dx(*),dy(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      DOUBLE PRECISION dtemp
      INTEGER i,ix,iy,m,mp1
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC mod
!     ..
      ddot = 0.0d0
      dtemp = 0.0d0
      IF (n.LE.0) RETURN
      IF (incx.EQ.1 .AND. incy.EQ.1) THEN
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
         m = mod(n,5)
         IF (m.NE.0) THEN
            DO i = 1,m
               dtemp = dtemp + dx(i)*dy(i)
            END DO
            IF (n.LT.5) THEN
               ddot=dtemp
            RETURN
            END IF
         END IF
         mp1 = m + 1
         DO i = mp1,n,5
          dtemp = dtemp + dx(i)*dy(i) + dx(i+1)*dy(i+1) + &
     &            dx(i+2)*dy(i+2) + dx(i+3)*dy(i+3) + dx(i+4)*dy(i+4)
         END DO
      ELSE
!
!        code for unequal increments or equal increments
!          not equal to 1
!
         ix = 1
         iy = 1
         IF (incx.LT.0) ix = (-n+1)*incx + 1
         IF (incy.LT.0) iy = (-n+1)*incy + 1
         DO i = 1,n
            dtemp = dtemp + dx(ix)*dy(iy)
            ix = ix + incx
            iy = iy + incy
         END DO
      END IF
      ddot = dtemp
      RETURN
!
!     End of DDOT
!

      END
!> \brief \b DLARF1F applies an elementary reflector to a general rectangular
!              matrix assuming v(1) = 1.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download DLARF1F + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarf1f.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarf1f.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarf1f.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLARF1F( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          SIDE
!       INTEGER            INCV, LDC, M, N
!       DOUBLE PRECISION   TAU
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   C( LDC, * ), V( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLARF1F applies a real elementary reflector H to a real m by n matrix
!> C, from either the left or the right. H is represented in the form
!>
!>       H = I - tau * v * v**T
!>
!> where tau is a real scalar and v is a real vector.
!>
!> If tau = 0, then H is taken to be the unit matrix.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'L': form  H * C
!>          = 'R': form  C * H
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix C.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix C.
!> \endverbatim
!>
!> \param[in] V
!> \verbatim
!>          V is DOUBLE PRECISION array, dimension
!>                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
!>                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
!>          The vector v in the representation of H. V is not used if
!>          TAU = 0. V(1) is not referenced or modified.
!> \endverbatim
!>
!> \param[in] INCV
!> \verbatim
!>          INCV is INTEGER
!>          The increment between elements of v. INCV <> 0.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION
!>          The value tau in the representation of H.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension (LDC,N)
!>          On entry, the m by n matrix C.
!>          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
!>          or C * H if SIDE = 'R'.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C. LDC >= max(1,M).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension
!>                         (N) if SIDE = 'L'
!>                      or (M) if SIDE = 'R'
!> \endverbatim
!
!  To take advantage of the fact that v(1) = 1, we do the following
!     v = [ 1 v_2 ]**T
!     If SIDE='L'
!           |-----|
!           | C_1 |
!        C =| C_2 |
!           |-----|
!        C_1\in\mathbb{R}^{1\times n}, C_2\in\mathbb{R}^{m-1\times n}
!        So we compute:
!        C = HC   = (I - \tau vv**T)C
!                 = C - \tau vv**T C
!        w = C**T v  = [ C_1**T C_2**T ] [ 1 v_2 ]**T
!                    = C_1**T + C_2**T v ( DGEMM then DAXPY )
!        C  = C - \tau vv**T C
!           = C - \tau vw**T
!        Giving us   C_1 = C_1 - \tau w**T ( DAXPY )
!                 and
!                    C_2 = C_2 - \tau v_2w**T ( DGER )
!     If SIDE='R'
!
!        C = [ C_1 C_2 ]
!        C_1\in\mathbb{R}^{m\times 1}, C_2\in\mathbb{R}^{m\times n-1}
!        So we compute: 
!        C = CH   = C(I - \tau vv**T)
!                 = C - \tau Cvv**T
!
!        w = Cv   = [ C_1 C_2 ] [ 1 v_2 ]**T
!                 = C_1 + C_2v_2 ( DGEMM then DAXPY )
!        C  = C - \tau Cvv**T
!           = C - \tau wv**T
!        Giving us   C_1 = C_1 - \tau w ( DAXPY )
!                 and
!                    C_2 = C_2 - \tau wv_2**T ( DGER )
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup larf
!
!  =====================================================================

      SUBROUTINE dlarf1f( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            INCV, LDC, M, N
      DOUBLE PRECISION   TAU
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   C( LDC, * ), V( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      parameter( one = 1.0d+0, zero = 0.0d+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            APPLYLEFT
      INTEGER            I, LASTV, LASTC
!     ..
!     .. External Subroutines ..
!      EXTERNAL           dgemv, dger, daxpy, dscal
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      INTEGER            ILADLR, ILADLC
!      EXTERNAL           lsame, iladlr, iladlc
!     ..
!     .. Executable Statements ..
!
      applyleft = lsame( side, 'L' )
      lastv = 1
      lastc = 0
      IF( tau.NE.zero ) THEN
!     Set up variables for scanning V.  LASTV begins pointing to the end
!     of V.
         IF( applyleft ) THEN
            lastv = m
         ELSE
            lastv = n
         END IF
         IF( incv.GT.0 ) THEN
            i = 1 + (lastv-1) * incv
         ELSE
            i = 1
         END IF
!     Look for the last non-zero row in V.
!        Since we are assuming that V(1) = 1, and it is not stored, so we
!        shouldn't access it.
         DO WHILE( lastv.GT.1 .AND. v( i ).EQ.zero )
            lastv = lastv - 1
            i = i - incv
         END DO
         IF( applyleft ) THEN
!     Scan for the last non-zero column in C(1:lastv,:).
            lastc = iladlc(lastv, n, c, ldc)
         ELSE
!     Scan for the last non-zero row in C(:,1:lastv).
            lastc = iladlr(m, lastv, c, ldc)
         END IF
      END IF
      IF( lastc.EQ.0 ) THEN
         RETURN
      END IF
      IF( applyleft ) THEN
!
!        Form  H * C
!
         ! Check if lastv = 1. This means v = 1, So we just need to comp
         ! C := HC = (1-\tau)C.
         IF( lastv.EQ.1 ) THEN
!
!           C(1,1:lastc) := ( 1 - tau ) * C(1,1:lastc)
!
            CALL dscal(lastc, one - tau, c, ldc)
         ELSE
!
!           w(1:lastc,1) := C(1:lastv,1:lastc)**T * v(1:lastv,1)
!
            ! w(1:lastc,1) := C(2:lastv,1:lastc)**T * v(2:lastv,1)
            CALL dgemv( 'Transpose', lastv-1, lastc, one, c(1+1,1), &
     &                  ldc, v(1+incv), incv, zero, work, 1)
            ! w(1:lastc,1) += C(1,1:lastc)**T * v(1,1) = C(1,1:lastc)**T
            CALL daxpy(lastc, one, c, ldc, work, 1)
!
!           C(1:lastv,1:lastc) := C(...) - tau * v(1:lastv,1) * w(1:lastc,1)**T
!
            ! C(1, 1:lastc)   := C(...) - tau * v(1,1) * w(1:lastc,1)**T
            !                  = C(...) - tau * w(1:lastc,1)**T
            CALL daxpy(lastc, -tau, work, 1, c, ldc)
            ! C(2:lastv,1:lastc) := C(...) - tau * v(2:lastv,1)*w(1:last
            CALL dger(lastv-1, lastc, -tau, v(1+incv), incv, work, 1, &
     &                  c(1+1,1), ldc)
         END IF
      ELSE
!
!        Form  C * H
!
         ! Check if n = 1. This means v = 1, so we just need to compute
         ! C := CH = C(1-\tau).
         IF( lastv.EQ.1 ) THEN
!
!           C(1:lastc,1) := ( 1 - tau ) * C(1:lastc,1)
!
            CALL dscal(lastc, one - tau, c, 1)
         ELSE
!
!           w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1)
!
            ! w(1:lastc,1) := C(1:lastc,2:lastv) * v(2:lastv,1)
            CALL dgemv( 'No transpose', lastc, lastv-1, one, &
     &         c(1,1+1), ldc, v(1+incv), incv, zero, work, 1 )
            ! w(1:lastc,1) += C(1:lastc,1) v(1,1) = C(1:lastc,1)
            CALL daxpy(lastc, one, c, 1, work, 1)
!
!           C(1:lastc,1:lastv) := C(...) - tau * w(1:lastc,1) * v(1:lastv,1)**T
!
            ! C(1:lastc,1)     := C(...) - tau * w(1:lastc,1) * v(1,1)**
            !                   = C(...) - tau * w(1:lastc,1)
            CALL daxpy(lastc, -tau, work, 1, c, 1)
            ! C(1:lastc,2:lastv) := C(...) - tau * w(1:lastc,1) * v(2:la
            CALL dger( lastc, lastv-1, -tau, work, 1, v(1+incv), &
     &                  incv, c(1,1+1), ldc )
         END IF
      END IF
      RETURN
!
!     End of DLARF1F
!

      END
!> \brief \b DGER
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION ALPHA
!       INTEGER INCX,INCY,LDA,M,N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION A(LDA,*),X(*),Y(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGER   performs the rank 1 operation
!>
!>    A := alpha*x*y**T + A,
!>
!> where alpha is a scalar, x is an m element vector, y is an n element
!> vector and A is an m by n matrix.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>           On entry, M specifies the number of rows of the matrix A.
!>           M must be at least zero.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry, N specifies the number of columns of the matrix A.
!>           N must be at least zero.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is DOUBLE PRECISION.
!>           On entry, ALPHA specifies the scalar alpha.
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension at least
!>           ( 1 + ( m - 1 )*abs( INCX ) ).
!>           Before entry, the incremented array X must contain the m
!>           element vector x.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>           On entry, INCX specifies the increment for the elements of
!>           X. INCX must not be zero.
!> \endverbatim
!>
!> \param[in] Y
!> \verbatim
!>          Y is DOUBLE PRECISION array, dimension at least
!>           ( 1 + ( n - 1 )*abs( INCY ) ).
!>           Before entry, the incremented array Y must contain the n
!>           element vector y.
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>           On entry, INCY specifies the increment for the elements of
!>           Y. INCY must not be zero.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension ( LDA, N )
!>           Before entry, the leading m by n part of the array A must
!>           contain the matrix of coefficients. On exit, A is
!>           overwritten by the updated matrix.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in the calling (sub) program. LDA must be at least
!>           max( 1, m ).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup ger
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Level 2 Blas routine.
!>
!>  -- Written on 22-October-1986.
!>     Jack Dongarra, Argonne National Lab.
!>     Jeremy Du Croz, Nag Central Office.
!>     Sven Hammarling, Nag Central Office.
!>     Richard Hanson, Sandia National Labs.
!> \endverbatim
!>
!  =====================================================================

      SUBROUTINE dger(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
!
!  -- Reference BLAS level2 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA
      INTEGER INCX,INCY,LDA,M,N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),X(*),Y(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO
      parameter(zero=0.0d+0)
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,J,JY,KX
!     ..
!     .. External Subroutines ..
!      EXTERNAL xerbla
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC max
!     ..
!
!     Test the input parameters.
!
      info = 0
      IF (m.LT.0) THEN
          info = 1
      ELSE IF (n.LT.0) THEN
          info = 2
      ELSE IF (incx.EQ.0) THEN
          info = 5
      ELSE IF (incy.EQ.0) THEN
          info = 7
      ELSE IF (lda.LT.max(1,m)) THEN
          info = 9
      END IF
      IF (info.NE.0) THEN
          CALL xerbla('DGER  ',info)
          RETURN
      END IF
!
!     Quick return if possible.
!
      IF ((m.EQ.0) .OR. (n.EQ.0) .OR. (alpha.EQ.zero)) RETURN
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
      IF (incy.GT.0) THEN
          jy = 1
      ELSE
          jy = 1 - (n-1)*incy
      END IF
      IF (incx.EQ.1) THEN
          DO 20 j = 1,n
              IF (y(jy).NE.zero) THEN
                  temp = alpha*y(jy)
                  DO 10 i = 1,m
                      a(i,j) = a(i,j) + x(i)*temp
10                CONTINUE
              END IF
              jy = jy + incy
20        CONTINUE
      ELSE
          IF (incx.GT.0) THEN
              kx = 1
          ELSE
              kx = 1 - (m-1)*incx
          END IF
          DO 40 j = 1,n
              IF (y(jy).NE.zero) THEN
                  temp = alpha*y(jy)
                  ix = kx
                  DO 30 i = 1,m
                      a(i,j) = a(i,j) + x(ix)*temp
                      ix = ix + incx
30                CONTINUE
              END IF
              jy = jy + incy
40        CONTINUE
      END IF
!
      RETURN
!
!     End of DGER
!

      END
!> \brief \b DAXPY
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION DA
!       INTEGER INCX,INCY,N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION DX(*),DY(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    DAXPY constant times a vector plus a vector.
!>    uses unrolled loops for increments equal to one.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         number of elements in input vector(s)
!> \endverbatim
!>
!> \param[in] DA
!> \verbatim
!>          DA is DOUBLE PRECISION
!>           On entry, DA specifies the scalar alpha.
!> \endverbatim
!>
!> \param[in] DX
!> \verbatim
!>          DX is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>         storage spacing between elements of DX
!> \endverbatim
!>
!> \param[in,out] DY
!> \verbatim
!>          DY is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>         storage spacing between elements of DY
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup axpy
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, linpack, 3/11/78.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================

      SUBROUTINE daxpy(N,DA,DX,INCX,DY,INCY)
!
!  -- Reference BLAS level1 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION DA
      INTEGER INCX,INCY,N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER I,IX,IY,M,MP1
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC mod
!     ..
      IF (n.LE.0) RETURN
      IF (da.EQ.0.0d0) RETURN
      IF (incx.EQ.1 .AND. incy.EQ.1) THEN
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
         m = mod(n,4)
         IF (m.NE.0) THEN
            DO i = 1,m
               dy(i) = dy(i) + da*dx(i)
            END DO
         END IF
         IF (n.LT.4) RETURN
         mp1 = m + 1
         DO i = mp1,n,4
            dy(i) = dy(i) + da*dx(i)
            dy(i+1) = dy(i+1) + da*dx(i+1)
            dy(i+2) = dy(i+2) + da*dx(i+2)
            dy(i+3) = dy(i+3) + da*dx(i+3)
         END DO
      ELSE
!
!        code for unequal increments or equal increments
!          not equal to 1
!
         ix = 1
         iy = 1
         IF (incx.LT.0) ix = (-n+1)*incx + 1
         IF (incy.LT.0) iy = (-n+1)*incy + 1
         DO i = 1,n
          dy(iy) = dy(iy) + da*dx(ix)
          ix = ix + incx
          iy = iy + incy
         END DO
      END IF
      RETURN
!
!     End of DAXPY
!

      END
!> \brief \b DGEBAL
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download DGEBAL + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgebal.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgebal.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgebal.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGEBAL( JOB, N, A, LDA, ILO, IHI, SCALE, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOB
!       INTEGER            IHI, ILO, INFO, LDA, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), SCALE( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGEBAL balances a general real matrix A.  This involves, first,
!> permuting A by a similarity transformation to isolate eigenvalues
!> in the first 1 to ILO-1 and last IHI+1 to N elements on the
!> diagonal; and second, applying a diagonal similarity transformation
!> to rows and columns ILO to IHI to make the rows and columns as
!> close in norm as possible.  Both steps are optional.
!>
!> Balancing may reduce the 1-norm of the matrix, and improve the
!> accuracy of the computed eigenvalues and/or eigenvectors.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOB
!> \verbatim
!>          JOB is CHARACTER*1
!>          Specifies the operations to be performed on A:
!>          = 'N':  none:  simply set ILO = 1, IHI = N, SCALE(I) = 1.0
!>                  for i = 1,...,N;
!>          = 'P':  permute only;
!>          = 'S':  scale only;
!>          = 'B':  both permute and scale.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the input matrix A.
!>          On exit,  A is overwritten by the balanced matrix.
!>          If JOB = 'N', A is not referenced.
!>          See Further Details.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] ILO
!> \verbatim
!>          ILO is INTEGER
!> \endverbatim
!> \param[out] IHI
!> \verbatim
!>          IHI is INTEGER
!>          ILO and IHI are set to integers such that on exit
!>          A(i,j) = 0 if i > j and j = 1,...,ILO-1 or I = IHI+1,...,N.
!>          If JOB = 'N' or 'S', ILO = 1 and IHI = N.
!> \endverbatim
!>
!> \param[out] SCALE
!> \verbatim
!>          SCALE is DOUBLE PRECISION array, dimension (N)
!>          Details of the permutations and scaling factors applied to
!>          A.  If P(j) is the index of the row and column interchanged
!>          with row and column j and D(j) is the scaling factor
!>          applied to row and column j, then
!>          SCALE(j) = P(j)    for j = 1,...,ILO-1
!>                   = D(j)    for j = ILO,...,IHI
!>                   = P(j)    for j = IHI+1,...,N.
!>          The order in which the interchanges are made is N to IHI+1,
!>          then 1 to ILO-1.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup gebal
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The permutations consist of row and column interchanges which put
!>  the matrix in the form
!>
!>             ( T1   X   Y  )
!>     P A P = (  0   B   Z  )
!>             (  0   0   T2 )
!>
!>  where T1 and T2 are upper triangular matrices whose eigenvalues lie
!>  along the diagonal.  The column indices ILO and IHI mark the starting
!>  and ending columns of the submatrix B. Balancing consists of applying
!>  a diagonal similarity transformation inv(D) * B * D to make the
!>  1-norms of each row of B and its corresponding column nearly equal.
!>  The output matrix is
!>
!>     ( T1     X*D          Y    )
!>     (  0  inv(D)*B*D  inv(D)*Z ).
!>     (  0      0           T2   )
!>
!>  Information about the permutations P and the diagonal matrix D is
!>  returned in the vector SCALE.
!>
!>  This subroutine is based on the EISPACK routine BALANC.
!>
!>  Modified by Tzu-Yi Chen, Computer Science Division, University of
!>    California at Berkeley, USA
!>
!>  Refactored by Evert Provoost, Department of Computer Science,
!>    KU Leuven, Belgium
!> \endverbatim
!>
!  =====================================================================

      SUBROUTINE dgebal( JOB, N, A, LDA, ILO, IHI, SCALE, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          JOB
      INTEGER            IHI, ILO, INFO, LDA, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), SCALE( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      parameter( zero = 0.0d+0, one = 1.0d+0 )
      DOUBLE PRECISION   SCLFAC
      parameter( sclfac = 2.0d+0 )
      DOUBLE PRECISION   FACTOR
      parameter( factor = 0.95d+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            NOCONV, CANSWAP
      INTEGER            I, ICA, IRA, J, K, L
      DOUBLE PRECISION   C, CA, F, G, R, RA, S, SFMAX1, SFMAX2, SFMIN1, &
     &                   SFMIN2
!     ..
!     .. External Functions ..
!      LOGICAL            DISNAN, LSAME
!      INTEGER            IDAMAX
!      DOUBLE PRECISION   DLAMCH, DNRM2
!      EXTERNAL           disnan, lsame, idamax, dlamch, &
!     &                   dnrm2
!     ..
!     .. External Subroutines ..
!      EXTERNAL           dscal, dswap, xerbla
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          abs, max, min
!     ..
!     Test the input parameters
!
      info = 0
      IF( .NOT.lsame( job, 'N' ) .AND. &
     &    .NOT.lsame( job, 'P' ) .AND. &
     &    .NOT.lsame( job, 'S' ) .AND. &
     &                .NOT.lsame( job, 'B' ) ) THEN
         info = -1
      ELSE IF( n.LT.0 ) THEN
         info = -2
      ELSE IF( lda.LT.max( 1, n ) ) THEN
         info = -4
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DGEBAL', -info )
         RETURN
      END IF
!
!     Quick returns.
!
      IF( n.EQ.0 ) THEN
         ilo = 1
         ihi = 0
         RETURN
      END IF
!
      IF( lsame( job, 'N' ) ) THEN
         DO i = 1, n
            scale( i ) = one
         END DO
         ilo = 1
         ihi = n
         RETURN
      END IF
!
!     Permutation to isolate eigenvalues if possible.
!
      k = 1
      l = n
!
      IF( .NOT.lsame( job, 'S' ) ) THEN
!
!        Row and column exchange.
!
         noconv = .true.
         DO WHILE( noconv )
!
!           Search for rows isolating an eigenvalue and push them down.
!
            noconv = .false.
            DO i = l, 1, -1
               canswap = .true.
               DO j = 1, l
                  IF( i.NE.j .AND. a( i, j ).NE.zero ) THEN
                     canswap = .false.
                     EXIT
                  END IF
               END DO
!
               IF( canswap ) THEN
                  scale( l ) = i
                  IF( i.NE.l ) THEN
                     CALL dswap( l, a( 1, i ), 1, a( 1, l ), 1 )
                     CALL dswap( n-k+1, a( i, k ), lda, a( l, k ), &
     &                           lda )
                  END IF
                  noconv = .true.
!
                  IF( l.EQ.1 ) THEN
                     ilo = 1
                     ihi = 1
                     RETURN
                  END IF
!
                  l = l - 1
               END IF
            END DO
!
         END DO

         noconv = .true.
         DO WHILE( noconv )
!
!           Search for columns isolating an eigenvalue and push them left.
!
            noconv = .false.
            DO j = k, l
               canswap = .true.
               DO i = k, l
                  IF( i.NE.j .AND. a( i, j ).NE.zero ) THEN
                     canswap = .false.
                     EXIT
                  END IF
               END DO
!
               IF( canswap ) THEN
                  scale( k ) = j
                  IF( j.NE.k ) THEN
                     CALL dswap( l, a( 1, j ), 1, a( 1, k ), 1 )
                     CALL dswap( n-k+1, a( j, k ), lda, a( k, k ), &
     &                           lda )
                  END IF
                  noconv = .true.
!
                  k = k + 1
               END IF
            END DO
!
         END DO
!
      END IF
!
!     Initialize SCALE for non-permuted submatrix.
!
      DO i = k, l
         scale( i ) = one
      END DO
!
!     If we only had to permute, we are done.
!
      IF( lsame( job, 'P' ) ) THEN
         ilo = k
         ihi = l
         RETURN
      END IF
!
!     Balance the submatrix in rows K to L.
!
!     Iterative loop for norm reduction.
!
      sfmin1 = dlamch( 'S' ) / dlamch( 'P' )
      sfmax1 = one / sfmin1
      sfmin2 = sfmin1*sclfac
      sfmax2 = one / sfmin2
!
      noconv = .true.
      DO WHILE( noconv )
         noconv = .false.
!
         DO i = k, l
!
            c = dnrm2( l-k+1, a( k, i ), 1 )
            r = dnrm2( l-k+1, a( i, k ), lda )
            ica = idamax( l, a( 1, i ), 1 )
            ca = abs( a( ica, i ) )
            ira = idamax( n-k+1, a( i, k ), lda )
            ra = abs( a( i, ira+k-1 ) )
!
!           Guard against zero C or R due to underflow.
!
            IF( c.EQ.zero .OR. r.EQ.zero ) cycle
!
!           Exit if NaN to avoid infinite loop
!
            IF( disnan( c+ca+r+ra ) ) THEN
               info = -3
               CALL xerbla( 'DGEBAL', -info )
               RETURN
            END IF
!
            g = r / sclfac
            f = one
            s = c + r
!
            DO WHILE( c.LT.g .AND. max( f, c, ca ).LT.sfmax2 .AND. &
     &                min( r, g, ra ).GT.sfmin2 )
               f = f*sclfac
               c = c*sclfac
               ca = ca*sclfac
               r = r / sclfac
               g = g / sclfac
               ra = ra / sclfac
            END DO
!
            g = c / sclfac
!
            DO WHILE( g.GE.r .AND. max( r, ra ).LT.sfmax2 .AND. &
     &                min( f, c, g, ca ).GT.sfmin2 )
               f = f / sclfac
               c = c / sclfac
               g = g / sclfac
               ca = ca / sclfac
               r = r*sclfac
               ra = ra*sclfac
            END DO
!
!           Now balance.
!
            IF( ( c+r ).GE.factor*s ) cycle
            IF( f.LT.one .AND. scale( i ).LT.one ) THEN
               IF( f*scale( i ).LE.sfmin1 ) cycle
            END IF
            IF( f.GT.one .AND. scale( i ).GT.one ) THEN
               IF( scale( i ).GE.sfmax1 / f ) cycle
            END IF
            g = one / f
            scale( i ) = scale( i )*f
            noconv = .true.
!
            CALL dscal( n-k+1, g, a( i, k ), lda )
            CALL dscal( l, f, a( 1, i ), 1 )
!
         END DO
!
      END DO
!
      ilo = k
      ihi = l
!
      RETURN
!
!     End of DGEBAL
!

      END
!> \brief \b DLAHQR computes the eigenvalues and Schur factorization of an upper Hessenberg matrix, using the double-shift/single-shift QR algorithm.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download DLAHQR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlahqr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlahqr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlahqr.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI,
!                          ILOZ, IHIZ, Z, LDZ, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N
!       LOGICAL            WANTT, WANTZ
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   H( LDH, * ), WI( * ), WR( * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    DLAHQR is an auxiliary routine called by DHSEQR to update the
!>    eigenvalues and Schur decomposition already computed by DHSEQR, by
!>    dealing with the Hessenberg submatrix in rows and columns ILO to
!>    IHI.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] WANTT
!> \verbatim
!>          WANTT is LOGICAL
!>          = .TRUE. : the full Schur form T is required;
!>          = .FALSE.: only eigenvalues are required.
!> \endverbatim
!>
!> \param[in] WANTZ
!> \verbatim
!>          WANTZ is LOGICAL
!>          = .TRUE. : the matrix of Schur vectors Z is required;
!>          = .FALSE.: Schur vectors are not required.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix H.  N >= 0.
!> \endverbatim
!>
!> \param[in] ILO
!> \verbatim
!>          ILO is INTEGER
!> \endverbatim
!>
!> \param[in] IHI
!> \verbatim
!>          IHI is INTEGER
!>          It is assumed that H is already upper quasi-triangular in
!>          rows and columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless
!>          ILO = 1). DLAHQR works primarily with the Hessenberg
!>          submatrix in rows and columns ILO to IHI, but applies
!>          transformations to all of H if WANTT is .TRUE..
!>          1 <= ILO <= max(1,IHI); IHI <= N.
!> \endverbatim
!>
!> \param[in,out] H
!> \verbatim
!>          H is DOUBLE PRECISION array, dimension (LDH,N)
!>          On entry, the upper Hessenberg matrix H.
!>          On exit, if INFO is zero and if WANTT is .TRUE., H is upper
!>          quasi-triangular in rows and columns ILO:IHI, with any
!>          2-by-2 diagonal blocks in standard form. If INFO is zero
!>          and WANTT is .FALSE., the contents of H are unspecified on
!>          exit.  The output state of H if INFO is nonzero is given
!>          below under the description of INFO.
!> \endverbatim
!>
!> \param[in] LDH
!> \verbatim
!>          LDH is INTEGER
!>          The leading dimension of the array H. LDH >= max(1,N).
!> \endverbatim
!>
!> \param[out] WR
!> \verbatim
!>          WR is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] WI
!> \verbatim
!>          WI is DOUBLE PRECISION array, dimension (N)
!>          The real and imaginary parts, respectively, of the computed
!>          eigenvalues ILO to IHI are stored in the corresponding
!>          elements of WR and WI. If two eigenvalues are computed as a
!>          complex conjugate pair, they are stored in consecutive
!>          elements of WR and WI, say the i-th and (i+1)th, with
!>          WI(i) > 0 and WI(i+1) < 0. If WANTT is .TRUE., the
!>          eigenvalues are stored in the same order as on the diagonal
!>          of the Schur form returned in H, with WR(i) = H(i,i), and, if
!>          H(i:i+1,i:i+1) is a 2-by-2 diagonal block,
!>          WI(i) = sqrt(H(i+1,i)*H(i,i+1)) and WI(i+1) = -WI(i).
!> \endverbatim
!>
!> \param[in] ILOZ
!> \verbatim
!>          ILOZ is INTEGER
!> \endverbatim
!>
!> \param[in] IHIZ
!> \verbatim
!>          IHIZ is INTEGER
!>          Specify the rows of Z to which transformations must be
!>          applied if WANTZ is .TRUE..
!>          1 <= ILOZ <= ILO; IHI <= IHIZ <= N.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is DOUBLE PRECISION array, dimension (LDZ,N)
!>          If WANTZ is .TRUE., on entry Z must contain the current
!>          matrix Z of transformations accumulated by DHSEQR, and on
!>          exit Z has been updated; transformations are applied only to
!>          the submatrix Z(ILOZ:IHIZ,ILO:IHI).
!>          If WANTZ is .FALSE., Z is not referenced.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of the array Z. LDZ >= max(1,N).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>           = 0:  successful exit
!>           > 0:  If INFO = i, DLAHQR failed to compute all the
!>                  eigenvalues ILO to IHI in a total of 30 iterations
!>                  per eigenvalue; elements i+1:ihi of WR and WI
!>                  contain those eigenvalues which have been
!>                  successfully computed.
!>
!>                  If INFO > 0 and WANTT is .FALSE., then on exit,
!>                  the remaining unconverged eigenvalues are the
!>                  eigenvalues of the upper Hessenberg matrix rows
!>                  and columns ILO through INFO of the final, output
!>                  value of H.
!>
!>                  If INFO > 0 and WANTT is .TRUE., then on exit
!>          (*)       (initial value of H)*U  = U*(final value of H)
!>                  where U is an orthogonal matrix.    The final
!>                  value of H is upper Hessenberg and triangular in
!>                  rows and columns INFO+1 through IHI.
!>
!>                  If INFO > 0 and WANTZ is .TRUE., then on exit
!>                      (final value of Z)  = (initial value of Z)*U
!>                  where U is the orthogonal matrix in (*)
!>                  (regardless of the value of WANTT.)
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup lahqr
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>     02-96 Based on modifications by
!>     David Day, Sandia National Laboratory, USA
!>
!>     12-04 Further modifications by
!>     Ralph Byers, University of Kansas, USA
!>     This is a modified version of DLAHQR from LAPACK version 3.0.
!>     It is (1) more robust against overflow and underflow and
!>     (2) adopts the more conservative Ahues & Tisseur stopping
!>     criterion (LAWN 122, 1997).
!> \endverbatim
!>
!  =====================================================================

!   crap etienne dlahqr

!> \brief \b ILADLC scans a matrix for its last non-zero column.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download ILADLC + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/iladlc.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/iladlc.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/iladlc.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       INTEGER FUNCTION ILADLC( M, N, A, LDA )
!
!       .. Scalar Arguments ..
!       INTEGER            M, N, LDA
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ILADLC scans A for its last non-zero column.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          The m by n matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >= max(1,M).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup ilalc
!
!  =====================================================================

      INTEGER FUNCTION iladlc( M, N, A, LDA )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            m, n, lda
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   a( lda, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION zero
      parameter( zero = 0.0d+0 )
!     ..
!     .. Local Scalars ..
      INTEGER i
!     ..
!     .. Executable Statements ..
!
!     Quick test for the common case where one corner is non-zero.
      IF( n.EQ.0 ) THEN
         iladlc = n
      ELSE IF( a(1, n).NE.zero .OR. a(m, n).NE.zero ) THEN
         iladlc = n
      ELSE
!     Now scan each column from the end, returning with the first non-zero.
         DO iladlc = n, 1, -1
            DO i = 1, m
               IF( a(i, iladlc).NE.zero ) RETURN
            END DO
         END DO
      END IF
      RETURN

      END
!> \brief \b DTRMM
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION ALPHA
!       INTEGER LDA,LDB,M,N
!       CHARACTER DIAG,SIDE,TRANSA,UPLO
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION A(LDA,*),B(LDB,*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DTRMM  performs one of the matrix-matrix operations
!>
!>    B := alpha*op( A )*B,   or   B := alpha*B*op( A ),
!>
!> where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
!> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!>
!>    op( A ) = A   or   op( A ) = A**T.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>           On entry,  SIDE specifies whether  op( A ) multiplies B from
!>           the left or right as follows:
!>
!>              SIDE = 'L' or 'l'   B := alpha*op( A )*B.
!>
!>              SIDE = 'R' or 'r'   B := alpha*B*op( A ).
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>           On entry, UPLO specifies whether the matrix A is an upper or
!>           lower triangular matrix as follows:
!>
!>              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!>
!>              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!> \endverbatim
!>
!> \param[in] TRANSA
!> \verbatim
!>          TRANSA is CHARACTER*1
!>           On entry, TRANSA specifies the form of op( A ) to be used in
!>           the matrix multiplication as follows:
!>
!>              TRANSA = 'N' or 'n'   op( A ) = A.
!>
!>              TRANSA = 'T' or 't'   op( A ) = A**T.
!>
!>              TRANSA = 'C' or 'c'   op( A ) = A**T.
!> \endverbatim
!>
!> \param[in] DIAG
!> \verbatim
!>          DIAG is CHARACTER*1
!>           On entry, DIAG specifies whether or not A is unit triangular
!>           as follows:
!>
!>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!>
!>              DIAG = 'N' or 'n'   A is not assumed to be unit
!>                                  triangular.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>           On entry, M specifies the number of rows of B. M must be at
!>           least zero.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry, N specifies the number of columns of B.  N must be
!>           at least zero.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is DOUBLE PRECISION.
!>           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!>           zero then  A is not referenced and  B need not be set before
!>           entry.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>           A is DOUBLE PRECISION array, dimension ( LDA, k ), where k is m
!>           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
!>           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
!>           upper triangular part of the array  A must contain the upper
!>           triangular matrix  and the strictly lower triangular part of
!>           A is not referenced.
!>           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
!>           lower triangular part of the array  A must contain the lower
!>           triangular matrix  and the strictly upper triangular part of
!>           A is not referenced.
!>           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
!>           A  are not referenced either,  but are assumed to be  unity.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!>           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
!>           then LDA must be at least max( 1, n ).
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension ( LDB, N )
!>           Before entry,  the leading  m by n part of the array  B must
!>           contain the matrix  B,  and  on exit  is overwritten  by the
!>           transformed matrix.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>           On entry, LDB specifies the first dimension of B as declared
!>           in  the  calling  (sub)  program.   LDB  must  be  at  least
!>           max( 1, m ).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup trmm
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Level 3 Blas routine.
!>
!>  -- Written on 8-February-1989.
!>     Jack Dongarra, Argonne National Laboratory.
!>     Iain Duff, AERE Harwell.
!>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!>     Sven Hammarling, Numerical Algorithms Group Ltd.
!> \endverbatim
!>
!  =====================================================================

      SUBROUTINE dtrmm(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
!
!  -- Reference BLAS level3 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA
      INTEGER LDA,LDB,M,N
      CHARACTER DIAG,SIDE,TRANSA,UPLO
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),B(LDB,*)
!     ..
!
!  =====================================================================
!
!     .. External Functions ..
!      LOGICAL LSAME
!      EXTERNAL lsame
!     ..
!     .. External Subroutines ..
!      EXTERNAL xerbla
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC max
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,J,K,NROWA
      LOGICAL LSIDE,NOUNIT,UPPER
!     ..
!     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      parameter(one=1.0d+0,zero=0.0d+0)
!     ..
!
!     Test the input parameters.
!
      lside = lsame(side,'L')
      IF (lside) THEN
          nrowa = m
      ELSE
          nrowa = n
      END IF
      nounit = lsame(diag,'N')
      upper = lsame(uplo,'U')
!
      info = 0
      IF ((.NOT.lside) .AND. (.NOT.lsame(side,'R'))) THEN
          info = 1
      ELSE IF ((.NOT.upper) .AND. (.NOT.lsame(uplo,'L'))) THEN
          info = 2
      ELSE IF ((.NOT.lsame(transa,'N')) .AND. &
     &         (.NOT.lsame(transa,'T')) .AND. &
     &         (.NOT.lsame(transa,'C'))) THEN
          info = 3
      ELSE IF ((.NOT.lsame(diag,'U')) .AND. &
     &         (.NOT.lsame(diag,'N'))) THEN
          info = 4
      ELSE IF (m.LT.0) THEN
          info = 5
      ELSE IF (n.LT.0) THEN
          info = 6
      ELSE IF (lda.LT.max(1,nrowa)) THEN
          info = 9
      ELSE IF (ldb.LT.max(1,m)) THEN
          info = 11
      END IF
      IF (info.NE.0) THEN
          CALL xerbla('DTRMM ',info)
          RETURN
      END IF
!
!     Quick return if possible.
!
      IF (m.EQ.0 .OR. n.EQ.0) RETURN
!
!     And when  alpha.eq.zero.
!
      IF (alpha.EQ.zero) THEN
          DO 20 j = 1,n
              DO 10 i = 1,m
                  b(i,j) = zero
10            CONTINUE
20        CONTINUE
          RETURN
      END IF
!
!     Start the operations.
!
      IF (lside) THEN
          IF (lsame(transa,'N')) THEN
!
!           Form  B := alpha*A*B.
!
              IF (upper) THEN
                  DO 50 j = 1,n
                      DO 40 k = 1,m
                          IF (b(k,j).NE.zero) THEN
                              temp = alpha*b(k,j)
                              DO 30 i = 1,k - 1
                                  b(i,j) = b(i,j) + temp*a(i,k)
30                            CONTINUE
                              IF (nounit) temp = temp*a(k,k)
                              b(k,j) = temp
                          END IF
40                    CONTINUE
50                CONTINUE
              ELSE
                  DO 80 j = 1,n
                      DO 70 k = m,1,-1
                          IF (b(k,j).NE.zero) THEN
                              temp = alpha*b(k,j)
                              b(k,j) = temp
                              IF (nounit) b(k,j) = b(k,j)*a(k,k)
                              DO 60 i = k + 1,m
                                  b(i,j) = b(i,j) + temp*a(i,k)
60                            CONTINUE
                          END IF
70                    CONTINUE
80                CONTINUE
              END IF
          ELSE
!
!           Form  B := alpha*A**T*B.
!
              IF (upper) THEN
                  DO 110 j = 1,n
                      DO 100 i = m,1,-1
                          temp = b(i,j)
                          IF (nounit) temp = temp*a(i,i)
                          DO 90 k = 1,i - 1
                              temp = temp + a(k,i)*b(k,j)
90                        CONTINUE
                          b(i,j) = alpha*temp
100                   CONTINUE
110               CONTINUE
              ELSE
                  DO 140 j = 1,n
                      DO 130 i = 1,m
                          temp = b(i,j)
                          IF (nounit) temp = temp*a(i,i)
                          DO 120 k = i + 1,m
                              temp = temp + a(k,i)*b(k,j)
120                       CONTINUE
                          b(i,j) = alpha*temp
130                   CONTINUE
140               CONTINUE
              END IF
          END IF
      ELSE
          IF (lsame(transa,'N')) THEN
!
!           Form  B := alpha*B*A.
!
              IF (upper) THEN
                  DO 180 j = n,1,-1
                      temp = alpha
                      IF (nounit) temp = temp*a(j,j)
                      DO 150 i = 1,m
                          b(i,j) = temp*b(i,j)
150                   CONTINUE
                      DO 170 k = 1,j - 1
                          IF (a(k,j).NE.zero) THEN
                              temp = alpha*a(k,j)
                              DO 160 i = 1,m
                                  b(i,j) = b(i,j) + temp*b(i,k)
160                           CONTINUE
                          END IF
170                   CONTINUE
180               CONTINUE
              ELSE
                  DO 220 j = 1,n
                      temp = alpha
                      IF (nounit) temp = temp*a(j,j)
                      DO 190 i = 1,m
                          b(i,j) = temp*b(i,j)
190                   CONTINUE
                      DO 210 k = j + 1,n
                          IF (a(k,j).NE.zero) THEN
                              temp = alpha*a(k,j)
                              DO 200 i = 1,m
                                  b(i,j) = b(i,j) + temp*b(i,k)
200                           CONTINUE
                          END IF
210                   CONTINUE
220               CONTINUE
              END IF
          ELSE
!
!           Form  B := alpha*B*A**T.
!
              IF (upper) THEN
                  DO 260 k = 1,n
                      DO 240 j = 1,k - 1
                          IF (a(j,k).NE.zero) THEN
                              temp = alpha*a(j,k)
                              DO 230 i = 1,m
                                  b(i,j) = b(i,j) + temp*b(i,k)
230                           CONTINUE
                          END IF
240                   CONTINUE
                      temp = alpha
                      IF (nounit) temp = temp*a(k,k)
                      IF (temp.NE.one) THEN
                          DO 250 i = 1,m
                              b(i,k) = temp*b(i,k)
250                       CONTINUE
                      END IF
260               CONTINUE
              ELSE
                  DO 300 k = n,1,-1
                      DO 280 j = k + 1,n
                          IF (a(j,k).NE.zero) THEN
                              temp = alpha*a(j,k)
                              DO 270 i = 1,m
                                  b(i,j) = b(i,j) + temp*b(i,k)
270                           CONTINUE
                          END IF
280                   CONTINUE
                      temp = alpha
                      IF (nounit) temp = temp*a(k,k)
                      IF (temp.NE.one) THEN
                          DO 290 i = 1,m
                              b(i,k) = temp*b(i,k)
290                       CONTINUE
                      END IF
300               CONTINUE
              END IF
          END IF
      END IF
!
      RETURN
!
!     End of DTRMM
!

      END
!> \brief \b DLAHQR computes the eigenvalues and Schur factorization of an upper Hessenberg matrix, using the double-shift/single-shift QR algorithm.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download DLAHQR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlahqr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlahqr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlahqr.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI,
!                          ILOZ, IHIZ, Z, LDZ, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N
!       LOGICAL            WANTT, WANTZ
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   H( LDH, * ), WI( * ), WR( * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    DLAHQR is an auxiliary routine called by DHSEQR to update the
!>    eigenvalues and Schur decomposition already computed by DHSEQR, by
!>    dealing with the Hessenberg submatrix in rows and columns ILO to
!>    IHI.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] WANTT
!> \verbatim
!>          WANTT is LOGICAL
!>          = .TRUE. : the full Schur form T is required;
!>          = .FALSE.: only eigenvalues are required.
!> \endverbatim
!>
!> \param[in] WANTZ
!> \verbatim
!>          WANTZ is LOGICAL
!>          = .TRUE. : the matrix of Schur vectors Z is required;
!>          = .FALSE.: Schur vectors are not required.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix H.  N >= 0.
!> \endverbatim
!>
!> \param[in] ILO
!> \verbatim
!>          ILO is INTEGER
!> \endverbatim
!>
!> \param[in] IHI
!> \verbatim
!>          IHI is INTEGER
!>          It is assumed that H is already upper quasi-triangular in
!>          rows and columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless
!>          ILO = 1). DLAHQR works primarily with the Hessenberg
!>          submatrix in rows and columns ILO to IHI, but applies
!>          transformations to all of H if WANTT is .TRUE..
!>          1 <= ILO <= max(1,IHI); IHI <= N.
!> \endverbatim
!>
!> \param[in,out] H
!> \verbatim
!>          H is DOUBLE PRECISION array, dimension (LDH,N)
!>          On entry, the upper Hessenberg matrix H.
!>          On exit, if INFO is zero and if WANTT is .TRUE., H is upper
!>          quasi-triangular in rows and columns ILO:IHI, with any
!>          2-by-2 diagonal blocks in standard form. If INFO is zero
!>          and WANTT is .FALSE., the contents of H are unspecified on
!>          exit.  The output state of H if INFO is nonzero is given
!>          below under the description of INFO.
!> \endverbatim
!>
!> \param[in] LDH
!> \verbatim
!>          LDH is INTEGER
!>          The leading dimension of the array H. LDH >= max(1,N).
!> \endverbatim
!>
!> \param[out] WR
!> \verbatim
!>          WR is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] WI
!> \verbatim
!>          WI is DOUBLE PRECISION array, dimension (N)
!>          The real and imaginary parts, respectively, of the computed
!>          eigenvalues ILO to IHI are stored in the corresponding
!>          elements of WR and WI. If two eigenvalues are computed as a
!>          complex conjugate pair, they are stored in consecutive
!>          elements of WR and WI, say the i-th and (i+1)th, with
!>          WI(i) > 0 and WI(i+1) < 0. If WANTT is .TRUE., the
!>          eigenvalues are stored in the same order as on the diagonal
!>          of the Schur form returned in H, with WR(i) = H(i,i), and, if
!>          H(i:i+1,i:i+1) is a 2-by-2 diagonal block,
!>          WI(i) = sqrt(H(i+1,i)*H(i,i+1)) and WI(i+1) = -WI(i).
!> \endverbatim
!>
!> \param[in] ILOZ
!> \verbatim
!>          ILOZ is INTEGER
!> \endverbatim
!>
!> \param[in] IHIZ
!> \verbatim
!>          IHIZ is INTEGER
!>          Specify the rows of Z to which transformations must be
!>          applied if WANTZ is .TRUE..
!>          1 <= ILOZ <= ILO; IHI <= IHIZ <= N.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is DOUBLE PRECISION array, dimension (LDZ,N)
!>          If WANTZ is .TRUE., on entry Z must contain the current
!>          matrix Z of transformations accumulated by DHSEQR, and on
!>          exit Z has been updated; transformations are applied only to
!>          the submatrix Z(ILOZ:IHIZ,ILO:IHI).
!>          If WANTZ is .FALSE., Z is not referenced.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of the array Z. LDZ >= max(1,N).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>           = 0:  successful exit
!>           > 0:  If INFO = i, DLAHQR failed to compute all the
!>                  eigenvalues ILO to IHI in a total of 30 iterations
!>                  per eigenvalue; elements i+1:ihi of WR and WI
!>                  contain those eigenvalues which have been
!>                  successfully computed.
!>
!>                  If INFO > 0 and WANTT is .FALSE., then on exit,
!>                  the remaining unconverged eigenvalues are the
!>                  eigenvalues of the upper Hessenberg matrix rows
!>                  and columns ILO through INFO of the final, output
!>                  value of H.
!>
!>                  If INFO > 0 and WANTT is .TRUE., then on exit
!>          (*)       (initial value of H)*U  = U*(final value of H)
!>                  where U is an orthogonal matrix.    The final
!>                  value of H is upper Hessenberg and triangular in
!>                  rows and columns INFO+1 through IHI.
!>
!>                  If INFO > 0 and WANTZ is .TRUE., then on exit
!>                      (final value of Z)  = (initial value of Z)*U
!>                  where U is the orthogonal matrix in (*)
!>                  (regardless of the value of WANTT.)
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup lahqr
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>     02-96 Based on modifications by
!>     David Day, Sandia National Laboratory, USA
!>
!>     12-04 Further modifications by
!>     Ralph Byers, University of Kansas, USA
!>     This is a modified version of DLAHQR from LAPACK version 3.0.
!>     It is (1) more robust against overflow and underflow and
!>     (2) adopts the more conservative Ahues & Tisseur stopping
!>     criterion (LAWN 122, 1997).
!> \endverbatim
!>
!  =====================================================================

      SUBROUTINE dlahqr( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, &
     &                   ILOZ, IHIZ, Z, LDZ, INFO )
      IMPLICIT NONE
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N
      LOGICAL            WANTT, WANTZ
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   H( LDH, * ), WI( * ), WR( * ), Z( LDZ, * )
!     ..
!
!  =========================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO
      parameter( zero = 0.0d0, one = 1.0d0, two = 2.0d0 )
      DOUBLE PRECISION   DAT1, DAT2
      parameter( dat1 = 3.0d0 / 4.0d0, dat2 = -0.4375d0 )
      INTEGER            KEXSH
      parameter( kexsh = 10 )
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION   AA, AB, BA, BB, CS, DET, H11, H12, H21, H21S, &
     &                   h22, rt1i, rt1r, rt2i, rt2r, rtdisc, s, safmax, &
     &                   safmin, smlnum, sn, sum, t1, t2, t3, tr, tst, &
     &                   ulp, v2, v3
      INTEGER            I, I1, I2, ITS, ITMAX, J, K, L, M, NH, NR, NZ, &
     &                   kdefl
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION   V( 3 )
!     ..
!     .. External Functions ..
!      DOUBLE PRECISION   DLAMCH
!      EXTERNAL           dlamch
!     ..
!     .. External Subroutines ..
!      EXTERNAL           dcopy, dlanv2, dlarfg, drot
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          abs, dble, max, min, sqrt
!     ..
!     .. Executable Statements ..
!
      info = 0
!
!     Quick return if possible
!
      IF( n.EQ.0 ) &
     &   RETURN
      IF( ilo.EQ.ihi ) THEN
         wr( ilo ) = h( ilo, ilo )
         wi( ilo ) = zero
         RETURN
      END IF
!
!     ==== clear out the trash ====
      DO 10 j = ilo, ihi - 3
         h( j+2, j ) = zero
         h( j+3, j ) = zero
10    CONTINUE
      IF( ilo.LE.ihi-2 ) &
     &   h( ihi, ihi-2 ) = zero
!
      nh = ihi - ilo + 1
      nz = ihiz - iloz + 1
!
!     Set machine-dependent constants for the stopping criterion.
!
      safmin = dlamch( 'SAFE MINIMUM' )
      safmax = one / safmin
      ulp = dlamch( 'PRECISION' )
      smlnum = safmin*( dble( nh ) / ulp )
!
!     I1 and I2 are the indices of the first row and last column of H
!     to which transformations must be applied. If eigenvalues only are
!     being computed, I1 and I2 are set inside the main loop.
!
      IF( wantt ) THEN
         i1 = 1
         i2 = n
      END IF
!
!     ITMAX is the total number of QR iterations allowed.
!
      itmax = 30 * max( 10, nh )
!
!     KDEFL counts the number of iterations since a deflation
!
      kdefl = 0
!
!     The main loop begins here. I is the loop index and decreases from
!     IHI to ILO in steps of 1 or 2. Each iteration of the loop works
!     with the active submatrix in rows and columns L to I.
!     Eigenvalues I+1 to IHI have already converged. Either L = ILO or
!     H(L,L-1) is negligible so that the matrix splits.
!
      i = ihi
20    CONTINUE
      l = ilo
      IF( i.LT.ilo ) &
     &   GO TO 160
!
!     Perform QR iterations on rows and columns ILO to I until a
!     submatrix of order 1 or 2 splits off at the bottom because a
!     subdiagonal element has become negligible.
!
      DO 140 its = 0, itmax
!
!        Look for a single small subdiagonal element.
!
         DO 30 k = i, l + 1, -1
            IF( abs( h( k, k-1 ) ).LE.smlnum ) &
     &         GO TO 40
            tst = abs( h( k-1, k-1 ) ) + abs( h( k, k ) )
            IF( tst.EQ.zero ) THEN
               IF( k-2.GE.ilo ) &
     &            tst = tst + abs( h( k-1, k-2 ) )
               IF( k+1.LE.ihi ) &
     &            tst = tst + abs( h( k+1, k ) )
            END IF
!           ==== The following is a conservative small subdiagonal
!           .    deflation  criterion due to Ahues & Tisseur (LAWN 122,
!           .    1997). It has better mathematical foundation and
!           .    improves accuracy in some cases.  ====
            IF( abs( h( k, k-1 ) ).LE.ulp*tst ) THEN
               ab = max( abs( h( k, k-1 ) ), abs( h( k-1, k ) ) )
               ba = min( abs( h( k, k-1 ) ), abs( h( k-1, k ) ) )
               aa = max( abs( h( k, k ) ), &
     &              abs( h( k-1, k-1 )-h( k, k ) ) )
               bb = min( abs( h( k, k ) ), &
     &              abs( h( k-1, k-1 )-h( k, k ) ) )
               s = aa + ab
               IF( ba*( ab / s ).LE.max( smlnum, &
     &             ulp*( bb*( aa / s ) ) ) )GO TO 40
            END IF
30       CONTINUE
40       CONTINUE
         l = k
         IF( l.GT.ilo ) THEN
!
!           H(L,L-1) is negligible
!
            h( l, l-1 ) = zero
         END IF
!
!        Exit from loop if a submatrix of order 1 or 2 has split off.
!
         IF( l.GE.i-1 ) &
     &      GO TO 150
         kdefl = kdefl + 1
!
!        Now the active submatrix is in rows and columns L to I. If
!        eigenvalues only are being computed, only the active submatrix
!        need be transformed.
!
         IF( .NOT.wantt ) THEN
            i1 = l
            i2 = i
         END IF
!
         IF( mod(kdefl,2*kexsh).EQ.0 ) THEN
!
!           Exceptional shift.
!
            s = abs( h( i, i-1 ) ) + abs( h( i-1, i-2 ) )
            h11 = dat1*s + h( i, i )
            h12 = dat2*s
            h21 = s
            h22 = h11
         ELSE IF( mod(kdefl,kexsh).EQ.0 ) THEN
!
!           Exceptional shift.
!
            s = abs( h( l+1, l ) ) + abs( h( l+2, l+1 ) )
            h11 = dat1*s + h( l, l )
            h12 = dat2*s
            h21 = s
            h22 = h11
         ELSE
!
!           Prepare to use Francis' double shift
!           (i.e. 2nd degree generalized Rayleigh quotient)
!
            h11 = h( i-1, i-1 )
            h21 = h( i, i-1 )
            h12 = h( i-1, i )
            h22 = h( i, i )
         END IF
         s = abs( h11 ) + abs( h12 ) + abs( h21 ) + abs( h22 )
         IF( s.EQ.zero ) THEN
            rt1r = zero
            rt1i = zero
            rt2r = zero
            rt2i = zero
         ELSE
            h11 = h11 / s
            h21 = h21 / s
            h12 = h12 / s
            h22 = h22 / s
            tr = ( h11+h22 ) / two
            det = ( h11-tr )*( h22-tr ) - h12*h21
            rtdisc = sqrt( abs( det ) )
            IF( det.GE.zero ) THEN
!
!              ==== complex conjugate shifts ====
!
               rt1r = tr*s
               rt2r = rt1r
               rt1i = rtdisc*s
               rt2i = -rt1i
            ELSE
!
!              ==== real shifts (use only one of them)  ====
!
               rt1r = tr + rtdisc
               rt2r = tr - rtdisc
               IF( abs( rt1r-h22 ).LE.abs( rt2r-h22 ) ) THEN
                  rt1r = rt1r*s
                  rt2r = rt1r
               ELSE
                  rt2r = rt2r*s
                  rt1r = rt2r
               END IF
               rt1i = zero
               rt2i = zero
            END IF
         END IF
!
!        Look for two consecutive small subdiagonal elements.
!
         DO 50 m = i - 2, l, -1
!           Determine the effect of starting the double-shift QR
!           iteration at row M, and see if this would make H(M,M-1)
!           negligible.  (The following uses scaling to avoid
!           overflows and most underflows.)
!
            h21s = h( m+1, m )
            s = abs( h( m, m )-rt2r ) + abs( rt2i ) + abs( h21s )
            h21s = h( m+1, m ) / s
            v( 1 ) = h21s*h( m, m+1 ) + ( h( m, m )-rt1r )* &
     &               ( ( h( m, m )-rt2r ) / s ) - rt1i*( rt2i / s )
            v( 2 ) = h21s*( h( m, m )+h( m+1, m+1 )-rt1r-rt2r )
            v( 3 ) = h21s*h( m+2, m+1 )
            s = abs( v( 1 ) ) + abs( v( 2 ) ) + abs( v( 3 ) )
            v( 1 ) = v( 1 ) / s
            v( 2 ) = v( 2 ) / s
            v( 3 ) = v( 3 ) / s
            IF( m.EQ.l ) &
     &         GO TO 60
            IF( abs( h( m, m-1 ) )*( abs( v( 2 ) )+abs( v( 3 ) ) ).LE. &
     &          ulp*abs( v( 1 ) )*( abs( h( m-1, m-1 ) )+abs( h( m, &
     &          m ) )+abs( h( m+1, m+1 ) ) ) )GO TO 60
50       CONTINUE
60       CONTINUE
!
!        Double-shift QR step
!
         DO 130 k = m, i - 1
!
!           The first iteration of this loop determines a reflection G
!           from the vector V and applies it from left and right to H,
!           thus creating a nonzero bulge below the subdiagonal.
!
!           Each subsequent iteration determines a reflection G to
!           restore the Hessenberg form in the (K-1)th column, and thus
!           chases the bulge one step toward the bottom of the active
!           submatrix. NR is the order of G.
!
            nr = min( 3, i-k+1 )
            IF( k.GT.m ) &
     &         CALL dcopy( nr, h( k, k-1 ), 1, v, 1 )
            CALL dlarfg( nr, v( 1 ), v( 2 ), 1, t1 )
            IF( k.GT.m ) THEN
               h( k, k-1 ) = v( 1 )
               h( k+1, k-1 ) = zero
               IF( k.LT.i-1 ) &
     &            h( k+2, k-1 ) = zero
            ELSE IF( m.GT.l ) THEN
!               ==== Use the following instead of
!               .    H( K, K-1 ) = -H( K, K-1 ) to
!               .    avoid a bug when v(2) and v(3)
!               .    underflow. ====
               h( k, k-1 ) = h( k, k-1 )*( one-t1 )
            END IF
            v2 = v( 2 )
            t2 = t1*v2
            IF( nr.EQ.3 ) THEN
               v3 = v( 3 )
               t3 = t1*v3
!
!              Apply G from the left to transform the rows of the matrix
!              in columns K to I2.
!
               DO 70 j = k, i2
                  sum = h( k, j ) + v2*h( k+1, j ) + v3*h( k+2, j )
                  h( k, j ) = h( k, j ) - sum*t1
                  h( k+1, j ) = h( k+1, j ) - sum*t2
                  h( k+2, j ) = h( k+2, j ) - sum*t3
70             CONTINUE
!
!              Apply G from the right to transform the columns of the
!              matrix in rows I1 to min(K+3,I).
!
               DO 80 j = i1, min( k+3, i )
                  sum = h( j, k ) + v2*h( j, k+1 ) + v3*h( j, k+2 )
                  h( j, k ) = h( j, k ) - sum*t1
                  h( j, k+1 ) = h( j, k+1 ) - sum*t2
                  h( j, k+2 ) = h( j, k+2 ) - sum*t3
80             CONTINUE
!
               IF( wantz ) THEN
!
!                 Accumulate transformations in the matrix Z
!
                  DO 90 j = iloz, ihiz
                     sum = z( j, k ) + v2*z( j, k+1 ) + v3*z( j, k+2 )
                     z( j, k ) = z( j, k ) - sum*t1
                     z( j, k+1 ) = z( j, k+1 ) - sum*t2
                     z( j, k+2 ) = z( j, k+2 ) - sum*t3
90                CONTINUE
               END IF
            ELSE IF( nr.EQ.2 ) THEN
!
!              Apply G from the left to transform the rows of the matrix
!              in columns K to I2.
!
               DO 100 j = k, i2
                  sum = h( k, j ) + v2*h( k+1, j )
                  h( k, j ) = h( k, j ) - sum*t1
                  h( k+1, j ) = h( k+1, j ) - sum*t2
100            CONTINUE
!
!              Apply G from the right to transform the columns of the
!              matrix in rows I1 to min(K+3,I).
!
               DO 110 j = i1, i
                  sum = h( j, k ) + v2*h( j, k+1 )
                  h( j, k ) = h( j, k ) - sum*t1
                  h( j, k+1 ) = h( j, k+1 ) - sum*t2
110            CONTINUE
!
               IF( wantz ) THEN
!
!                 Accumulate transformations in the matrix Z
!
                  DO 120 j = iloz, ihiz
                     sum = z( j, k ) + v2*z( j, k+1 )
                     z( j, k ) = z( j, k ) - sum*t1
                     z( j, k+1 ) = z( j, k+1 ) - sum*t2
120               CONTINUE
               END IF
            END IF
130      CONTINUE
!
140   CONTINUE
!
!     Failure to converge in remaining number of iterations
!
      info = i
      RETURN
!
150   CONTINUE
!
      IF( l.EQ.i ) THEN
!
!        H(I,I-1) is negligible: one eigenvalue has converged.
!
         wr( i ) = h( i, i )
         wi( i ) = zero
      ELSE IF( l.EQ.i-1 ) THEN
!
!        H(I-1,I-2) is negligible: a pair of eigenvalues have converged.
!
!        Transform the 2-by-2 submatrix to standard Schur form,
!        and compute and store the eigenvalues.
!
         CALL dlanv2( h( i-1, i-1 ), h( i-1, i ), h( i, i-1 ), &
     &                h( i, i ), wr( i-1 ), wi( i-1 ), wr( i ), wi( i ), &
     &                cs, sn )
!
         IF( wantt ) THEN
!
!           Apply the transformation to the rest of H.
!
            IF( i2.GT.i ) &
     &         CALL drot( i2-i, h( i-1, i+1 ), ldh, h( i, i+1 ), ldh, &
     &                    cs, sn )
            CALL drot( i-i1-1, h( i1, i-1 ), 1, h( i1, i ), 1, cs, &
     &                 sn )
         END IF
         IF( wantz ) THEN
!
!           Apply the transformation to Z.
!
            CALL drot( nz, z( iloz, i-1 ), 1, z( iloz, i ), 1, cs, &
     &                 sn )
         END IF
      END IF
!     reset deflation counter
      kdefl = 0
!
!     return to start of the main loop with new value of I.
!
      i = l - 1
      GO TO 20
!
160   CONTINUE
      RETURN
!
!     End of DLAHQR
!

      END
!> \brief \b DLASCL multiplies a general rectangular matrix by a real scalar defined as cto/cfrom.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download DLASCL + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlascl.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlascl.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlascl.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          TYPE
!       INTEGER            INFO, KL, KU, LDA, M, N
!       DOUBLE PRECISION   CFROM, CTO
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLASCL multiplies the M by N real matrix A by the real scalar
!> CTO/CFROM.  This is done without over/underflow as long as the final
!> result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that
!> A may be full, upper triangular, lower triangular, upper Hessenberg,
!> or banded.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TYPE
!> \verbatim
!>          TYPE is CHARACTER*1
!>          TYPE indices the storage type of the input matrix.
!>          = 'G':  A is a full matrix.
!>          = 'L':  A is a lower triangular matrix.
!>          = 'U':  A is an upper triangular matrix.
!>          = 'H':  A is an upper Hessenberg matrix.
!>          = 'B':  A is a symmetric band matrix with lower bandwidth KL
!>                  and upper bandwidth KU and with the only the lower
!>                  half stored.
!>          = 'Q':  A is a symmetric band matrix with lower bandwidth KL
!>                  and upper bandwidth KU and with the only the upper
!>                  half stored.
!>          = 'Z':  A is a band matrix with lower bandwidth KL and upper
!>                  bandwidth KU. See DGBTRF for storage details.
!> \endverbatim
!>
!> \param[in] KL
!> \verbatim
!>          KL is INTEGER
!>          The lower bandwidth of A.  Referenced only if TYPE = 'B',
!>          'Q' or 'Z'.
!> \endverbatim
!>
!> \param[in] KU
!> \verbatim
!>          KU is INTEGER
!>          The upper bandwidth of A.  Referenced only if TYPE = 'B',
!>          'Q' or 'Z'.
!> \endverbatim
!>
!> \param[in] CFROM
!> \verbatim
!>          CFROM is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] CTO
!> \verbatim
!>          CTO is DOUBLE PRECISION
!>
!>          The matrix A is multiplied by CTO/CFROM. A(I,J) is computed
!>          without over/underflow if the final result CTO*A(I,J)/CFROM
!>          can be represented without over/underflow.  CFROM must be
!>          nonzero.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          The matrix to be multiplied by CTO/CFROM.  See TYPE for the
!>          storage type.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.
!>          If TYPE = 'G', 'L', 'U', 'H', LDA >= max(1,M);
!>             TYPE = 'B', LDA >= KL+1;
!>             TYPE = 'Q', LDA >= KU+1;
!>             TYPE = 'Z', LDA >= 2*KL+KU+1.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          0  - successful exit
!>          <0 - if INFO = -i, the i-th argument had an illegal value.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup lascl
!
!  =====================================================================


!> \brief \b DLAQR0 computes the eigenvalues of a Hessenberg matrix, and optionally the matrices from the Schur decomposition.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download DLAQR0 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaqr0.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaqr0.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaqr0.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAQR0( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI,
!                          ILOZ, IHIZ, Z, LDZ, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, LWORK, N
!       LOGICAL            WANTT, WANTZ
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   H( LDH, * ), WI( * ), WORK( * ), WR( * ),
!      $                   Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    DLAQR0 computes the eigenvalues of a Hessenberg matrix H
!>    and, optionally, the matrices T and Z from the Schur decomposition
!>    H = Z T Z**T, where T is an upper quasi-triangular matrix (the
!>    Schur form), and Z is the orthogonal matrix of Schur vectors.
!>
!>    Optionally Z may be postmultiplied into an input orthogonal
!>    matrix Q so that this routine can give the Schur factorization
!>    of a matrix A which has been reduced to the Hessenberg form H
!>    by the orthogonal matrix Q:  A = Q*H*Q**T = (QZ)*T*(QZ)**T.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] WANTT
!> \verbatim
!>          WANTT is LOGICAL
!>          = .TRUE. : the full Schur form T is required;
!>          = .FALSE.: only eigenvalues are required.
!> \endverbatim
!>
!> \param[in] WANTZ
!> \verbatim
!>          WANTZ is LOGICAL
!>          = .TRUE. : the matrix of Schur vectors Z is required;
!>          = .FALSE.: Schur vectors are not required.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           The order of the matrix H.  N >= 0.
!> \endverbatim
!>
!> \param[in] ILO
!> \verbatim
!>          ILO is INTEGER
!> \endverbatim
!>
!> \param[in] IHI
!> \verbatim
!>          IHI is INTEGER
!>           It is assumed that H is already upper triangular in rows
!>           and columns 1:ILO-1 and IHI+1:N and, if ILO > 1,
!>           H(ILO,ILO-1) is zero. ILO and IHI are normally set by a
!>           previous call to DGEBAL, and then passed to DGEHRD when the
!>           matrix output by DGEBAL is reduced to Hessenberg form.
!>           Otherwise, ILO and IHI should be set to 1 and N,
!>           respectively.  If N > 0, then 1 <= ILO <= IHI <= N.
!>           If N = 0, then ILO = 1 and IHI = 0.
!> \endverbatim
!>
!> \param[in,out] H
!> \verbatim
!>          H is DOUBLE PRECISION array, dimension (LDH,N)
!>           On entry, the upper Hessenberg matrix H.
!>           On exit, if INFO = 0 and WANTT is .TRUE., then H contains
!>           the upper quasi-triangular matrix T from the Schur
!>           decomposition (the Schur form); 2-by-2 diagonal blocks
!>           (corresponding to complex conjugate pairs of eigenvalues)
!>           are returned in standard form, with H(i,i) = H(i+1,i+1)
!>           and H(i+1,i)*H(i,i+1) < 0. If INFO = 0 and WANTT is
!>           .FALSE., then the contents of H are unspecified on exit.
!>           (The output value of H when INFO > 0 is given under the
!>           description of INFO below.)
!>
!>           This subroutine may explicitly set H(i,j) = 0 for i > j and
!>           j = 1, 2, ... ILO-1 or j = IHI+1, IHI+2, ... N.
!> \endverbatim
!>
!> \param[in] LDH
!> \verbatim
!>          LDH is INTEGER
!>           The leading dimension of the array H. LDH >= max(1,N).
!> \endverbatim
!>
!> \param[out] WR
!> \verbatim
!>          WR is DOUBLE PRECISION array, dimension (IHI)
!> \endverbatim
!>
!> \param[out] WI
!> \verbatim
!>          WI is DOUBLE PRECISION array, dimension (IHI)
!>           The real and imaginary parts, respectively, of the computed
!>           eigenvalues of H(ILO:IHI,ILO:IHI) are stored in WR(ILO:IHI)
!>           and WI(ILO:IHI). If two eigenvalues are computed as a
!>           complex conjugate pair, they are stored in consecutive
!>           elements of WR and WI, say the i-th and (i+1)th, with
!>           WI(i) > 0 and WI(i+1) < 0. If WANTT is .TRUE., then
!>           the eigenvalues are stored in the same order as on the
!>           diagonal of the Schur form returned in H, with
!>           WR(i) = H(i,i) and, if H(i:i+1,i:i+1) is a 2-by-2 diagonal
!>           block, WI(i) = sqrt(-H(i+1,i)*H(i,i+1)) and
!>           WI(i+1) = -WI(i).
!> \endverbatim
!>
!> \param[in] ILOZ
!> \verbatim
!>          ILOZ is INTEGER
!> \endverbatim
!>
!> \param[in] IHIZ
!> \verbatim
!>          IHIZ is INTEGER
!>           Specify the rows of Z to which transformations must be
!>           applied if WANTZ is .TRUE..
!>           1 <= ILOZ <= ILO; IHI <= IHIZ <= N.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is DOUBLE PRECISION array, dimension (LDZ,IHI)
!>           If WANTZ is .FALSE., then Z is not referenced.
!>           If WANTZ is .TRUE., then Z(ILO:IHI,ILOZ:IHIZ) is
!>           replaced by Z(ILO:IHI,ILOZ:IHIZ)*U where U is the
!>           orthogonal Schur factor of H(ILO:IHI,ILO:IHI).
!>           (The output value of Z when INFO > 0 is given under
!>           the description of INFO below.)
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>           The leading dimension of the array Z.  if WANTZ is .TRUE.
!>           then LDZ >= MAX(1,IHIZ).  Otherwise, LDZ >= 1.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension LWORK
!>           On exit, if LWORK = -1, WORK(1) returns an estimate of
!>           the optimal value for LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>           The dimension of the array WORK.  LWORK >= max(1,N)
!>           is sufficient, but LWORK typically as large as 6*N may
!>           be required for optimal performance.  A workspace query
!>           to determine the optimal workspace size is recommended.
!>
!>           If LWORK = -1, then DLAQR0 does a workspace query.
!>           In this case, DLAQR0 checks the input parameters and
!>           estimates the optimal workspace size for the given
!>           values of N, ILO and IHI.  The estimate is returned
!>           in WORK(1).  No error message related to LWORK is
!>           issued by XERBLA.  Neither H nor Z are accessed.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>             = 0:  successful exit
!>             > 0:  if INFO = i, DLAQR0 failed to compute all of
!>                the eigenvalues.  Elements 1:ilo-1 and i+1:n of WR
!>                and WI contain those eigenvalues which have been
!>                successfully computed.  (Failures are rare.)
!>
!>                If INFO > 0 and WANT is .FALSE., then on exit,
!>                the remaining unconverged eigenvalues are the eigen-
!>                values of the upper Hessenberg matrix rows and
!>                columns ILO through INFO of the final, output
!>                value of H.
!>
!>                If INFO > 0 and WANTT is .TRUE., then on exit
!>
!>           (*)  (initial value of H)*U  = U*(final value of H)
!>
!>                where U is an orthogonal matrix.  The final
!>                value of H is upper Hessenberg and quasi-triangular
!>                in rows and columns INFO+1 through IHI.
!>
!>                If INFO > 0 and WANTZ is .TRUE., then on exit
!>
!>                  (final value of Z(ILO:IHI,ILOZ:IHIZ)
!>                   =  (initial value of Z(ILO:IHI,ILOZ:IHIZ)*U
!>
!>                where U is the orthogonal matrix in (*) (regard-
!>                less of the value of WANTT.)
!>
!>                If INFO > 0 and WANTZ is .FALSE., then Z is not
!>                accessed.
!> \endverbatim
!
!> \par Contributors:
!  ==================
!>
!>       Karen Braman and Ralph Byers, Department of Mathematics,
!>       University of Kansas, USA
!
!> \par References:
!  ================
!>
!>       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
!>       Algorithm Part I: Maintaining Well Focused Shifts, and Level 3
!>       Performance, SIAM Journal of Matrix Analysis, volume 23, pages
!>       929--947, 2002.
!> \n
!>       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
!>       Algorithm Part II: Aggressive Early Deflation, SIAM Journal
!>       of Matrix Analysis, volume 23, pages 948--973, 2002.
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup laqr0
!
!  =====================================================================

      SUBROUTINE dlaqr0( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, &
     &                   ILOZ, IHIZ, Z, LDZ, WORK, LWORK, INFO )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, LWORK, N
      LOGICAL            WANTT, WANTZ
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   H( LDH, * ), WI( * ), WORK( * ), WR( * ), &
     &                   z( ldz, * )
!     ..
!
!  ================================================================
!
!     .. Parameters ..
!
!     ==== Matrices of order NTINY or smaller must be processed by
!     .    DLAHQR because of insufficient subdiagonal scratch space.
!     .    (This is a hard limit.) ====
      INTEGER            NTINY
      parameter( ntiny = 15 )
!
!     ==== Exceptional deflation windows:  try to cure rare
!     .    slow convergence by varying the size of the
!     .    deflation window after KEXNW iterations. ====
      INTEGER            KEXNW
      parameter( kexnw = 5 )
!
!     ==== Exceptional shifts: try to cure rare slow convergence
!     .    with ad-hoc exceptional shifts every KEXSH iterations.
!     .    ====
      INTEGER            KEXSH
      parameter( kexsh = 6 )
!
!     ==== The constants WILK1 and WILK2 are used to form the
!     .    exceptional shifts. ====
      DOUBLE PRECISION   WILK1, WILK2
      parameter( wilk1 = 0.75d0, wilk2 = -0.4375d0 )
      DOUBLE PRECISION   ZERO, ONE
      parameter( zero = 0.0d0, one = 1.0d0 )
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION   AA, BB, CC, CS, DD, SN, SS, SWAP
      INTEGER            I, INF, IT, ITMAX, K, KACC22, KBOT, KDU, KS, &
     &                   kt, ktop, ku, kv, kwh, kwtop, kwv, ld, ls, &
     &                   lwkopt, ndec, ndfl, nh, nho, nibble, nmin, ns, &
     &                   nsmax, nsr, nve, nw, nwmax, nwr, nwupbd
      LOGICAL            SORTED
      CHARACTER          JBCMPZ*2
!     ..
!     .. External Functions ..
!      INTEGER            ILAENV
!      EXTERNAL           ilaenv
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION   ZDUM( 1, 1 )
!     ..
!     .. External Subroutines ..
!      EXTERNAL           dlacpy, dlahqr, dlanv2, dlaqr3, dlaqr4, &
!     &                   dlaqr5
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          abs, dble, int, max, min, mod
!     ..
!     .. Executable Statements ..
      info = 0
!
!     ==== Quick return for N = 0: nothing to do. ====
!
      IF( n.EQ.0 ) THEN
         work( 1 ) = one
         RETURN
      END IF
!
      IF( n.LE.ntiny ) THEN
!
!        ==== Tiny matrices must use DLAHQR. ====
!
         lwkopt = 1
         IF( lwork.NE.-1 ) &
     &      CALL dlahqr( wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, &
     &                   iloz, ihiz, z, ldz, info )
      ELSE
!
!        ==== Use small bulge multi-shift QR with aggressive early
!        .    deflation on larger-than-tiny matrices. ====
!
!        ==== Hope for the best. ====
!
         info = 0
!
!        ==== Set up job flags for ILAENV. ====
!
         IF( wantt ) THEN
            jbcmpz( 1: 1 ) = 'S'
         ELSE
            jbcmpz( 1: 1 ) = 'E'
         END IF
         IF( wantz ) THEN
            jbcmpz( 2: 2 ) = 'V'
         ELSE
            jbcmpz( 2: 2 ) = 'N'
         END IF
!
!        ==== NWR = recommended deflation window size.  At this
!        .    point,  N .GT. NTINY = 15, so there is enough
!        .    subdiagonal workspace for NWR.GE.2 as required.
!        .    (In fact, there is enough subdiagonal space for
!        .    NWR.GE.4.) ====
!
         nwr = ilaenv( 13, 'DLAQR0', jbcmpz, n, ilo, ihi, lwork )
         nwr = max( 2, nwr )
         nwr = min( ihi-ilo+1, ( n-1 ) / 3, nwr )
!
!        ==== NSR = recommended number of simultaneous shifts.
!        .    At this point N .GT. NTINY = 15, so there is at
!        .    enough subdiagonal workspace for NSR to be even
!        .    and greater than or equal to two as required. ====
!
         nsr = ilaenv( 15, 'DLAQR0', jbcmpz, n, ilo, ihi, lwork )
         nsr = min( nsr, ( n-3 ) / 6, ihi-ilo )
         nsr = max( 2, nsr-mod( nsr, 2 ) )
!
!        ==== Estimate optimal workspace ====
!
!        ==== Workspace query call to DLAQR3 ====
!
         CALL dlaqr3( wantt, wantz, n, ilo, ihi, nwr+1, h, ldh, iloz, &
     &                ihiz, z, ldz, ls, ld, wr, wi, h, ldh, n, h, ldh, &
     &                n, h, ldh, work, -1 )
!
!        ==== Optimal workspace = MAX(DLAQR5, DLAQR3) ====
!
         lwkopt = max( 3*nsr / 2, int( work( 1 ) ) )
!
!        ==== Quick return in case of workspace query. ====
!
         IF( lwork.EQ.-1 ) THEN
            work( 1 ) = dble( lwkopt )
            RETURN
         END IF
!
!        ==== DLAHQR/DLAQR0 crossover point ====
!
         nmin = ilaenv( 12, 'DLAQR0', jbcmpz, n, ilo, ihi, lwork )
         nmin = max( ntiny, nmin )
!
!        ==== Nibble crossover point ====
!
         nibble = ilaenv( 14, 'DLAQR0', jbcmpz, n, ilo, ihi, lwork )
         nibble = max( 0, nibble )
!
!        ==== Accumulate reflections during ttswp?  Use block
!        .    2-by-2 structure during matrix-matrix multiply? ====
!
         kacc22 = ilaenv( 16, 'DLAQR0', jbcmpz, n, ilo, ihi, lwork )
         kacc22 = max( 0, kacc22 )
         kacc22 = min( 2, kacc22 )
!
!        ==== NWMAX = the largest possible deflation window for
!        .    which there is sufficient workspace. ====
!
         nwmax = min( ( n-1 ) / 3, lwork / 2 )
         nw = nwmax
!
!        ==== NSMAX = the Largest number of simultaneous shifts
!        .    for which there is sufficient workspace. ====
!
         nsmax = min( ( n-3 ) / 6, 2*lwork / 3 )
         nsmax = nsmax - mod( nsmax, 2 )
!
!        ==== NDFL: an iteration count restarted at deflation. ====
!
         ndfl = 1
!
!        ==== ITMAX = iteration limit ====
!
         itmax = max( 30, 2*kexsh )*max( 10, ( ihi-ilo+1 ) )
!
!        ==== Last row and column in the active block ====
!
         kbot = ihi
!
!        ==== Main Loop ====
!
         DO 80 it = 1, itmax
!
!           ==== Done when KBOT falls below ILO ====
!
            IF( kbot.LT.ilo ) &
     &         GO TO 90
!
!           ==== Locate active block ====
!
            DO 10 k = kbot, ilo + 1, -1
               IF( h( k, k-1 ).EQ.zero ) &
     &            GO TO 20
10          CONTINUE
            k = ilo
20          CONTINUE
            ktop = k
!
!           ==== Select deflation window size:
!           .    Typical Case:
!           .      If possible and advisable, nibble the entire
!           .      active block.  If not, use size MIN(NWR,NWMAX)
!           .      or MIN(NWR+1,NWMAX) depending upon which has
!           .      the smaller corresponding subdiagonal entry
!           .      (a heuristic).
!           .
!           .    Exceptional Case:
!           .      If there have been no deflations in KEXNW or
!           .      more iterations, then vary the deflation window
!           .      size.   At first, because, larger windows are,
!           .      in general, more powerful than smaller ones,
!           .      rapidly increase the window to the maximum possible.
!           .      Then, gradually reduce the window size. ====
!
            nh = kbot - ktop + 1
            nwupbd = min( nh, nwmax )
            IF( ndfl.LT.kexnw ) THEN
               nw = min( nwupbd, nwr )
            ELSE
               nw = min( nwupbd, 2*nw )
            END IF
            IF( nw.LT.nwmax ) THEN
               IF( nw.GE.nh-1 ) THEN
                  nw = nh
               ELSE
                  kwtop = kbot - nw + 1
                  IF( abs( h( kwtop, kwtop-1 ) ).GT. &
     &                abs( h( kwtop-1, kwtop-2 ) ) )nw = nw + 1
               END IF
            END IF
            IF( ndfl.LT.kexnw ) THEN
               ndec = -1
            ELSE IF( ndec.GE.0 .OR. nw.GE.nwupbd ) THEN
               ndec = ndec + 1
               IF( nw-ndec.LT.2 ) &
     &            ndec = 0
               nw = nw - ndec
            END IF
!
!           ==== Aggressive early deflation:
!           .    split workspace under the subdiagonal into
!           .      - an nw-by-nw work array V in the lower
!           .        left-hand-corner,
!           .      - an NW-by-at-least-NW-but-more-is-better
!           .        (NW-by-NHO) horizontal work array along
!           .        the bottom edge,
!           .      - an at-least-NW-but-more-is-better (NHV-by-NW)
!           .        vertical work array along the left-hand-edge.
!           .        ====
!
            kv = n - nw + 1
            kt = nw + 1
            nho = ( n-nw-1 ) - kt + 1
            kwv = nw + 2
            nve = ( n-nw ) - kwv + 1
!
!           ==== Aggressive early deflation ====
!
            CALL dlaqr3( wantt, wantz, n, ktop, kbot, nw, h, ldh, &
     &                   iloz, &
     &                   ihiz, z, ldz, ls, ld, wr, wi, h( kv, 1 ), ldh, &
     &                   nho, h( kv, kt ), ldh, nve, h( kwv, 1 ), ldh, &
     &                   work, lwork )
!
!           ==== Adjust KBOT accounting for new deflations. ====
!
            kbot = kbot - ld
!
!           ==== KS points to the shifts. ====
!
            ks = kbot - ls + 1
!
!           ==== Skip an expensive QR sweep if there is a (partly
!           .    heuristic) reason to expect that many eigenvalues
!           .    will deflate without it.  Here, the QR sweep is
!           .    skipped if many eigenvalues have just been deflated
!           .    or if the remaining active block is small.
!
            IF( ( ld.EQ.0 ) .OR. ( ( 100*ld.LE.nw*nibble ) .AND. ( kbot- &
     &          ktop+1.GT.min( nmin, nwmax ) ) ) ) THEN
!
!              ==== NS = nominal number of simultaneous shifts.
!              .    This may be lowered (slightly) if DLAQR3
!              .    did not provide that many shifts. ====
!
               ns = min( nsmax, nsr, max( 2, kbot-ktop ) )
               ns = ns - mod( ns, 2 )
!
!              ==== If there have been no deflations
!              .    in a multiple of KEXSH iterations,
!              .    then try exceptional shifts.
!              .    Otherwise use shifts provided by
!              .    DLAQR3 above or from the eigenvalues
!              .    of a trailing principal submatrix. ====
!
               IF( mod( ndfl, kexsh ).EQ.0 ) THEN
                  ks = kbot - ns + 1
                  DO 30 i = kbot, max( ks+1, ktop+2 ), -2
                     ss = abs( h( i, i-1 ) ) + abs( h( i-1, i-2 ) )
                     aa = wilk1*ss + h( i, i )
                     bb = ss
                     cc = wilk2*ss
                     dd = aa
                     CALL dlanv2( aa, bb, cc, dd, wr( i-1 ), &
     &                            wi( i-1 ), &
     &                            wr( i ), wi( i ), cs, sn )
30                CONTINUE
                  IF( ks.EQ.ktop ) THEN
                     wr( ks+1 ) = h( ks+1, ks+1 )
                     wi( ks+1 ) = zero
                     wr( ks ) = wr( ks+1 )
                     wi( ks ) = wi( ks+1 )
                  END IF
               ELSE
!
!                 ==== Got NS/2 or fewer shifts? Use DLAQR4 or
!                 .    DLAHQR on a trailing principal submatrix to
!                 .    get more. (Since NS.LE.NSMAX.LE.(N-3)/6,
!                 .    there is enough space below the subdiagonal
!                 .    to fit an NS-by-NS scratch array.) ====
!
                  IF( kbot-ks+1.LE.ns / 2 ) THEN
                     ks = kbot - ns + 1
                     kt = n - ns + 1
                     CALL dlacpy( 'A', ns, ns, h( ks, ks ), ldh, &
     &                            h( kt, 1 ), ldh )
                     IF( ns.GT.nmin ) THEN
                        CALL dlaqr4( .false., .false., ns, 1, ns, &
     &                               h( kt, 1 ), ldh, wr( ks ), &
     &                               wi( ks ), 1, 1, zdum, 1, work, &
     &                               lwork, inf )
                     ELSE
                        CALL dlahqr( .false., .false., ns, 1, ns, &
     &                               h( kt, 1 ), ldh, wr( ks ), &
     &                               wi( ks ), 1, 1, zdum, 1, inf )
                     END IF
                     ks = ks + inf
!
!                    ==== In case of a rare QR failure use
!                    .    eigenvalues of the trailing 2-by-2
!                    .    principal submatrix.  ====
!
                     IF( ks.GE.kbot ) THEN
                        aa = h( kbot-1, kbot-1 )
                        cc = h( kbot, kbot-1 )
                        bb = h( kbot-1, kbot )
                        dd = h( kbot, kbot )
                        CALL dlanv2( aa, bb, cc, dd, wr( kbot-1 ), &
     &                               wi( kbot-1 ), wr( kbot ), &
     &                               wi( kbot ), cs, sn )
                        ks = kbot - 1
                     END IF
                  END IF
!
                  IF( kbot-ks+1.GT.ns ) THEN
!
!                    ==== Sort the shifts (Helps a little)
!                    .    Bubble sort keeps complex conjugate
!                    .    pairs together. ====
!
                     sorted = .false.
                     DO 50 k = kbot, ks + 1, -1
                        IF( sorted ) &
     &                     GO TO 60
                        sorted = .true.
                        DO 40 i = ks, k - 1
                           IF( abs( wr( i ) )+abs( wi( i ) ).LT. &
     &                         abs( wr( i+1 ) )+abs( wi( i+1 ) ) ) THEN
                              sorted = .false.
!
                              swap = wr( i )
                              wr( i ) = wr( i+1 )
                              wr( i+1 ) = swap
!
                              swap = wi( i )
                              wi( i ) = wi( i+1 )
                              wi( i+1 ) = swap
                           END IF
40                      CONTINUE
50                   CONTINUE
60                   CONTINUE
                  END IF
!
!                 ==== Shuffle shifts into pairs of real shifts
!                 .    and pairs of complex conjugate shifts
!                 .    assuming complex conjugate shifts are
!                 .    already adjacent to one another. (Yes,
!                 .    they are.)  ====
!
                  DO 70 i = kbot, ks + 2, -2
                     IF( wi( i ).NE.-wi( i-1 ) ) THEN
!
                        swap = wr( i )
                        wr( i ) = wr( i-1 )
                        wr( i-1 ) = wr( i-2 )
                        wr( i-2 ) = swap
!
                        swap = wi( i )
                        wi( i ) = wi( i-1 )
                        wi( i-1 ) = wi( i-2 )
                        wi( i-2 ) = swap
                     END IF
70                CONTINUE
               END IF
!
!              ==== If there are only two shifts and both are
!              .    real, then use only one.  ====
!
               IF( kbot-ks+1.EQ.2 ) THEN
                  IF( wi( kbot ).EQ.zero ) THEN
                     IF( abs( wr( kbot )-h( kbot, kbot ) ).LT. &
     &                   abs( wr( kbot-1 )-h( kbot, kbot ) ) ) THEN
                        wr( kbot-1 ) = wr( kbot )
                     ELSE
                        wr( kbot ) = wr( kbot-1 )
                     END IF
                  END IF
               END IF
!
!              ==== Use up to NS of the the smallest magnitude
!              .    shifts.  If there aren't NS shifts available,
!              .    then use them all, possibly dropping one to
!              .    make the number of shifts even. ====
!
               ns = min( ns, kbot-ks+1 )
               ns = ns - mod( ns, 2 )
               ks = kbot - ns + 1
!
!              ==== Small-bulge multi-shift QR sweep:
!              .    split workspace under the subdiagonal into
!              .    - a KDU-by-KDU work array U in the lower
!              .      left-hand-corner,
!              .    - a KDU-by-at-least-KDU-but-more-is-better
!              .      (KDU-by-NHo) horizontal work array WH along
!              .      the bottom edge,
!              .    - and an at-least-KDU-but-more-is-better-by-KDU
!              .      (NVE-by-KDU) vertical work WV arrow along
!              .      the left-hand-edge. ====
!
               kdu = 2*ns
               ku = n - kdu + 1
               kwh = kdu + 1
               nho = ( n-kdu+1-4 ) - ( kdu+1 ) + 1
               kwv = kdu + 4
               nve = n - kdu - kwv + 1
!
!              ==== Small-bulge multi-shift QR sweep ====
!
               CALL dlaqr5( wantt, wantz, kacc22, n, ktop, kbot, ns, &
     &                      wr( ks ), wi( ks ), h, ldh, iloz, ihiz, z, &
     &                      ldz, work, 3, h( ku, 1 ), ldh, nve, &
     &                      h( kwv, 1 ), ldh, nho, h( ku, kwh ), ldh )
            END IF
!
!           ==== Note progress (or the lack of it). ====
!
            IF( ld.GT.0 ) THEN
               ndfl = 1
            ELSE
               ndfl = ndfl + 1
            END IF
!
!           ==== End of main loop ====
80       CONTINUE
!
!        ==== Iteration limit exceeded.  Set INFO to show where
!        .    the problem occurred and exit. ====
!
         info = kbot
90       CONTINUE
      END IF
!
!     ==== Return the optimal value of LWORK. ====
!
      work( 1 ) = dble( lwkopt )
!
!     ==== End of DLAQR0 ====
!

      END
!> \brief \b DLALN2 solves a 1-by-1 or 2-by-2 linear system of equations of the specified form.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download DLALN2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaln2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaln2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaln2.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLALN2( LTRANS, NA, NW, SMIN, CA, A, LDA, D1, D2, B,
!                          LDB, WR, WI, X, LDX, SCALE, XNORM, INFO )
!
!       .. Scalar Arguments ..
!       LOGICAL            LTRANS
!       INTEGER            INFO, LDA, LDB, LDX, NA, NW
!       DOUBLE PRECISION   CA, D1, D2, SCALE, SMIN, WI, WR, XNORM
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLALN2 solves a system of the form  (ca A - w D ) X = s B
!> or (ca A**T - w D) X = s B   with possible scaling ("s") and
!> perturbation of A.  (A**T means A-transpose.)
!>
!> A is an NA x NA real matrix, ca is a real scalar, D is an NA x NA
!> real diagonal matrix, w is a real or complex value, and X and B are
!> NA x 1 matrices -- real if w is real, complex if w is complex.  NA
!> may be 1 or 2.
!>
!> If w is complex, X and B are represented as NA x 2 matrices,
!> the first column of each being the real part and the second
!> being the imaginary part.
!>
!> "s" is a scaling factor (<= 1), computed by DLALN2, which is
!> so chosen that X can be computed without overflow.  X is further
!> scaled if necessary to assure that norm(ca A - w D)*norm(X) is less
!> than overflow.
!>
!> If both singular values of (ca A - w D) are less than SMIN,
!> SMIN*identity will be used instead of (ca A - w D).  If only one
!> singular value is less than SMIN, one element of (ca A - w D) will be
!> perturbed enough to make the smallest singular value roughly SMIN.
!> If both singular values are at least SMIN, (ca A - w D) will not be
!> perturbed.  In any case, the perturbation will be at most some small
!> multiple of max( SMIN, ulp*norm(ca A - w D) ).  The singular values
!> are computed by infinity-norm approximations, and thus will only be
!> correct to a factor of 2 or so.
!>
!> Note: all input quantities are assumed to be smaller than overflow
!> by a reasonable factor.  (See BIGNUM.)
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] LTRANS
!> \verbatim
!>          LTRANS is LOGICAL
!>          =.TRUE.:  A-transpose will be used.
!>          =.FALSE.: A will be used (not transposed.)
!> \endverbatim
!>
!> \param[in] NA
!> \verbatim
!>          NA is INTEGER
!>          The size of the matrix A.  It may (only) be 1 or 2.
!> \endverbatim
!>
!> \param[in] NW
!> \verbatim
!>          NW is INTEGER
!>          1 if "w" is real, 2 if "w" is complex.  It may only be 1
!>          or 2.
!> \endverbatim
!>
!> \param[in] SMIN
!> \verbatim
!>          SMIN is DOUBLE PRECISION
!>          The desired lower bound on the singular values of A.  This
!>          should be a safe distance away from underflow or overflow,
!>          say, between (underflow/machine precision) and  (machine
!>          precision * overflow ).  (See BIGNUM and ULP.)
!> \endverbatim
!>
!> \param[in] CA
!> \verbatim
!>          CA is DOUBLE PRECISION
!>          The coefficient c, which A is multiplied by.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,NA)
!>          The NA x NA matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of A.  It must be at least NA.
!> \endverbatim
!>
!> \param[in] D1
!> \verbatim
!>          D1 is DOUBLE PRECISION
!>          The 1,1 element in the diagonal matrix D.
!> \endverbatim
!>
!> \param[in] D2
!> \verbatim
!>          D2 is DOUBLE PRECISION
!>          The 2,2 element in the diagonal matrix D.  Not used if NA=1.
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (LDB,NW)
!>          The NA x NW matrix B (right-hand side).  If NW=2 ("w" is
!>          complex), column 1 contains the real part of B and column 2
!>          contains the imaginary part.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of B.  It must be at least NA.
!> \endverbatim
!>
!> \param[in] WR
!> \verbatim
!>          WR is DOUBLE PRECISION
!>          The real part of the scalar "w".
!> \endverbatim
!>
!> \param[in] WI
!> \verbatim
!>          WI is DOUBLE PRECISION
!>          The imaginary part of the scalar "w".  Not used if NW=1.
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension (LDX,NW)
!>          The NA x NW matrix X (unknowns), as computed by DLALN2.
!>          If NW=2 ("w" is complex), on exit, column 1 will contain
!>          the real part of X and column 2 will contain the imaginary
!>          part.
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of X.  It must be at least NA.
!> \endverbatim
!>
!> \param[out] SCALE
!> \verbatim
!>          SCALE is DOUBLE PRECISION
!>          The scale factor that B must be multiplied by to insure
!>          that overflow does not occur when computing X.  Thus,
!>          (ca A - w D) X  will be SCALE*B, not B (ignoring
!>          perturbations of A.)  It will be at most 1.
!> \endverbatim
!>
!> \param[out] XNORM
!> \verbatim
!>          XNORM is DOUBLE PRECISION
!>          The infinity-norm of X, when X is regarded as an NA x NW
!>          real matrix.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          An error flag.  It will be set to zero if no error occurs,
!>          a negative number if an argument is in error, or a positive
!>          number if  ca A - w D  had to be perturbed.
!>          The possible values are:
!>          = 0: No error occurred, and (ca A - w D) did not have to be
!>                 perturbed.
!>          = 1: (ca A - w D) had to be perturbed to make its smallest
!>               (or only) singular value greater than SMIN.
!>          NOTE: In the interests of speed, this routine does not
!>                check the inputs for errors.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup laln2
!
!  =====================================================================

      SUBROUTINE dlaln2( LTRANS, NA, NW, SMIN, CA, A, LDA, D1, D2, B, &
     &                   LDB, WR, WI, X, LDX, SCALE, XNORM, INFO )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      LOGICAL            LTRANS
      INTEGER            INFO, LDA, LDB, LDX, NA, NW
      DOUBLE PRECISION   CA, D1, D2, SCALE, SMIN, WI, WR, XNORM
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), X( LDX, * )
!     ..
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      parameter( zero = 0.0d0, one = 1.0d0 )
      DOUBLE PRECISION   TWO
      parameter( two = 2.0d0 )
!     ..
!     .. Local Scalars ..
      INTEGER            ICMAX, J
      DOUBLE PRECISION   BBND, BI1, BI2, BIGNUM, BNORM, BR1, BR2, CI21, &
     &                   ci22, cmax, cnorm, cr21, cr22, csi, csr, li21, &
     &                   lr21, smini, smlnum, temp, u22abs, ui11, ui11r, &
     &                   ui12, ui12s, ui22, ur11, ur11r, ur12, ur12s, &
     &                   ur22, xi1, xi2, xr1, xr2
!     ..
!     .. Local Arrays ..
      LOGICAL            RSWAP( 4 ), ZSWAP( 4 )
      INTEGER            IPIVOT( 4, 4 )
      DOUBLE PRECISION   CI( 2, 2 ), CIV( 4 ), CR( 2, 2 ), CRV( 4 )
!     ..
!     .. External Functions ..
!      DOUBLE PRECISION   DLAMCH
!      EXTERNAL           dlamch
!     ..
!     .. External Subroutines ..
!      EXTERNAL           dladiv
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          abs, max
!     ..
!     .. Equivalences ..
      equivalence( ci( 1, 1 ), civ( 1 ) ), &
     &                   ( cr( 1, 1 ), crv( 1 ) )
!     ..
!     .. Data statements ..
      DATA               zswap / .false., .false., .true., .true. /
      DATA               rswap / .false., .true., .false., .true. /
      DATA               ipivot / 1, 2, 3, 4, 2, 1, 4, 3, 3, 4, 1, 2, 4, &
     &                   3, 2, 1 /
!     ..
!     .. Executable Statements ..
!
!     Compute BIGNUM
!
      smlnum = two*dlamch( 'Safe minimum' )
      bignum = one / smlnum
      smini = max( smin, smlnum )
!
!     Don't check for input errors
!
      info = 0
!
!     Standard Initializations
!
      scale = one
!
      IF( na.EQ.1 ) THEN
!
!        1 x 1  (i.e., scalar) system   C X = B
!
         IF( nw.EQ.1 ) THEN
!
!           Real 1x1 system.
!
!           C = ca A - w D
!
            csr = ca*a( 1, 1 ) - wr*d1
            cnorm = abs( csr )
!
!           If | C | < SMINI, use C = SMINI
!
            IF( cnorm.LT.smini ) THEN
               csr = smini
               cnorm = smini
               info = 1
            END IF
!
!           Check scaling for  X = B / C
!
            bnorm = abs( b( 1, 1 ) )
            IF( cnorm.LT.one .AND. bnorm.GT.one ) THEN
               IF( bnorm.GT.bignum*cnorm ) &
     &            scale = one / bnorm
            END IF
!
!           Compute X
!
            x( 1, 1 ) = ( b( 1, 1 )*scale ) / csr
            xnorm = abs( x( 1, 1 ) )
         ELSE
!
!           Complex 1x1 system (w is complex)
!
!           C = ca A - w D
!
            csr = ca*a( 1, 1 ) - wr*d1
            csi = -wi*d1
            cnorm = abs( csr ) + abs( csi )
!
!           If | C | < SMINI, use C = SMINI
!
            IF( cnorm.LT.smini ) THEN
               csr = smini
               csi = zero
               cnorm = smini
               info = 1
            END IF
!
!           Check scaling for  X = B / C
!
            bnorm = abs( b( 1, 1 ) ) + abs( b( 1, 2 ) )
            IF( cnorm.LT.one .AND. bnorm.GT.one ) THEN
               IF( bnorm.GT.bignum*cnorm ) &
     &            scale = one / bnorm
            END IF
!
!           Compute X
!
            CALL dladiv( scale*b( 1, 1 ), scale*b( 1, 2 ), csr, csi, &
     &                   x( 1, 1 ), x( 1, 2 ) )
            xnorm = abs( x( 1, 1 ) ) + abs( x( 1, 2 ) )
         END IF
!
      ELSE
!
!        2x2 System
!
!        Compute the real part of  C = ca A - w D  (or  ca A**T - w D )
!
         cr( 1, 1 ) = ca*a( 1, 1 ) - wr*d1
         cr( 2, 2 ) = ca*a( 2, 2 ) - wr*d2
         IF( ltrans ) THEN
            cr( 1, 2 ) = ca*a( 2, 1 )
            cr( 2, 1 ) = ca*a( 1, 2 )
         ELSE
            cr( 2, 1 ) = ca*a( 2, 1 )
            cr( 1, 2 ) = ca*a( 1, 2 )
         END IF
!
         IF( nw.EQ.1 ) THEN
!
!           Real 2x2 system  (w is real)
!
!           Find the largest element in C
!
            cmax = zero
            icmax = 0
!
            DO 10 j = 1, 4
               IF( abs( crv( j ) ).GT.cmax ) THEN
                  cmax = abs( crv( j ) )
                  icmax = j
               END IF
10          CONTINUE
!
!           If norm(C) < SMINI, use SMINI*identity.
!
            IF( cmax.LT.smini ) THEN
               bnorm = max( abs( b( 1, 1 ) ), abs( b( 2, 1 ) ) )
               IF( smini.LT.one .AND. bnorm.GT.one ) THEN
                  IF( bnorm.GT.bignum*smini ) &
     &               scale = one / bnorm
               END IF
               temp = scale / smini
               x( 1, 1 ) = temp*b( 1, 1 )
               x( 2, 1 ) = temp*b( 2, 1 )
               xnorm = temp*bnorm
               info = 1
               RETURN
            END IF
!
!           Gaussian elimination with complete pivoting.
!
            ur11 = crv( icmax )
            cr21 = crv( ipivot( 2, icmax ) )
            ur12 = crv( ipivot( 3, icmax ) )
            cr22 = crv( ipivot( 4, icmax ) )
            ur11r = one / ur11
            lr21 = ur11r*cr21
            ur22 = cr22 - ur12*lr21
!
!           If smaller pivot < SMINI, use SMINI
!
            IF( abs( ur22 ).LT.smini ) THEN
               ur22 = smini
               info = 1
            END IF
            IF( rswap( icmax ) ) THEN
               br1 = b( 2, 1 )
               br2 = b( 1, 1 )
            ELSE
               br1 = b( 1, 1 )
               br2 = b( 2, 1 )
            END IF
            br2 = br2 - lr21*br1
            bbnd = max( abs( br1*( ur22*ur11r ) ), abs( br2 ) )
            IF( bbnd.GT.one .AND. abs( ur22 ).LT.one ) THEN
               IF( bbnd.GE.bignum*abs( ur22 ) ) &
     &            scale = one / bbnd
            END IF
!
            xr2 = ( br2*scale ) / ur22
            xr1 = ( scale*br1 )*ur11r - xr2*( ur11r*ur12 )
            IF( zswap( icmax ) ) THEN
               x( 1, 1 ) = xr2
               x( 2, 1 ) = xr1
            ELSE
               x( 1, 1 ) = xr1
               x( 2, 1 ) = xr2
            END IF
            xnorm = max( abs( xr1 ), abs( xr2 ) )
!
!           Further scaling if  norm(A) norm(X) > overflow
!
            IF( xnorm.GT.one .AND. cmax.GT.one ) THEN
               IF( xnorm.GT.bignum / cmax ) THEN
                  temp = cmax / bignum
                  x( 1, 1 ) = temp*x( 1, 1 )
                  x( 2, 1 ) = temp*x( 2, 1 )
                  xnorm = temp*xnorm
                  scale = temp*scale
               END IF
            END IF
         ELSE
!
!           Complex 2x2 system  (w is complex)
!
!           Find the largest element in C
!
            ci( 1, 1 ) = -wi*d1
            ci( 2, 1 ) = zero
            ci( 1, 2 ) = zero
            ci( 2, 2 ) = -wi*d2
            cmax = zero
            icmax = 0
!
            DO 20 j = 1, 4
               IF( abs( crv( j ) )+abs( civ( j ) ).GT.cmax ) THEN
                  cmax = abs( crv( j ) ) + abs( civ( j ) )
                  icmax = j
               END IF
20          CONTINUE
!
!           If norm(C) < SMINI, use SMINI*identity.
!
            IF( cmax.LT.smini ) THEN
               bnorm = max( abs( b( 1, 1 ) )+abs( b( 1, 2 ) ), &
     &                 abs( b( 2, 1 ) )+abs( b( 2, 2 ) ) )
               IF( smini.LT.one .AND. bnorm.GT.one ) THEN
                  IF( bnorm.GT.bignum*smini ) &
     &               scale = one / bnorm
               END IF
               temp = scale / smini
               x( 1, 1 ) = temp*b( 1, 1 )
               x( 2, 1 ) = temp*b( 2, 1 )
               x( 1, 2 ) = temp*b( 1, 2 )
               x( 2, 2 ) = temp*b( 2, 2 )
               xnorm = temp*bnorm
               info = 1
               RETURN
            END IF
!
!           Gaussian elimination with complete pivoting.
!
            ur11 = crv( icmax )
            ui11 = civ( icmax )
            cr21 = crv( ipivot( 2, icmax ) )
            ci21 = civ( ipivot( 2, icmax ) )
            ur12 = crv( ipivot( 3, icmax ) )
            ui12 = civ( ipivot( 3, icmax ) )
            cr22 = crv( ipivot( 4, icmax ) )
            ci22 = civ( ipivot( 4, icmax ) )
            IF( icmax.EQ.1 .OR. icmax.EQ.4 ) THEN
!
!              Code when off-diagonals of pivoted C are real
!
               IF( abs( ur11 ).GT.abs( ui11 ) ) THEN
                  temp = ui11 / ur11
                  ur11r = one / ( ur11*( one+temp**2 ) )
                  ui11r = -temp*ur11r
               ELSE
                  temp = ur11 / ui11
                  ui11r = -one / ( ui11*( one+temp**2 ) )
                  ur11r = -temp*ui11r
               END IF
               lr21 = cr21*ur11r
               li21 = cr21*ui11r
               ur12s = ur12*ur11r
               ui12s = ur12*ui11r
               ur22 = cr22 - ur12*lr21
               ui22 = ci22 - ur12*li21
            ELSE
!
!              Code when diagonals of pivoted C are real
!
               ur11r = one / ur11
               ui11r = zero
               lr21 = cr21*ur11r
               li21 = ci21*ur11r
               ur12s = ur12*ur11r
               ui12s = ui12*ur11r
               ur22 = cr22 - ur12*lr21 + ui12*li21
               ui22 = -ur12*li21 - ui12*lr21
            END IF
            u22abs = abs( ur22 ) + abs( ui22 )
!
!           If smaller pivot < SMINI, use SMINI
!
            IF( u22abs.LT.smini ) THEN
               ur22 = smini
               ui22 = zero
               info = 1
            END IF
            IF( rswap( icmax ) ) THEN
               br2 = b( 1, 1 )
               br1 = b( 2, 1 )
               bi2 = b( 1, 2 )
               bi1 = b( 2, 2 )
            ELSE
               br1 = b( 1, 1 )
               br2 = b( 2, 1 )
               bi1 = b( 1, 2 )
               bi2 = b( 2, 2 )
            END IF
            br2 = br2 - lr21*br1 + li21*bi1
            bi2 = bi2 - li21*br1 - lr21*bi1
            bbnd = max( ( abs( br1 )+abs( bi1 ) )* &
     &             ( u22abs*( abs( ur11r )+abs( ui11r ) ) ), &
     &             abs( br2 )+abs( bi2 ) )
            IF( bbnd.GT.one .AND. u22abs.LT.one ) THEN
               IF( bbnd.GE.bignum*u22abs ) THEN
                  scale = one / bbnd
                  br1 = scale*br1
                  bi1 = scale*bi1
                  br2 = scale*br2
                  bi2 = scale*bi2
               END IF
            END IF
!
            CALL dladiv( br2, bi2, ur22, ui22, xr2, xi2 )
            xr1 = ur11r*br1 - ui11r*bi1 - ur12s*xr2 + ui12s*xi2
            xi1 = ui11r*br1 + ur11r*bi1 - ui12s*xr2 - ur12s*xi2
            IF( zswap( icmax ) ) THEN
               x( 1, 1 ) = xr2
               x( 2, 1 ) = xr1
               x( 1, 2 ) = xi2
               x( 2, 2 ) = xi1
            ELSE
               x( 1, 1 ) = xr1
               x( 2, 1 ) = xr2
               x( 1, 2 ) = xi1
               x( 2, 2 ) = xi2
            END IF
            xnorm = max( abs( xr1 )+abs( xi1 ), abs( xr2 )+abs( xi2 ) )
!
!           Further scaling if  norm(A) norm(X) > overflow
!
            IF( xnorm.GT.one .AND. cmax.GT.one ) THEN
               IF( xnorm.GT.bignum / cmax ) THEN
                  temp = cmax / bignum
                  x( 1, 1 ) = temp*x( 1, 1 )
                  x( 2, 1 ) = temp*x( 2, 1 )
                  x( 1, 2 ) = temp*x( 1, 2 )
                  x( 2, 2 ) = temp*x( 2, 2 )
                  xnorm = temp*xnorm
                  scale = temp*scale
               END IF
            END IF
         END IF
      END IF
!
      RETURN
!
!     End of DLALN2
!

      END
!> \brief \b DLARF applies an elementary reflector to a general rectangular matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download DLARF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarf.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarf.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarf.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          SIDE
!       INTEGER            INCV, LDC, M, N
!       DOUBLE PRECISION   TAU
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   C( LDC, * ), V( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLARF applies a real elementary reflector H to a real m by n matrix
!> C, from either the left or the right. H is represented in the form
!>
!>       H = I - tau * v * v**T
!>
!> where tau is a real scalar and v is a real vector.
!>
!> If tau = 0, then H is taken to be the unit matrix.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'L': form  H * C
!>          = 'R': form  C * H
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix C.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix C.
!> \endverbatim
!>
!> \param[in] V
!> \verbatim
!>          V is DOUBLE PRECISION array, dimension
!>                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
!>                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
!>          The vector v in the representation of H. V is not used if
!>          TAU = 0.
!> \endverbatim
!>
!> \param[in] INCV
!> \verbatim
!>          INCV is INTEGER
!>          The increment between elements of v. INCV <> 0.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION
!>          The value tau in the representation of H.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension (LDC,N)
!>          On entry, the m by n matrix C.
!>          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
!>          or C * H if SIDE = 'R'.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C. LDC >= max(1,M).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension
!>                         (N) if SIDE = 'L'
!>                      or (M) if SIDE = 'R'
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup larf
!
!  =====================================================================

      SUBROUTINE dlarf( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            INCV, LDC, M, N
      DOUBLE PRECISION   TAU
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   C( LDC, * ), V( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      parameter( one = 1.0d+0, zero = 0.0d+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            APPLYLEFT
      INTEGER            I, LASTV, LASTC
!     ..
!     .. External Subroutines ..
!      EXTERNAL           dgemv, dger
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      INTEGER            ILADLR, ILADLC
!      EXTERNAL           lsame, iladlr, iladlc
!     ..
!     .. Executable Statements ..
!
      applyleft = lsame( side, 'L' )
      lastv = 0
      lastc = 0
      IF( tau.NE.zero ) THEN
!     Set up variables for scanning V.  LASTV begins pointing to the end
!     of V.
         IF( applyleft ) THEN
            lastv = m
         ELSE
            lastv = n
         END IF
         IF( incv.GT.0 ) THEN
            i = 1 + (lastv-1) * incv
         ELSE
            i = 1
         END IF
!     Look for the last non-zero row in V.
         DO WHILE( lastv.GT.0 .AND. v( i ).EQ.zero )
            lastv = lastv - 1
            i = i - incv
         END DO
         IF( applyleft ) THEN
!     Scan for the last non-zero column in C(1:lastv,:).
            lastc = iladlc(lastv, n, c, ldc)
         ELSE
!     Scan for the last non-zero row in C(:,1:lastv).
            lastc = iladlr(m, lastv, c, ldc)
         END IF
      END IF
!     Note that lastc.eq.0 renders the BLAS operations null; no special
!     case is needed at this level.
      IF( applyleft ) THEN
!
!        Form  H * C
!
         IF( lastv.GT.0 ) THEN
!
!           w(1:lastc,1) := C(1:lastv,1:lastc)**T * v(1:lastv,1)
!
            CALL dgemv( 'Transpose', lastv, lastc, one, c, ldc, v, &
     &                  incv, &
     &           zero, work, 1 )
!
!           C(1:lastv,1:lastc) := C(...) - v(1:lastv,1) * w(1:lastc,1)**T
!
            CALL dger( lastv, lastc, -tau, v, incv, work, 1, c, ldc )
         END IF
      ELSE
!
!        Form  C * H
!
         IF( lastv.GT.0 ) THEN
!
!           w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1)
!
            CALL dgemv( 'No transpose', lastc, lastv, one, c, ldc, &
     &           v, incv, zero, work, 1 )
!
!           C(1:lastc,1:lastv) := C(...) - w(1:lastc,1) * v(1:lastv,1)**T
!
            CALL dger( lastc, lastv, -tau, work, 1, v, incv, c, ldc )
         END IF
      END IF
      RETURN
!
!     End of DLARF
!

      END

!> \brief \b DLAPY2 returns sqrt(x2+y2).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download DLAPY2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlapy2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlapy2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlapy2.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION DLAPY2( X, Y )
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION   X, Y
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
!> overflow and unnecessary underflow.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] X
!> \verbatim
!>          X is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] Y
!> \verbatim
!>          Y is DOUBLE PRECISION
!>          X and Y specify the values x and y.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup lapy2
!
!  =====================================================================

      DOUBLE PRECISION FUNCTION dlapy2( X, Y )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   x, y
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   zero
      parameter( zero = 0.0d0 )
      DOUBLE PRECISION   one
      parameter( one = 1.0d0 )
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION   w, xabs, yabs, z, hugeval
      LOGICAL            x_is_nan, y_is_nan
!     ..
!     .. External Functions ..
!      LOGICAL            disnan
!      EXTERNAL           disnan
!     ..
!     .. External Subroutines ..
 !     DOUBLE PRECISION   dlamch
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          abs, max, min, sqrt
!     ..
!     .. Executable Statements ..
!
      x_is_nan = disnan( x )
      y_is_nan = disnan( y )
      IF ( x_is_nan ) dlapy2 = x
      IF ( y_is_nan ) dlapy2 = y
      hugeval = dlamch( 'Overflow' )
!
      IF ( .NOT.( x_is_nan.OR.y_is_nan ) ) THEN
         xabs = abs( x )
         yabs = abs( y )
         w = max( xabs, yabs )
         z = min( xabs, yabs )
         IF( z.EQ.zero .OR. w.GT.hugeval ) THEN
            dlapy2 = w
         ELSE
            dlapy2 = w*sqrt( one+( z / w )**2 )
         END IF
      END IF
      RETURN
!
!     End of DLAPY2
!

      END


!> \brief \b DLAHR2 reduces the specified number of first columns of a general rectangular matrix A so that elements below the specified subdiagonal are zero, and returns auxiliary matrices which are needed to apply the transformation to the unreduced part of A.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download DLAHR2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlahr2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlahr2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlahr2.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAHR2( N, K, NB, A, LDA, TAU, T, LDT, Y, LDY )
!
!       .. Scalar Arguments ..
!       INTEGER            K, LDA, LDT, LDY, N, NB
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION  A( LDA, * ), T( LDT, NB ), TAU( NB ),
!      $                   Y( LDY, NB )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAHR2 reduces the first NB columns of A real general n-BY-(n-k+1)
!> matrix A so that elements below the k-th subdiagonal are zero. The
!> reduction is performed by an orthogonal similarity transformation
!> Q**T * A * Q. The routine returns the matrices V and T which determine
!> Q as a block reflector I - V*T*V**T, and also the matrix Y = A * V * T.
!>
!> This is an auxiliary routine called by DGEHRD.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The offset for the reduction. Elements below the k-th
!>          subdiagonal in the first NB columns are reduced to zero.
!>          K < N.
!> \endverbatim
!>
!> \param[in] NB
!> \verbatim
!>          NB is INTEGER
!>          The number of columns to be reduced.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N-K+1)
!>          On entry, the n-by-(n-k+1) general matrix A.
!>          On exit, the elements on and above the k-th subdiagonal in
!>          the first NB columns are overwritten with the corresponding
!>          elements of the reduced matrix; the elements below the k-th
!>          subdiagonal, with the array TAU, represent the matrix Q as a
!>          product of elementary reflectors. The other columns of A are
!>          unchanged. See Further Details.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION array, dimension (NB)
!>          The scalar factors of the elementary reflectors. See Further
!>          Details.
!> \endverbatim
!>
!> \param[out] T
!> \verbatim
!>          T is DOUBLE PRECISION array, dimension (LDT,NB)
!>          The upper triangular matrix T.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of the array T.  LDT >= NB.
!> \endverbatim
!>
!> \param[out] Y
!> \verbatim
!>          Y is DOUBLE PRECISION array, dimension (LDY,NB)
!>          The n-by-nb matrix Y.
!> \endverbatim
!>
!> \param[in] LDY
!> \verbatim
!>          LDY is INTEGER
!>          The leading dimension of the array Y. LDY >= N.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup lahr2
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The matrix Q is represented as a product of nb elementary reflectors
!>
!>     Q = H(1) H(2) . . . H(nb).
!>
!>  Each H(i) has the form
!>
!>     H(i) = I - tau * v * v**T
!>
!>  where tau is a real scalar, and v is a real vector with
!>  v(1:i+k-1) = 0, v(i+k) = 1; v(i+k+1:n) is stored on exit in
!>  A(i+k+1:n,i), and tau in TAU(i).
!>
!>  The elements of the vectors v together form the (n-k+1)-by-nb matrix
!>  V which is needed, with T and Y, to apply the transformation to the
!>  unreduced part of the matrix, using an update of the form:
!>  A := (I - V*T*V**T) * (A - Y*V**T).
!>
!>  The contents of A on exit are illustrated by the following example
!>  with n = 7, k = 3 and nb = 2:
!>
!>     ( a   a   a   a   a )
!>     ( a   a   a   a   a )
!>     ( a   a   a   a   a )
!>     ( h   h   a   a   a )
!>     ( v1  h   a   a   a )
!>     ( v1  v2  a   a   a )
!>     ( v1  v2  a   a   a )
!>
!>  where a denotes an element of the original matrix A, h denotes a
!>  modified element of the upper Hessenberg matrix H, and vi denotes an
!>  element of the vector defining H(i).
!>
!>  This subroutine is a slight modification of LAPACK-3.0's DLAHRD
!>  incorporating improvements proposed by Quintana-Orti and Van de
!>  Gejin. Note that the entries of A(1:K,2:NB) differ from those
!>  returned by the original LAPACK-3.0's DLAHRD routine. (This
!>  subroutine is not backward compatible with LAPACK-3.0's DLAHRD.)
!> \endverbatim
!
!> \par References:
!  ================
!>
!>  Gregorio Quintana-Orti and Robert van de Geijn, "Improving the
!>  performance of reduction to Hessenberg form," ACM Transactions on
!>  Mathematical Software, 32(2):180-194, June 2006.
!>
!  =====================================================================

      SUBROUTINE dlahr2( N, K, NB, A, LDA, TAU, T, LDT, Y, LDY )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            K, LDA, LDT, LDY, N, NB
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION  A( LDA, * ), T( LDT, NB ), TAU( NB ), &
     &                   Y( LDY, NB )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      parameter( zero = 0.0d+0, &
     &                     one = 1.0d+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION  EI
!     ..
!     .. External Subroutines ..
!      EXTERNAL           daxpy, dcopy, dgemm, dgemv, dlacpy, &
!     &                   dlarfg, dscal, dtrmm, dtrmv
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          min
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF( n.LE.1 ) &
     &   RETURN
!
      DO 10 i = 1, nb
         IF( i.GT.1 ) THEN
!
!           Update A(K+1:N,I)
!
!           Update I-th column of A - Y * V**T
!
            CALL dgemv( 'NO TRANSPOSE', n-k, i-1, -one, y(k+1,1), &
     &                  ldy, &
     &                  a( k+i-1, 1 ), lda, one, a( k+1, i ), 1 )
!
!           Apply I - V * T**T * V**T to this column (call it b) from the
!           left, using the last column of T as workspace
!
!           Let  V = ( V1 )   and   b = ( b1 )   (first I-1 rows)
!                    ( V2 )             ( b2 )
!
!           where V1 is unit lower triangular
!
!           w := V1**T * b1
!
            CALL dcopy( i-1, a( k+1, i ), 1, t( 1, nb ), 1 )
            CALL dtrmv( 'Lower', 'Transpose', 'UNIT', &
     &                  i-1, a( k+1, 1 ), &
     &                  lda, t( 1, nb ), 1 )
!
!           w := w + V2**T * b2
!
            CALL dgemv( 'Transpose', n-k-i+1, i-1, &
     &                  one, a( k+i, 1 ), &
     &                  lda, a( k+i, i ), 1, one, t( 1, nb ), 1 )
!
!           w := T**T * w
!
            CALL dtrmv( 'Upper', 'Transpose', 'NON-UNIT', &
     &                  i-1, t, ldt, &
     &                  t( 1, nb ), 1 )
!
!           b2 := b2 - V2*w
!
            CALL dgemv( 'NO TRANSPOSE', n-k-i+1, i-1, -one, &
     &                  a( k+i, 1 ), &
     &                  lda, t( 1, nb ), 1, one, a( k+i, i ), 1 )
!
!           b1 := b1 - V1*w
!
            CALL dtrmv( 'Lower', 'NO TRANSPOSE', &
     &                  'UNIT', i-1, &
     &                  a( k+1, 1 ), lda, t( 1, nb ), 1 )
            CALL daxpy( i-1, -one, t( 1, nb ), 1, a( k+1, i ), 1 )
!
            a( k+i-1, i-1 ) = ei
         END IF
!
!        Generate the elementary reflector H(I) to annihilate
!        A(K+I+1:N,I)
!
         CALL dlarfg( n-k-i+1, a( k+i, i ), a( min( k+i+1, n ), i ), &
     &                1, &
     &                tau( i ) )
         ei = a( k+i, i )
         a( k+i, i ) = one
!
!        Compute  Y(K+1:N,I)
!
         CALL dgemv( 'NO TRANSPOSE', n-k, n-k-i+1, &
     &               one, a( k+1, i+1 ), &
     &               lda, a( k+i, i ), 1, zero, y( k+1, i ), 1 )
         CALL dgemv( 'Transpose', n-k-i+1, i-1, &
     &               one, a( k+i, 1 ), lda, &
     &               a( k+i, i ), 1, zero, t( 1, i ), 1 )
         CALL dgemv( 'NO TRANSPOSE', n-k, i-1, -one, &
     &               y( k+1, 1 ), ldy, &
     &               t( 1, i ), 1, one, y( k+1, i ), 1 )
         CALL dscal( n-k, tau( i ), y( k+1, i ), 1 )
!
!        Compute T(1:I,I)
!
         CALL dscal( i-1, -tau( i ), t( 1, i ), 1 )
         CALL dtrmv( 'Upper', 'No Transpose', 'NON-UNIT', &
     &               i-1, t, ldt, &
     &               t( 1, i ), 1 )
         t( i, i ) = tau( i )
!
10    CONTINUE
      a( k+nb, nb ) = ei
!
!     Compute Y(1:K,1:NB)
!
      CALL dlacpy( 'ALL', k, nb, a( 1, 2 ), lda, y, ldy )
      CALL dtrmm( 'RIGHT', 'Lower', 'NO TRANSPOSE', &
     &            'UNIT', k, nb, &
     &            one, a( k+1, 1 ), lda, y, ldy )
      IF( n.GT.k+nb ) &
     &   CALL dgemm( 'NO TRANSPOSE', 'NO TRANSPOSE', k, &
     &               nb, n-k-nb, one, &
     &               a( 1, 2+nb ), lda, a( k+1+nb, 1 ), lda, one, y, &
     &               ldy )
      CALL dtrmm( 'RIGHT', 'Upper', 'NO TRANSPOSE', &
     &            'NON-UNIT', k, nb, &
     &            one, t, ldt, y, ldy )
!
      RETURN
!
!     End of DLAHR2
!

      END

!> \brief \b DLAQR5 performs a single small-bulge multi-shift QR sweep.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download DLAQR5 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaqr5.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaqr5.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaqr5.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAQR5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NSHFTS,
!                          SR, SI, H, LDH, ILOZ, IHIZ, Z, LDZ, V, LDV, U,
!                          LDU, NV, WV, LDWV, NH, WH, LDWH )
!
!       .. Scalar Arguments ..
!       INTEGER            IHIZ, ILOZ, KACC22, KBOT, KTOP, LDH, LDU, LDV,
!      $                   LDWH, LDWV, LDZ, N, NH, NSHFTS, NV
!       LOGICAL            WANTT, WANTZ
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   H( LDH, * ), SI( * ), SR( * ), U( LDU, * ),
!      $                   V( LDV, * ), WH( LDWH, * ), WV( LDWV, * ),
!      $                   Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    DLAQR5, called by DLAQR0, performs a
!>    single small-bulge multi-shift QR sweep.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] WANTT
!> \verbatim
!>          WANTT is LOGICAL
!>             WANTT = .true. if the quasi-triangular Schur factor
!>             is being computed.  WANTT is set to .false. otherwise.
!> \endverbatim
!>
!> \param[in] WANTZ
!> \verbatim
!>          WANTZ is LOGICAL
!>             WANTZ = .true. if the orthogonal Schur factor is being
!>             computed.  WANTZ is set to .false. otherwise.
!> \endverbatim
!>
!> \param[in] KACC22
!> \verbatim
!>          KACC22 is INTEGER with value 0, 1, or 2.
!>             Specifies the computation mode of far-from-diagonal
!>             orthogonal updates.
!>        = 0: DLAQR5 does not accumulate reflections and does not
!>             use matrix-matrix multiply to update far-from-diagonal
!>             matrix entries.
!>        = 1: DLAQR5 accumulates reflections and uses matrix-matrix
!>             multiply to update the far-from-diagonal matrix entries.
!>        = 2: Same as KACC22 = 1. This option used to enable exploiting
!>             the 2-by-2 structure during matrix multiplications, but
!>             this is no longer supported.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>             N is the order of the Hessenberg matrix H upon which this
!>             subroutine operates.
!> \endverbatim
!>
!> \param[in] KTOP
!> \verbatim
!>          KTOP is INTEGER
!> \endverbatim
!>
!> \param[in] KBOT
!> \verbatim
!>          KBOT is INTEGER
!>             These are the first and last rows and columns of an
!>             isolated diagonal block upon which the QR sweep is to be
!>             applied. It is assumed without a check that
!>                       either KTOP = 1  or   H(KTOP,KTOP-1) = 0
!>             and
!>                       either KBOT = N  or   H(KBOT+1,KBOT) = 0.
!> \endverbatim
!>
!> \param[in] NSHFTS
!> \verbatim
!>          NSHFTS is INTEGER
!>             NSHFTS gives the number of simultaneous shifts.  NSHFTS
!>             must be positive and even.
!> \endverbatim
!>
!> \param[in,out] SR
!> \verbatim
!>          SR is DOUBLE PRECISION array, dimension (NSHFTS)
!> \endverbatim
!>
!> \param[in,out] SI
!> \verbatim
!>          SI is DOUBLE PRECISION array, dimension (NSHFTS)
!>             SR contains the real parts and SI contains the imaginary
!>             parts of the NSHFTS shifts of origin that define the
!>             multi-shift QR sweep.  On output SR and SI may be
!>             reordered.
!> \endverbatim
!>
!> \param[in,out] H
!> \verbatim
!>          H is DOUBLE PRECISION array, dimension (LDH,N)
!>             On input H contains a Hessenberg matrix.  On output a
!>             multi-shift QR sweep with shifts SR(J)+i*SI(J) is applied
!>             to the isolated diagonal block in rows and columns KTOP
!>             through KBOT.
!> \endverbatim
!>
!> \param[in] LDH
!> \verbatim
!>          LDH is INTEGER
!>             LDH is the leading dimension of H just as declared in the
!>             calling procedure.  LDH >= MAX(1,N).
!> \endverbatim
!>
!> \param[in] ILOZ
!> \verbatim
!>          ILOZ is INTEGER
!> \endverbatim
!>
!> \param[in] IHIZ
!> \verbatim
!>          IHIZ is INTEGER
!>             Specify the rows of Z to which transformations must be
!>             applied if WANTZ is .TRUE.. 1 <= ILOZ <= IHIZ <= N
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is DOUBLE PRECISION array, dimension (LDZ,IHIZ)
!>             If WANTZ = .TRUE., then the QR Sweep orthogonal
!>             similarity transformation is accumulated into
!>             Z(ILOZ:IHIZ,ILOZ:IHIZ) from the right.
!>             If WANTZ = .FALSE., then Z is unreferenced.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>             LDA is the leading dimension of Z just as declared in
!>             the calling procedure. LDZ >= N.
!> \endverbatim
!>
!> \param[out] V
!> \verbatim
!>          V is DOUBLE PRECISION array, dimension (LDV,NSHFTS/2)
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is INTEGER
!>             LDV is the leading dimension of V as declared in the
!>             calling procedure.  LDV >= 3.
!> \endverbatim
!>
!> \param[out] U
!> \verbatim
!>          U is DOUBLE PRECISION array, dimension (LDU,2*NSHFTS)
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>             LDU is the leading dimension of U just as declared in the
!>             in the calling subroutine.  LDU >= 2*NSHFTS.
!> \endverbatim
!>
!> \param[in] NV
!> \verbatim
!>          NV is INTEGER
!>             NV is the number of rows in WV agailable for workspace.
!>             NV >= 1.
!> \endverbatim
!>
!> \param[out] WV
!> \verbatim
!>          WV is DOUBLE PRECISION array, dimension (LDWV,2*NSHFTS)
!> \endverbatim
!>
!> \param[in] LDWV
!> \verbatim
!>          LDWV is INTEGER
!>             LDWV is the leading dimension of WV as declared in the
!>             in the calling subroutine.  LDWV >= NV.
!> \endverbatim
!
!> \param[in] NH
!> \verbatim
!>          NH is INTEGER
!>             NH is the number of columns in array WH available for
!>             workspace. NH >= 1.
!> \endverbatim
!>
!> \param[out] WH
!> \verbatim
!>          WH is DOUBLE PRECISION array, dimension (LDWH,NH)
!> \endverbatim
!>
!> \param[in] LDWH
!> \verbatim
!>          LDWH is INTEGER
!>             Leading dimension of WH just as declared in the
!>             calling procedure.  LDWH >= 2*NSHFTS.
!> \endverbatim
!>
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup laqr5
!
!> \par Contributors:
!  ==================
!>
!>       Karen Braman and Ralph Byers, Department of Mathematics,
!>       University of Kansas, USA
!>
!>       Lars Karlsson, Daniel Kressner, and Bruno Lang
!>
!>       Thijs Steel, Department of Computer science,
!>       KU Leuven, Belgium
!
!> \par References:
!  ================
!>
!>       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
!>       Algorithm Part I: Maintaining Well Focused Shifts, and Level 3
!>       Performance, SIAM Journal of Matrix Analysis, volume 23, pages
!>       929--947, 2002.
!>
!>       Lars Karlsson, Daniel Kressner, and Bruno Lang, Optimally packed
!>       chains of bulges in multishift QR algorithms.
!>       ACM Trans. Math. Softw. 40, 2, Article 12 (February 2014).
!>
!  =====================================================================

      SUBROUTINE dlaqr5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NSHFTS, &
     &                   SR, SI, H, LDH, ILOZ, IHIZ, Z, LDZ, V, LDV, U, &
     &                   LDU, NV, WV, LDWV, NH, WH, LDWH )
      IMPLICIT NONE
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            IHIZ, ILOZ, KACC22, KBOT, KTOP, LDH, LDU, LDV, &
     &                   LDWH, LDWV, LDZ, N, NH, NSHFTS, NV
      LOGICAL            WANTT, WANTZ
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   H( LDH, * ), SI( * ), SR( * ), U( LDU, * ), &
     &                   V( LDV, * ), WH( LDWH, * ), WV( LDWV, * ), &
     &                   z( ldz, * )
!     ..
!
!  ================================================================
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0d0, one = 1.0d0 )
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION   ALPHA, BETA, H11, H12, H21, H22, REFSUM, &
     &                   SAFMAX, SAFMIN, SCL, SMLNUM, SWAP, T1, T2, &
     &                   t3, tst1, tst2, ulp
      INTEGER            I, I2, I4, INCOL, J, JBOT, JCOL, JLEN, &
     &                   JROW, JTOP, K, K1, KDU, KMS, KRCOL, &
     &                   m, m22, mbot, mtop, nbmps, ndcol, &
     &                   ns, nu
      LOGICAL            ACCUM, BMP22
!     ..
!     .. External Functions ..
!      DOUBLE PRECISION   DLAMCH
!      EXTERNAL           DLAMCH
!     ..
!     .. Intrinsic Functions ..
!
      INTRINSIC          abs, dble, max, min, mod
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION   VT( 3 )
!     ..
!     .. External Subroutines ..
!      EXTERNAL           dgemm, dlacpy, dlaqr1, dlarfg, dlaset, &
!     &                   dtrmm
!     ..
!     .. Executable Statements ..
!
!     ==== If there are no shifts, then there is nothing to do. ====
!
      IF( nshfts.LT.2 ) &
     &   RETURN
!
!     ==== If the active block is empty or 1-by-1, then there
!     .    is nothing to do. ====
!
      IF( ktop.GE.kbot ) &
     &   RETURN
!
!     ==== Shuffle shifts into pairs of real shifts and pairs
!     .    of complex conjugate shifts assuming complex
!     .    conjugate shifts are already adjacent to one
!     .    another. ====
!
      DO 10 i = 1, nshfts - 2, 2
         IF( si( i ).NE.-si( i+1 ) ) THEN
!
            swap = sr( i )
            sr( i ) = sr( i+1 )
            sr( i+1 ) = sr( i+2 )
            sr( i+2 ) = swap
!
            swap = si( i )
            si( i ) = si( i+1 )
            si( i+1 ) = si( i+2 )
            si( i+2 ) = swap
         END IF
10    CONTINUE
!
!     ==== NSHFTS is supposed to be even, but if it is odd,
!     .    then simply reduce it by one.  The shuffle above
!     .    ensures that the dropped shift is real and that
!     .    the remaining shifts are paired. ====
!
      ns = nshfts - mod( nshfts, 2 )
!
!     ==== Machine constants for deflation ====
!
      safmin = dlamch( 'SAFE MINIMUM' )
      safmax = one / safmin
      ulp = dlamch( 'PRECISION' )
      smlnum = safmin*( dble( n ) / ulp )
!
!     ==== Use accumulated reflections to update far-from-diagonal
!     .    entries ? ====
!
      accum = ( kacc22.EQ.1 ) .OR. ( kacc22.EQ.2 )
!
!     ==== clear trash ====
!
      IF( ktop+2.LE.kbot ) &
     &   h( ktop+2, ktop ) = zero
!
!     ==== NBMPS = number of 2-shift bulges in the chain ====
!
      nbmps = ns / 2
!
!     ==== KDU = width of slab ====
!
      kdu = 4*nbmps
!
!     ==== Create and chase chains of NBMPS bulges ====
!
      DO 180 incol = ktop - 2*nbmps + 1, kbot - 2, 2*nbmps
!
!        JTOP = Index from which updates from the right start.
!
         IF( accum ) THEN
            jtop = max( ktop, incol )
         ELSE IF( wantt ) THEN
            jtop = 1
         ELSE
            jtop = ktop
         END IF
!
         ndcol = incol + kdu
         IF( accum ) &
     &      CALL dlaset( 'ALL', kdu, kdu, zero, one, u, ldu )
!
!        ==== Near-the-diagonal bulge chase.  The following loop
!        .    performs the near-the-diagonal part of a small bulge
!        .    multi-shift QR sweep.  Each 4*NBMPS column diagonal
!        .    chunk extends from column INCOL to column NDCOL
!        .    (including both column INCOL and column NDCOL). The
!        .    following loop chases a 2*NBMPS+1 column long chain of
!        .    NBMPS bulges 2*NBMPS columns to the right.  (INCOL
!        .    may be less than KTOP and and NDCOL may be greater than
!        .    KBOT indicating phantom columns from which to chase
!        .    bulges before they are actually introduced or to which
!        .    to chase bulges beyond column KBOT.)  ====
!
         DO 145 krcol = incol, min( incol+2*nbmps-1, kbot-2 )
!
!           ==== Bulges number MTOP to MBOT are active double implicit
!           .    shift bulges.  There may or may not also be small
!           .    2-by-2 bulge, if there is room.  The inactive bulges
!           .    (if any) must wait until the active bulges have moved
!           .    down the diagonal to make room.  The phantom matrix
!           .    paradigm described above helps keep track.  ====
!
            mtop = max( 1, ( ktop-krcol ) / 2+1 )
            mbot = min( nbmps, ( kbot-krcol-1 ) / 2 )
            m22 = mbot + 1
            bmp22 = ( mbot.LT.nbmps ) .AND. ( krcol+2*( m22-1 ) ).EQ. &
     &              ( kbot-2 )
!
!           ==== Generate reflections to chase the chain right
!           .    one column.  (The minimum value of K is KTOP-1.) ====
!
            IF ( bmp22 ) THEN
!
!              ==== Special case: 2-by-2 reflection at bottom treated
!              .    separately ====
!
               k = krcol + 2*( m22-1 )
               IF( k.EQ.ktop-1 ) THEN
                  CALL dlaqr1( 2, h( k+1, k+1 ), ldh, sr( 2*m22-1 ), &
     &                         si( 2*m22-1 ), sr( 2*m22 ), si( 2*m22 ), &
     &                         v( 1, m22 ) )
                  beta = v( 1, m22 )
                  CALL dlarfg( 2, beta, v( 2, m22 ), 1, v( 1, m22 ) )
               ELSE
                  beta = h( k+1, k )
                  v( 2, m22 ) = h( k+2, k )
                  CALL dlarfg( 2, beta, v( 2, m22 ), 1, v( 1, m22 ) )
                  h( k+1, k ) = beta
                  h( k+2, k ) = zero
               END IF

!
!              ==== Perform update from right within 
!              .    computational window. ====
!
               t1 = v( 1, m22 )
               t2 = t1*v( 2, m22 )
               DO 30 j = jtop, min( kbot, k+3 )
                  refsum = h( j, k+1 ) + v( 2, m22 )*h( j, k+2 )
                  h( j, k+1 ) = h( j, k+1 ) - refsum*t1
                  h( j, k+2 ) = h( j, k+2 ) - refsum*t2
30             CONTINUE
!
!              ==== Perform update from left within 
!              .    computational window. ====
!
               IF( accum ) THEN
                  jbot = min( ndcol, kbot )
               ELSE IF( wantt ) THEN
                  jbot = n
               ELSE
                  jbot = kbot
               END IF
               t1 = v( 1, m22 )
               t2 = t1*v( 2, m22 )
               DO 40 j = k+1, jbot
                  refsum = h( k+1, j ) + v( 2, m22 )*h( k+2, j )
                  h( k+1, j ) = h( k+1, j ) - refsum*t1
                  h( k+2, j ) = h( k+2, j ) - refsum*t2
40             CONTINUE
!
!              ==== The following convergence test requires that
!              .    the tradition small-compared-to-nearby-diagonals
!              .    criterion and the Ahues & Tisseur (LAWN 122, 1997)
!              .    criteria both be satisfied.  The latter improves
!              .    accuracy in some examples. Falling back on an
!              .    alternate convergence criterion when TST1 or TST2
!              .    is zero (as done here) is traditional but probably
!              .    unnecessary. ====
!
               IF( k.GE.ktop ) THEN
                  IF( h( k+1, k ).NE.zero ) THEN
                     tst1 = abs( h( k, k ) ) + abs( h( k+1, k+1 ) )
                     IF( tst1.EQ.zero ) THEN
                        IF( k.GE.ktop+1 ) &
     &                     tst1 = tst1 + abs( h( k, k-1 ) )
                        IF( k.GE.ktop+2 ) &
     &                     tst1 = tst1 + abs( h( k, k-2 ) )
                        IF( k.GE.ktop+3 ) &
     &                     tst1 = tst1 + abs( h( k, k-3 ) )
                        IF( k.LE.kbot-2 ) &
     &                     tst1 = tst1 + abs( h( k+2, k+1 ) )
                        IF( k.LE.kbot-3 ) &
     &                     tst1 = tst1 + abs( h( k+3, k+1 ) )
                        IF( k.LE.kbot-4 ) &
     &                     tst1 = tst1 + abs( h( k+4, k+1 ) )
                     END IF
                     IF( abs( h( k+1, k ) ) &
     &                   .LE.max( smlnum, ulp*tst1 ) ) THEN
                        h12 = max( abs( h( k+1, k ) ), &
     &                             abs( h( k, k+1 ) ) )
                        h21 = min( abs( h( k+1, k ) ), &
     &                             abs( h( k, k+1 ) ) )
                        h11 = max( abs( h( k+1, k+1 ) ), &
     &                             abs( h( k, k )-h( k+1, k+1 ) ) )
                        h22 = min( abs( h( k+1, k+1 ) ), &
     &                        abs( h( k, k )-h( k+1, k+1 ) ) )
                        scl = h11 + h12
                        tst2 = h22*( h11 / scl )
!
                        IF( tst2.EQ.zero .OR. h21*( h12 / scl ).LE. &
     &                      max( smlnum, ulp*tst2 ) ) THEN
                           h( k+1, k ) = zero
                        END IF
                     END IF
                  END IF
               END IF
!
!              ==== Accumulate orthogonal transformations. ====
!
               IF( accum ) THEN
                  kms = k - incol
                  t1 = v( 1, m22 )
                  t2 = t1*v( 2, m22 )
                  DO 50 j = max( 1, ktop-incol ), kdu
                     refsum = u( j, kms+1 ) + v( 2, m22 )*u( j, kms+2 )
                     u( j, kms+1 ) = u( j, kms+1 ) - refsum*t1
                     u( j, kms+2 ) = u( j, kms+2 ) - refsum*t2
50                   CONTINUE
               ELSE IF( wantz ) THEN
                  t1 = v( 1, m22 )
                  t2 = t1*v( 2, m22 )
                  DO 60 j = iloz, ihiz
                     refsum = z( j, k+1 )+v( 2, m22 )*z( j, k+2 )
                     z( j, k+1 ) = z( j, k+1 ) - refsum*t1
                     z( j, k+2 ) = z( j, k+2 ) - refsum*t2
60                CONTINUE
               END IF
            END IF
!
!           ==== Normal case: Chain of 3-by-3 reflections ====
!
            DO 80 m = mbot, mtop, -1
               k = krcol + 2*( m-1 )
               IF( k.EQ.ktop-1 ) THEN
                  CALL dlaqr1( 3, h( ktop, ktop ), ldh, sr( 2*m-1 ), &
     &                         si( 2*m-1 ), sr( 2*m ), si( 2*m ), &
     &                         v( 1, m ) )
                  alpha = v( 1, m )
                  CALL dlarfg( 3, alpha, v( 2, m ), 1, v( 1, m ) )
               ELSE
!
!                 ==== Perform delayed transformation of row below
!                 .    Mth bulge. Exploit fact that first two elements
!                 .    of row are actually zero. ====
!
                  t1 = v( 1, m )
                  t2 = t1*v( 2, m )
                  t3 = t1*v( 3, m )
                  refsum = v( 3, m )*h( k+3, k+2 )
                  h( k+3, k   ) = -refsum*t1
                  h( k+3, k+1 ) = -refsum*t2
                  h( k+3, k+2 ) = h( k+3, k+2 ) - refsum*t3
!
!                 ==== Calculate reflection to move
!                 .    Mth bulge one step. ====
!
                  beta      = h( k+1, k )
                  v( 2, m ) = h( k+2, k )
                  v( 3, m ) = h( k+3, k )
                  CALL dlarfg( 3, beta, v( 2, m ), 1, v( 1, m ) )
!
!                 ==== A Bulge may collapse because of vigilant
!                 .    deflation or destructive underflow.  In the
!                 .    underflow case, try the two-small-subdiagonals
!                 .    trick to try to reinflate the bulge.  ====
!
                  IF( h( k+3, k ).NE.zero .OR. h( k+3, k+1 ).NE. &
     &                zero .OR. h( k+3, k+2 ).EQ.zero ) THEN
!
!                    ==== Typical case: not collapsed (yet). ====
!
                     h( k+1, k ) = beta
                     h( k+2, k ) = zero
                     h( k+3, k ) = zero
                  ELSE
!
!                    ==== Atypical case: collapsed.  Attempt to
!                    .    reintroduce ignoring H(K+1,K) and H(K+2,K).
!                    .    If the fill resulting from the new
!                    .    reflector is too large, then abandon it.
!                    .    Otherwise, use the new one. ====
!
                     CALL dlaqr1( 3, h( k+1, k+1 ), ldh, sr( 2*m-1 ), &
     &                            si( 2*m-1 ), sr( 2*m ), si( 2*m ), &
     &                            vt )
                     alpha = vt( 1 )
                     CALL dlarfg( 3, alpha, vt( 2 ), 1, vt( 1 ) )
                     t1 = vt( 1 )
                     t2 = t1*vt( 2 )
                     t3 = t1*vt( 3 )
                     refsum = h( k+1, k ) + vt( 2 )*h( k+2, k )
!
                     IF( abs( h( k+2, k )-refsum*t2 )+ &
     &                   abs( refsum*t3 ).GT.ulp* &
     &                   ( abs( h( k, k ) )+abs( h( k+1, &
     &                   k+1 ) )+abs( h( k+2, k+2 ) ) ) ) THEN
!
!                       ==== Starting a new bulge here would
!                       .    create non-negligible fill.  Use
!                       .    the old one with trepidation. ====
!
                        h( k+1, k ) = beta
                        h( k+2, k ) = zero
                        h( k+3, k ) = zero
                     ELSE
!
!                       ==== Starting a new bulge here would
!                       .    create only negligible fill.
!                       .    Replace the old reflector with
!                       .    the new one. ====
!
                        h( k+1, k ) = h( k+1, k ) - refsum*t1
                        h( k+2, k ) = zero
                        h( k+3, k ) = zero
                        v( 1, m ) = vt( 1 )
                        v( 2, m ) = vt( 2 )
                        v( 3, m ) = vt( 3 )
                     END IF
                  END IF
               END IF
!
!              ====  Apply reflection from the right and
!              .     the first column of update from the left.
!              .     These updates are required for the vigilant
!              .     deflation check. We still delay most of the
!              .     updates from the left for efficiency. ====      
!
               t1 = v( 1, m )
               t2 = t1*v( 2, m )
               t3 = t1*v( 3, m )
               DO 70 j = jtop, min( kbot, k+3 )
                  refsum = h( j, k+1 ) + v( 2, m )*h( j, k+2 ) &
     &                     + v( 3, m )*h( j, k+3 )
                  h( j, k+1 ) = h( j, k+1 ) - refsum*t1
                  h( j, k+2 ) = h( j, k+2 ) - refsum*t2
                  h( j, k+3 ) = h( j, k+3 ) - refsum*t3
70             CONTINUE
!
!              ==== Perform update from left for subsequent
!              .    column. ====
!
               refsum = h( k+1, k+1 ) + v( 2, m )*h( k+2, k+1 ) &
     &                  + v( 3, m )*h( k+3, k+1 )
               h( k+1, k+1 ) = h( k+1, k+1 ) - refsum*t1
               h( k+2, k+1 ) = h( k+2, k+1 ) - refsum*t2
               h( k+3, k+1 ) = h( k+3, k+1 ) - refsum*t3
!
!              ==== The following convergence test requires that
!              .    the tradition small-compared-to-nearby-diagonals
!              .    criterion and the Ahues & Tisseur (LAWN 122, 1997)
!              .    criteria both be satisfied.  The latter improves
!              .    accuracy in some examples. Falling back on an
!              .    alternate convergence criterion when TST1 or TST2
!              .    is zero (as done here) is traditional but probably
!              .    unnecessary. ====
!
               IF( k.LT.ktop) &
     &              cycle
               IF( h( k+1, k ).NE.zero ) THEN
                  tst1 = abs( h( k, k ) ) + abs( h( k+1, k+1 ) )
                  IF( tst1.EQ.zero ) THEN
                     IF( k.GE.ktop+1 ) &
     &                  tst1 = tst1 + abs( h( k, k-1 ) )
                     IF( k.GE.ktop+2 ) &
     &                  tst1 = tst1 + abs( h( k, k-2 ) )
                     IF( k.GE.ktop+3 ) &
     &                  tst1 = tst1 + abs( h( k, k-3 ) )
                     IF( k.LE.kbot-2 ) &
     &                  tst1 = tst1 + abs( h( k+2, k+1 ) )
                     IF( k.LE.kbot-3 ) &
     &                  tst1 = tst1 + abs( h( k+3, k+1 ) )
                     IF( k.LE.kbot-4 ) &
     &                  tst1 = tst1 + abs( h( k+4, k+1 ) )
                  END IF
                  IF( abs( h( k+1, k ) ).LE.max( smlnum, ulp*tst1 ) ) &
     &                 THEN
                     h12 = max( abs( h( k+1, k ) ), abs( h( k, k+1 ) ) )
                     h21 = min( abs( h( k+1, k ) ), abs( h( k, k+1 ) ) )
                     h11 = max( abs( h( k+1, k+1 ) ), &
     &                     abs( h( k, k )-h( k+1, k+1 ) ) )
                     h22 = min( abs( h( k+1, k+1 ) ), &
     &                     abs( h( k, k )-h( k+1, k+1 ) ) )
                     scl = h11 + h12
                     tst2 = h22*( h11 / scl )
!
                     IF( tst2.EQ.zero .OR. h21*( h12 / scl ).LE. &
     &                   max( smlnum, ulp*tst2 ) ) THEN
                        h( k+1, k ) = zero
                     END IF
                  END IF
               END IF
80          CONTINUE
!
!           ==== Multiply H by reflections from the left ====
!
            IF( accum ) THEN
               jbot = min( ndcol, kbot )
            ELSE IF( wantt ) THEN
               jbot = n
            ELSE
               jbot = kbot
            END IF
!
            DO 100 m = mbot, mtop, -1
               k = krcol + 2*( m-1 )
               t1 = v( 1, m )
               t2 = t1*v( 2, m )
               t3 = t1*v( 3, m )
               DO 90 j = max( ktop, krcol + 2*m ), jbot
                  refsum = h( k+1, j ) + v( 2, m )*h( k+2, j ) &
     &                     + v( 3, m )*h( k+3, j )
                  h( k+1, j ) = h( k+1, j ) - refsum*t1
                  h( k+2, j ) = h( k+2, j ) - refsum*t2
                  h( k+3, j ) = h( k+3, j ) - refsum*t3
90             CONTINUE
100         CONTINUE
!
!           ==== Accumulate orthogonal transformations. ====
!
            IF( accum ) THEN
!
!              ==== Accumulate U. (If needed, update Z later
!              .    with an efficient matrix-matrix
!              .    multiply.) ====
!
               DO 120 m = mbot, mtop, -1
                  k = krcol + 2*( m-1 )
                  kms = k - incol
                  i2 = max( 1, ktop-incol )
                  i2 = max( i2, kms-(krcol-incol)+1 )
                  i4 = min( kdu, krcol + 2*( mbot-1 ) - incol + 5 )
                  t1 = v( 1, m )
                  t2 = t1*v( 2, m )
                  t3 = t1*v( 3, m )
                  DO 110 j = i2, i4
                     refsum = u( j, kms+1 ) + v( 2, m )*u( j, kms+2 ) &
     &                        + v( 3, m )*u( j, kms+3 )
                     u( j, kms+1 ) = u( j, kms+1 ) - refsum*t1
                     u( j, kms+2 ) = u( j, kms+2 ) - refsum*t2
                     u( j, kms+3 ) = u( j, kms+3 ) - refsum*t3
110               CONTINUE
120            CONTINUE
            ELSE IF( wantz ) THEN
!
!              ==== U is not accumulated, so update Z
!              .    now by multiplying by reflections
!              .    from the right. ====
!
               DO 140 m = mbot, mtop, -1
                  k = krcol + 2*( m-1 )
                  t1 = v( 1, m )
                  t2 = t1*v( 2, m )
                  t3 = t1*v( 3, m )
                  DO 130 j = iloz, ihiz
                     refsum = z( j, k+1 ) + v( 2, m )*z( j, k+2 ) &
     &                        + v( 3, m )*z( j, k+3 )
                     z( j, k+1 ) = z( j, k+1 ) - refsum*t1
                     z( j, k+2 ) = z( j, k+2 ) - refsum*t2
                     z( j, k+3 ) = z( j, k+3 ) - refsum*t3
130               CONTINUE
140            CONTINUE
            END IF
!
!           ==== End of near-the-diagonal bulge chase. ====
!
145      CONTINUE
!
!        ==== Use U (if accumulated) to update far-from-diagonal
!        .    entries in H.  If required, use U to update Z as
!        .    well. ====
!
         IF( accum ) THEN
            IF( wantt ) THEN
               jtop = 1
               jbot = n
            ELSE
               jtop = ktop
               jbot = kbot
            END IF
            k1 = max( 1, ktop-incol )
            nu = ( kdu-max( 0, ndcol-kbot ) ) - k1 + 1
!
!           ==== Horizontal Multiply ====
!
            DO 150 jcol = min( ndcol, kbot ) + 1, jbot, nh
               jlen = min( nh, jbot-jcol+1 )
               CALL dgemm( 'C', 'N', nu, jlen, nu, one, u( k1, k1 ), &
     &                        ldu, h( incol+k1, jcol ), ldh, zero, wh, &
     &                        ldwh )
               CALL dlacpy( 'ALL', nu, jlen, wh, ldwh, &
     &                         h( incol+k1, jcol ), ldh )
150         CONTINUE
!
!           ==== Vertical multiply ====
!
            DO 160 jrow = jtop, max( ktop, incol ) - 1, nv
               jlen = min( nv, max( ktop, incol )-jrow )
               CALL dgemm( 'N', 'N', jlen, nu, nu, one, &
     &                     h( jrow, incol+k1 ), ldh, u( k1, k1 ), &
     &                     ldu, zero, wv, ldwv )
               CALL dlacpy( 'ALL', jlen, nu, wv, ldwv, &
     &                      h( jrow, incol+k1 ), ldh )
160         CONTINUE
!
!           ==== Z multiply (also vertical) ====
!
            IF( wantz ) THEN
               DO 170 jrow = iloz, ihiz, nv
                  jlen = min( nv, ihiz-jrow+1 )
                  CALL dgemm( 'N', 'N', jlen, nu, nu, one, &
     &                        z( jrow, incol+k1 ), ldz, u( k1, k1 ), &
     &                        ldu, zero, wv, ldwv )
                  CALL dlacpy( 'ALL', jlen, nu, wv, ldwv, &
     &                         z( jrow, incol+k1 ), ldz )
170            CONTINUE
            END IF
         END IF
180   CONTINUE
!
!     ==== End of DLAQR5 ====
!

      END

!> \brief \b DGEMV
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION ALPHA,BETA
!       INTEGER INCX,INCY,LDA,M,N
!       CHARACTER TRANS
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION A(LDA,*),X(*),Y(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGEMV  performs one of the matrix-vector operations
!>
!>    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
!>
!> where alpha and beta are scalars, x and y are vectors and A is an
!> m by n matrix.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>           On entry, TRANS specifies the operation to be performed as
!>           follows:
!>
!>              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
!>
!>              TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y.
!>
!>              TRANS = 'C' or 'c'   y := alpha*A**T*x + beta*y.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>           On entry, M specifies the number of rows of the matrix A.
!>           M must be at least zero.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry, N specifies the number of columns of the matrix A.
!>           N must be at least zero.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is DOUBLE PRECISION.
!>           On entry, ALPHA specifies the scalar alpha.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension ( LDA, N )
!>           Before entry, the leading m by n part of the array A must
!>           contain the matrix of coefficients.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in the calling (sub) program. LDA must be at least
!>           max( 1, m ).
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension at least
!>           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
!>           and at least
!>           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
!>           Before entry, the incremented array X must contain the
!>           vector x.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>           On entry, INCX specifies the increment for the elements of
!>           X. INCX must not be zero.
!> \endverbatim
!>
!> \param[in] BETA
!> \verbatim
!>          BETA is DOUBLE PRECISION.
!>           On entry, BETA specifies the scalar beta. When BETA is
!>           supplied as zero then Y need not be set on input.
!> \endverbatim
!>
!> \param[in,out] Y
!> \verbatim
!>          Y is DOUBLE PRECISION array, dimension at least
!>           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
!>           and at least
!>           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
!>           Before entry with BETA non-zero, the incremented array Y
!>           must contain the vector y. On exit, Y is overwritten by the
!>           updated vector y.
!>           If either m or n is zero, then Y not referenced and the function
!>           performs a quick return.
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>           On entry, INCY specifies the increment for the elements of
!>           Y. INCY must not be zero.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup gemv
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Level 2 Blas routine.
!>  The vector and matrix arguments are not referenced when N = 0, or M = 0
!>
!>  -- Written on 22-October-1986.
!>     Jack Dongarra, Argonne National Lab.
!>     Jeremy Du Croz, Nag Central Office.
!>     Sven Hammarling, Nag Central Office.
!>     Richard Hanson, Sandia National Labs.
!> \endverbatim
!>
!  =====================================================================

      SUBROUTINE dgemv(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
!
!  -- Reference BLAS level2 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BETA
      INTEGER INCX,INCY,LDA,M,N
      CHARACTER TRANS
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),X(*),Y(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      parameter(one=1.0d+0,zero=0.0d+0)
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY,LENX,LENY
!     ..
!     .. External Functions ..
!      LOGICAL LSAME
!      EXTERNAL lsame
!     ..
!     .. External Subroutines ..
!      EXTERNAL xerbla
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC max
!     ..
!
!     Test the input parameters.
!
      info = 0
      IF (.NOT.lsame(trans,'N') .AND. .NOT.lsame(trans,'T') .AND. &
     &    .NOT.lsame(trans,'C')) THEN
          info = 1
      ELSE IF (m.LT.0) THEN
          info = 2
      ELSE IF (n.LT.0) THEN
          info = 3
      ELSE IF (lda.LT.max(1,m)) THEN
          info = 6
      ELSE IF (incx.EQ.0) THEN
          info = 8
      ELSE IF (incy.EQ.0) THEN
          info = 11
      END IF
      IF (info.NE.0) THEN
          CALL xerbla('DGEMV ',info)
          RETURN
      END IF
!
!     Quick return if possible.
!
      IF ((m.EQ.0) .OR. (n.EQ.0) .OR. &
     &    ((alpha.EQ.zero).AND. (beta.EQ.one))) RETURN
!
!     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
!     up the start points in  X  and  Y.
!
      IF (lsame(trans,'N')) THEN
          lenx = n
          leny = m
      ELSE
          lenx = m
          leny = n
      END IF
      IF (incx.GT.0) THEN
          kx = 1
      ELSE
          kx = 1 - (lenx-1)*incx
      END IF
      IF (incy.GT.0) THEN
          ky = 1
      ELSE
          ky = 1 - (leny-1)*incy
      END IF
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
!     First form  y := beta*y.
!
      IF (beta.NE.one) THEN
          IF (incy.EQ.1) THEN
              IF (beta.EQ.zero) THEN
                  DO 10 i = 1,leny
                      y(i) = zero
10                CONTINUE
              ELSE
                  DO 20 i = 1,leny
                      y(i) = beta*y(i)
20                CONTINUE
              END IF
          ELSE
              iy = ky
              IF (beta.EQ.zero) THEN
                  DO 30 i = 1,leny
                      y(iy) = zero
                      iy = iy + incy
30                CONTINUE
              ELSE
                  DO 40 i = 1,leny
                      y(iy) = beta*y(iy)
                      iy = iy + incy
40                CONTINUE
              END IF
          END IF
      END IF
      IF (alpha.EQ.zero) RETURN
      IF (lsame(trans,'N')) THEN
!
!        Form  y := alpha*A*x + y.
!
          jx = kx
          IF (incy.EQ.1) THEN
              DO 60 j = 1,n
                  temp = alpha*x(jx)
                  DO 50 i = 1,m
                      y(i) = y(i) + temp*a(i,j)
50                CONTINUE
                  jx = jx + incx
60            CONTINUE
          ELSE
              DO 80 j = 1,n
                  temp = alpha*x(jx)
                  iy = ky
                  DO 70 i = 1,m
                      y(iy) = y(iy) + temp*a(i,j)
                      iy = iy + incy
70                CONTINUE
                  jx = jx + incx
80            CONTINUE
          END IF
      ELSE
!
!        Form  y := alpha*A**T*x + y.
!
          jy = ky
          IF (incx.EQ.1) THEN
              DO 100 j = 1,n
                  temp = zero
                  DO 90 i = 1,m
                      temp = temp + a(i,j)*x(i)
90                CONTINUE
                  y(jy) = y(jy) + alpha*temp
                  jy = jy + incy
100           CONTINUE
          ELSE
              DO 120 j = 1,n
                  temp = zero
                  ix = kx
                  DO 110 i = 1,m
                      temp = temp + a(i,j)*x(ix)
                      ix = ix + incx
110               CONTINUE
                  y(jy) = y(jy) + alpha*temp
                  jy = jy + incy
120           CONTINUE
          END IF
      END IF
!
      RETURN
!
!     End of DGEMV
!

      END
!> \brief \b DISNAN tests input for NaN.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download DISNAN + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/disnan.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/disnan.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/disnan.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       LOGICAL FUNCTION DISNAN( DIN )
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION, INTENT(IN) :: DIN
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DISNAN returns .TRUE. if its argument is NaN, and .FALSE.
!> otherwise.  To be replaced by the Fortran 2003 intrinsic in the
!> future.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] DIN
!> \verbatim
!>          DIN is DOUBLE PRECISION
!>          Input to test for NaN.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup isnan
!
!  =====================================================================

      LOGICAL FUNCTION disnan( DIN )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION, INTENT(IN) :: din
!     ..
!
!  =====================================================================
!
!  .. External Functions ..
!      LOGICAL dlaisnan
!      EXTERNAL dlaisnan
!  ..
!  .. Executable Statements ..
      disnan = dlaisnan(din,din)
      RETURN

      END
!> \brief \b DCOPY
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,INCY,N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION DX(*),DY(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    DCOPY copies a vector, x, to a vector, y.
!>    uses unrolled loops for increments equal to 1.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         number of elements in input vector(s)
!> \endverbatim
!>
!> \param[in] DX
!> \verbatim
!>          DX is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>         storage spacing between elements of DX
!> \endverbatim
!>
!> \param[out] DY
!> \verbatim
!>          DY is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>         storage spacing between elements of DY
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup copy
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, linpack, 3/11/78.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================

      SUBROUTINE dcopy(N,DX,INCX,DY,INCY)
!
!  -- Reference BLAS level1 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER I,IX,IY,M,MP1
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC mod
!     ..
      IF (n.LE.0) RETURN
      IF (incx.EQ.1 .AND. incy.EQ.1) THEN
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
         m = mod(n,7)
         IF (m.NE.0) THEN
            DO i = 1,m
               dy(i) = dx(i)
            END DO
            IF (n.LT.7) RETURN
         END IF
         mp1 = m + 1
         DO i = mp1,n,7
            dy(i) = dx(i)
            dy(i+1) = dx(i+1)
            dy(i+2) = dx(i+2)
            dy(i+3) = dx(i+3)
            dy(i+4) = dx(i+4)
            dy(i+5) = dx(i+5)
            dy(i+6) = dx(i+6)
         END DO
      ELSE
!
!        code for unequal increments or equal increments
!          not equal to 1
!
         ix = 1
         iy = 1
         IF (incx.LT.0) ix = (-n+1)*incx + 1
         IF (incy.LT.0) iy = (-n+1)*incy + 1
         DO i = 1,n
            dy(iy) = dx(ix)
            ix = ix + incx
            iy = iy + incy
         END DO
      END IF
      RETURN
!
!     End of DCOPY
!

      END
!> \brief \b DLARFG generates an elementary reflector (Householder matrix).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download DLARFG + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarfg.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarfg.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarfg.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLARFG( N, ALPHA, X, INCX, TAU )
!
!       .. Scalar Arguments ..
!       INTEGER            INCX, N
!       DOUBLE PRECISION   ALPHA, TAU
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   X( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLARFG generates a real elementary reflector H of order n, such
!> that
!>
!>       H * ( alpha ) = ( beta ),   H**T * H = I.
!>           (   x   )   (   0  )
!>
!> where alpha and beta are scalars, and x is an (n-1)-element real
!> vector. H is represented in the form
!>
!>       H = I - tau * ( 1 ) * ( 1 v**T ) ,
!>                     ( v )
!>
!> where tau is a real scalar and v is a real (n-1)-element
!> vector.
!>
!> If the elements of x are all zero, then tau = 0 and H is taken to be
!> the unit matrix.
!>
!> Otherwise  1 <= tau <= 2.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the elementary reflector.
!> \endverbatim
!>
!> \param[in,out] ALPHA
!> \verbatim
!>          ALPHA is DOUBLE PRECISION
!>          On entry, the value alpha.
!>          On exit, it is overwritten with the value beta.
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension
!>                         (1+(N-2)*abs(INCX))
!>          On entry, the vector x.
!>          On exit, it is overwritten with the vector v.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>          The increment between elements of X. INCX > 0.
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION
!>          The value tau.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup larfg
!
!  =====================================================================

      SUBROUTINE dlarfg( N, ALPHA, X, INCX, TAU )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   ALPHA, TAU
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      parameter( one = 1.0d+0, zero = 0.0d+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            J, KNT
      DOUBLE PRECISION   BETA, RSAFMN, SAFMIN, XNORM
!     ..
!     .. External Functions ..
!      DOUBLE PRECISION   DLAMCH, DLAPY2, DNRM2
!      EXTERNAL           dlamch, dlapy2, dnrm2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          abs, sign
!     ..
!     .. External Subroutines ..
!      EXTERNAL           dscal
!     ..
!     .. Executable Statements ..
!
      IF( n.LE.1 ) THEN
         tau = zero
         RETURN
      END IF
!
      xnorm = dnrm2( n-1, x, incx )
!
      IF( xnorm.EQ.zero ) THEN
!
!        H  =  I
!
         tau = zero
      ELSE
!
!        general case
!
         beta = -sign( dlapy2( alpha, xnorm ), alpha )
         safmin = dlamch( 'S' ) / dlamch( 'E' )
         knt = 0
         IF( abs( beta ).LT.safmin ) THEN
!
!           XNORM, BETA may be inaccurate; scale X and recompute them
!
            rsafmn = one / safmin
10          CONTINUE
            knt = knt + 1
            CALL dscal( n-1, rsafmn, x, incx )
            beta = beta*rsafmn
            alpha = alpha*rsafmn
            IF( (abs( beta ).LT.safmin) .AND. (knt .LT. 20) ) &
     &         GO TO 10
!
!           New BETA is at most 1, at least SAFMIN
!
            xnorm = dnrm2( n-1, x, incx )
            beta = -sign( dlapy2( alpha, xnorm ), alpha )
         END IF
         tau = ( beta-alpha ) / beta
         CALL dscal( n-1, one / ( alpha-beta ), x, incx )
!
!        If ALPHA is subnormal, it may lose relative accuracy
!
         DO 20 j = 1, knt
            beta = beta*safmin
20       CONTINUE
         alpha = beta
      END IF
!
      RETURN
!
!     End of DLARFG
!

      END
!> \brief \b DLARFB applies a block reflector or its transpose to a general rectangular matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download DLARFB + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarfb.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarfb.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarfb.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV,
!                          T, LDT, C, LDC, WORK, LDWORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIRECT, SIDE, STOREV, TRANS
!       INTEGER            K, LDC, LDT, LDV, LDWORK, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   C( LDC, * ), T( LDT, * ), V( LDV, * ),
!      $                   WORK( LDWORK, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLARFB applies a real block reflector H or its transpose H**T to a
!> real m by n matrix C, from either the left or the right.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'L': apply H or H**T from the Left
!>          = 'R': apply H or H**T from the Right
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          = 'N': apply H (No transpose)
!>          = 'T': apply H**T (Transpose)
!> \endverbatim
!>
!> \param[in] DIRECT
!> \verbatim
!>          DIRECT is CHARACTER*1
!>          Indicates how H is formed from a product of elementary
!>          reflectors
!>          = 'F': H = H(1) H(2) . . . H(k) (Forward)
!>          = 'B': H = H(k) . . . H(2) H(1) (Backward)
!> \endverbatim
!>
!> \param[in] STOREV
!> \verbatim
!>          STOREV is CHARACTER*1
!>          Indicates how the vectors which define the elementary
!>          reflectors are stored:
!>          = 'C': Columnwise
!>          = 'R': Rowwise
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix C.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix C.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The order of the matrix T (= the number of elementary
!>          reflectors whose product defines the block reflector).
!>          If SIDE = 'L', M >= K >= 0;
!>          if SIDE = 'R', N >= K >= 0.
!> \endverbatim
!>
!> \param[in] V
!> \verbatim
!>          V is DOUBLE PRECISION array, dimension
!>                                (LDV,K) if STOREV = 'C'
!>                                (LDV,M) if STOREV = 'R' and SIDE = 'L'
!>                                (LDV,N) if STOREV = 'R' and SIDE = 'R'
!>          The matrix V. See Further Details.
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is INTEGER
!>          The leading dimension of the array V.
!>          If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M);
!>          if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N);
!>          if STOREV = 'R', LDV >= K.
!> \endverbatim
!>
!> \param[in] T
!> \verbatim
!>          T is DOUBLE PRECISION array, dimension (LDT,K)
!>          The triangular k by k matrix T in the representation of the
!>          block reflector.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of the array T. LDT >= K.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension (LDC,N)
!>          On entry, the m by n matrix C.
!>          On exit, C is overwritten by H*C or H**T*C or C*H or C*H**T.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C. LDC >= max(1,M).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (LDWORK,K)
!> \endverbatim
!>
!> \param[in] LDWORK
!> \verbatim
!>          LDWORK is INTEGER
!>          The leading dimension of the array WORK.
!>          If SIDE = 'L', LDWORK >= max(1,N);
!>          if SIDE = 'R', LDWORK >= max(1,M).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup larfb
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The shape of the matrix V and the storage of the vectors which define
!>  the H(i) is best illustrated by the following example with n = 5 and
!>  k = 3. The triangular part of V (including its diagonal) is not
!>  referenced.
!>
!>  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':
!>
!>               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
!>                   ( v1  1    )                     (     1 v2 v2 v2 )
!>                   ( v1 v2  1 )                     (        1 v3 v3 )
!>                   ( v1 v2 v3 )
!>                   ( v1 v2 v3 )
!>
!>  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':
!>
!>               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
!>                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
!>                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
!>                   (     1 v3 )
!>                   (        1 )
!> \endverbatim
!>
!  =====================================================================

      SUBROUTINE dlarfb( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, &
     &                   LDV, &
     &                   T, LDT, C, LDC, WORK, LDWORK )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          DIRECT, SIDE, STOREV, TRANS
      INTEGER            K, LDC, LDT, LDV, LDWORK, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   C( LDC, * ), T( LDT, * ), V( LDV, * ), &
     &                   WORK( LDWORK, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0d+0 )
!     ..
!     .. Local Scalars ..
      CHARACTER          TRANST
      INTEGER            I, J
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
!      EXTERNAL           dcopy, dgemm, dtrmm
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF( m.LE.0 .OR. n.LE.0 ) &
     &   RETURN
!
      IF( lsame( trans, 'N' ) ) THEN
         transt = 'T'
      ELSE
         transt = 'N'
      END IF
!
      IF( lsame( storev, 'C' ) ) THEN
!
         IF( lsame( direct, 'F' ) ) THEN
!
!           Let  V =  ( V1 )    (first K rows)
!                     ( V2 )
!           where  V1  is unit lower triangular.
!
            IF( lsame( side, 'L' ) ) THEN
!
!              Form  H * C  or  H**T * C  where  C = ( C1 )
!                                                    ( C2 )
!
!              W := C**T * V  =  (C1**T * V1 + C2**T * V2)  (stored in WORK)
!
!              W := C1**T
!
               DO 10 j = 1, k
                  CALL dcopy( n, c( j, 1 ), ldc, work( 1, j ), 1 )
10             CONTINUE
!
!              W := W * V1
!
               CALL dtrmm( 'Right', 'Lower', 'No transpose', 'Unit', &
     &                     n, &
     &                     k, one, v, ldv, work, ldwork )
               IF( m.GT.k ) THEN
!
!                 W := W + C2**T * V2
!
                  CALL dgemm( 'Transpose', 'No transpose', n, k, m-k, &
     &                        one, c( k+1, 1 ), ldc, v( k+1, 1 ), ldv, &
     &                        one, work, ldwork )
               END IF
!
!              W := W * T**T  or  W * T
!
               CALL dtrmm( 'Right', 'Upper', transt, 'Non-unit', n, &
     &                     k, &
     &                     one, t, ldt, work, ldwork )
!
!              C := C - V * W**T
!
               IF( m.GT.k ) THEN
!
!                 C2 := C2 - V2 * W**T
!
                  CALL dgemm( 'No transpose', 'Transpose', m-k, n, k, &
     &                        -one, v( k+1, 1 ), ldv, work, ldwork, one, &
     &                        c( k+1, 1 ), ldc )
               END IF
!
!              W := W * V1**T
!
               CALL dtrmm( 'Right', 'Lower', 'Transpose', 'Unit', n, &
     &                     k, &
     &                     one, v, ldv, work, ldwork )
!
!              C1 := C1 - W**T
!
               DO 30 j = 1, k
                  DO 20 i = 1, n
                     c( j, i ) = c( j, i ) - work( i, j )
20                CONTINUE
30             CONTINUE
!
            ELSE IF( lsame( side, 'R' ) ) THEN
!
!              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
!
!              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
!
!              W := C1
!
               DO 40 j = 1, k
                  CALL dcopy( m, c( 1, j ), 1, work( 1, j ), 1 )
40             CONTINUE
!
!              W := W * V1
!
               CALL dtrmm( 'Right', 'Lower', 'No transpose', 'Unit', &
     &                     m, &
     &                     k, one, v, ldv, work, ldwork )
               IF( n.GT.k ) THEN
!
!                 W := W + C2 * V2
!
                  CALL dgemm( 'No transpose', 'No transpose', m, k, &
     &                        n-k, &
     &                        one, c( 1, k+1 ), ldc, v( k+1, 1 ), ldv, &
     &                        one, work, ldwork )
               END IF
!
!              W := W * T  or  W * T**T
!
               CALL dtrmm( 'Right', 'Upper', trans, 'Non-unit', m, k, &
     &                     one, t, ldt, work, ldwork )
!
!              C := C - W * V**T
!
               IF( n.GT.k ) THEN
!
!                 C2 := C2 - W * V2**T
!
                  CALL dgemm( 'No transpose', 'Transpose', m, n-k, k, &
     &                        -one, work, ldwork, v( k+1, 1 ), ldv, one, &
     &                        c( 1, k+1 ), ldc )
               END IF
!
!              W := W * V1**T
!
               CALL dtrmm( 'Right', 'Lower', 'Transpose', 'Unit', m, &
     &                     k, &
     &                     one, v, ldv, work, ldwork )
!
!              C1 := C1 - W
!
               DO 60 j = 1, k
                  DO 50 i = 1, m
                     c( i, j ) = c( i, j ) - work( i, j )
50                CONTINUE
60             CONTINUE
            END IF
!
         ELSE
!
!           Let  V =  ( V1 )
!                     ( V2 )    (last K rows)
!           where  V2  is unit upper triangular.
!
            IF( lsame( side, 'L' ) ) THEN
!
!              Form  H * C  or  H**T * C  where  C = ( C1 )
!                                                    ( C2 )
!
!              W := C**T * V  =  (C1**T * V1 + C2**T * V2)  (stored in WORK)
!
!              W := C2**T
!
               DO 70 j = 1, k
                  CALL dcopy( n, c( m-k+j, 1 ), ldc, work( 1, j ), &
     &                        1 )
70             CONTINUE
!
!              W := W * V2
!
               CALL dtrmm( 'Right', 'Upper', 'No transpose', 'Unit', &
     &                     n, &
     &                     k, one, v( m-k+1, 1 ), ldv, work, ldwork )
               IF( m.GT.k ) THEN
!
!                 W := W + C1**T * V1
!
                  CALL dgemm( 'Transpose', 'No transpose', n, k, m-k, &
     &                        one, c, ldc, v, ldv, one, work, ldwork )
               END IF
!
!              W := W * T**T  or  W * T
!
               CALL dtrmm( 'Right', 'Lower', transt, 'Non-unit', n, &
     &                     k, &
     &                     one, t, ldt, work, ldwork )
!
!              C := C - V * W**T
!
               IF( m.GT.k ) THEN
!
!                 C1 := C1 - V1 * W**T
!
                  CALL dgemm( 'No transpose', 'Transpose', m-k, n, k, &
     &                        -one, v, ldv, work, ldwork, one, c, ldc )
               END IF
!
!              W := W * V2**T
!
               CALL dtrmm( 'Right', 'Upper', 'Transpose', 'Unit', n, &
     &                     k, &
     &                     one, v( m-k+1, 1 ), ldv, work, ldwork )
!
!              C2 := C2 - W**T
!
               DO 90 j = 1, k
                  DO 80 i = 1, n
                     c( m-k+j, i ) = c( m-k+j, i ) - work( i, j )
80                CONTINUE
90             CONTINUE
!
            ELSE IF( lsame( side, 'R' ) ) THEN
!
!              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
!
!              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
!
!              W := C2
!
               DO 100 j = 1, k
                  CALL dcopy( m, c( 1, n-k+j ), 1, work( 1, j ), 1 )
100            CONTINUE
!
!              W := W * V2
!
               CALL dtrmm( 'Right', 'Upper', 'No transpose', 'Unit', &
     &                     m, &
     &                     k, one, v( n-k+1, 1 ), ldv, work, ldwork )
               IF( n.GT.k ) THEN
!
!                 W := W + C1 * V1
!
                  CALL dgemm( 'No transpose', 'No transpose', m, k, &
     &                        n-k, &
     &                        one, c, ldc, v, ldv, one, work, ldwork )
               END IF
!
!              W := W * T  or  W * T**T
!
               CALL dtrmm( 'Right', 'Lower', trans, 'Non-unit', m, k, &
     &                     one, t, ldt, work, ldwork )
!
!              C := C - W * V**T
!
               IF( n.GT.k ) THEN
!
!                 C1 := C1 - W * V1**T
!
                  CALL dgemm( 'No transpose', 'Transpose', m, n-k, k, &
     &                        -one, work, ldwork, v, ldv, one, c, ldc )
               END IF
!
!              W := W * V2**T
!
               CALL dtrmm( 'Right', 'Upper', 'Transpose', 'Unit', m, &
     &                     k, &
     &                     one, v( n-k+1, 1 ), ldv, work, ldwork )
!
!              C2 := C2 - W
!
               DO 120 j = 1, k
                  DO 110 i = 1, m
                     c( i, n-k+j ) = c( i, n-k+j ) - work( i, j )
110               CONTINUE
120            CONTINUE
            END IF
         END IF
!
      ELSE IF( lsame( storev, 'R' ) ) THEN
!
         IF( lsame( direct, 'F' ) ) THEN
!
!           Let  V =  ( V1  V2 )    (V1: first K columns)
!           where  V1  is unit upper triangular.
!
            IF( lsame( side, 'L' ) ) THEN
!
!              Form  H * C  or  H**T * C  where  C = ( C1 )
!                                                    ( C2 )
!
!              W := C**T * V**T  =  (C1**T * V1**T + C2**T * V2**T) (stored in WORK)
!
!              W := C1**T
!
               DO 130 j = 1, k
                  CALL dcopy( n, c( j, 1 ), ldc, work( 1, j ), 1 )
130            CONTINUE
!
!              W := W * V1**T
!
               CALL dtrmm( 'Right', 'Upper', 'Transpose', 'Unit', n, &
     &                     k, &
     &                     one, v, ldv, work, ldwork )
               IF( m.GT.k ) THEN
!
!                 W := W + C2**T * V2**T
!
                  CALL dgemm( 'Transpose', 'Transpose', n, k, m-k, &
     &                        one, &
     &                        c( k+1, 1 ), ldc, v( 1, k+1 ), ldv, one, &
     &                        work, ldwork )
               END IF
!
!              W := W * T**T  or  W * T
!
               CALL dtrmm( 'Right', 'Upper', transt, 'Non-unit', n, &
     &                     k, &
     &                     one, t, ldt, work, ldwork )
!
!              C := C - V**T * W**T
!
               IF( m.GT.k ) THEN
!
!                 C2 := C2 - V2**T * W**T
!
                  CALL dgemm( 'Transpose', 'Transpose', m-k, n, k, &
     &                        -one, &
     &                        v( 1, k+1 ), ldv, work, ldwork, one, &
     &                        c( k+1, 1 ), ldc )
               END IF
!
!              W := W * V1
!
               CALL dtrmm( 'Right', 'Upper', 'No transpose', 'Unit', &
     &                     n, &
     &                     k, one, v, ldv, work, ldwork )
!
!              C1 := C1 - W**T
!
               DO 150 j = 1, k
                  DO 140 i = 1, n
                     c( j, i ) = c( j, i ) - work( i, j )
140               CONTINUE
150            CONTINUE
!
            ELSE IF( lsame( side, 'R' ) ) THEN
!
!              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
!
!              W := C * V**T  =  (C1*V1**T + C2*V2**T)  (stored in WORK)
!
!              W := C1
!
               DO 160 j = 1, k
                  CALL dcopy( m, c( 1, j ), 1, work( 1, j ), 1 )
160            CONTINUE
!
!              W := W * V1**T
!
               CALL dtrmm( 'Right', 'Upper', 'Transpose', 'Unit', m, &
     &                     k, &
     &                     one, v, ldv, work, ldwork )
               IF( n.GT.k ) THEN
!
!                 W := W + C2 * V2**T
!
                  CALL dgemm( 'No transpose', 'Transpose', m, k, n-k, &
     &                        one, c( 1, k+1 ), ldc, v( 1, k+1 ), ldv, &
     &                        one, work, ldwork )
               END IF
!
!              W := W * T  or  W * T**T
!
               CALL dtrmm( 'Right', 'Upper', trans, 'Non-unit', m, k, &
     &                     one, t, ldt, work, ldwork )
!
!              C := C - W * V
!
               IF( n.GT.k ) THEN
!
!                 C2 := C2 - W * V2
!
                  CALL dgemm( 'No transpose', 'No transpose', m, n-k, &
     &                        k, &
     &                        -one, work, ldwork, v( 1, k+1 ), ldv, one, &
     &                        c( 1, k+1 ), ldc )
               END IF
!
!              W := W * V1
!
               CALL dtrmm( 'Right', 'Upper', 'No transpose', 'Unit', &
     &                     m, &
     &                     k, one, v, ldv, work, ldwork )
!
!              C1 := C1 - W
!
               DO 180 j = 1, k
                  DO 170 i = 1, m
                     c( i, j ) = c( i, j ) - work( i, j )
170               CONTINUE
180            CONTINUE
!
            END IF
!
         ELSE
!
!           Let  V =  ( V1  V2 )    (V2: last K columns)
!           where  V2  is unit lower triangular.
!
            IF( lsame( side, 'L' ) ) THEN
!
!              Form  H * C  or  H**T * C  where  C = ( C1 )
!                                                    ( C2 )
!
!              W := C**T * V**T  =  (C1**T * V1**T + C2**T * V2**T) (stored in WORK)
!
!              W := C2**T
!
               DO 190 j = 1, k
                  CALL dcopy( n, c( m-k+j, 1 ), ldc, work( 1, j ), &
     &                        1 )
190            CONTINUE
!
!              W := W * V2**T
!
               CALL dtrmm( 'Right', 'Lower', 'Transpose', 'Unit', n, &
     &                     k, &
     &                     one, v( 1, m-k+1 ), ldv, work, ldwork )
               IF( m.GT.k ) THEN
!
!                 W := W + C1**T * V1**T
!
                  CALL dgemm( 'Transpose', 'Transpose', n, k, m-k, &
     &                        one, &
     &                        c, ldc, v, ldv, one, work, ldwork )
               END IF
!
!              W := W * T**T  or  W * T
!
               CALL dtrmm( 'Right', 'Lower', transt, 'Non-unit', n, &
     &                     k, &
     &                     one, t, ldt, work, ldwork )
!
!              C := C - V**T * W**T
!
               IF( m.GT.k ) THEN
!
!                 C1 := C1 - V1**T * W**T
!
                  CALL dgemm( 'Transpose', 'Transpose', m-k, n, k, &
     &                        -one, &
     &                        v, ldv, work, ldwork, one, c, ldc )
               END IF
!
!              W := W * V2
!
               CALL dtrmm( 'Right', 'Lower', 'No transpose', 'Unit', &
     &                     n, &
     &                     k, one, v( 1, m-k+1 ), ldv, work, ldwork )
!
!              C2 := C2 - W**T
!
               DO 210 j = 1, k
                  DO 200 i = 1, n
                     c( m-k+j, i ) = c( m-k+j, i ) - work( i, j )
200               CONTINUE
210            CONTINUE
!
            ELSE IF( lsame( side, 'R' ) ) THEN
!
!              Form  C * H  or  C * H'  where  C = ( C1  C2 )
!
!              W := C * V**T  =  (C1*V1**T + C2*V2**T)  (stored in WORK)
!
!              W := C2
!
               DO 220 j = 1, k
                  CALL dcopy( m, c( 1, n-k+j ), 1, work( 1, j ), 1 )
220            CONTINUE
!
!              W := W * V2**T
!
               CALL dtrmm( 'Right', 'Lower', 'Transpose', 'Unit', m, &
     &                     k, &
     &                     one, v( 1, n-k+1 ), ldv, work, ldwork )
               IF( n.GT.k ) THEN
!
!                 W := W + C1 * V1**T
!
                  CALL dgemm( 'No transpose', 'Transpose', m, k, n-k, &
     &                        one, c, ldc, v, ldv, one, work, ldwork )
               END IF
!
!              W := W * T  or  W * T**T
!
               CALL dtrmm( 'Right', 'Lower', trans, 'Non-unit', m, k, &
     &                     one, t, ldt, work, ldwork )
!
!              C := C - W * V
!
               IF( n.GT.k ) THEN
!
!                 C1 := C1 - W * V1
!
                  CALL dgemm( 'No transpose', 'No transpose', m, n-k, &
     &                        k, &
     &                        -one, work, ldwork, v, ldv, one, c, ldc )
               END IF
!
!              W := W * V2
!
               CALL dtrmm( 'Right', 'Lower', 'No transpose', 'Unit', &
     &                     m, &
     &                     k, one, v( 1, n-k+1 ), ldv, work, ldwork )
!
!              C1 := C1 - W
!
               DO 240 j = 1, k
                  DO 230 i = 1, m
                     c( i, n-k+j ) = c( i, n-k+j ) - work( i, j )
230               CONTINUE
240            CONTINUE
!
            END IF
!
         END IF
      END IF
!
      RETURN
!
!     End of DLARFB
!

      END
!> \brief \b DLAISNAN tests input for NaN by comparing two arguments for inequality.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download DLAISNAN + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaisnan.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaisnan.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaisnan.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       LOGICAL FUNCTION DLAISNAN( DIN1, DIN2 )
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION, INTENT(IN) :: DIN1, DIN2
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> This routine is not for general use.  It exists solely to avoid
!> over-optimization in DISNAN.
!>
!> DLAISNAN checks for NaNs by comparing its two arguments for
!> inequality.  NaN is the only floating-point value where NaN != NaN
!> returns .TRUE.  To check for NaNs, pass the same variable as both
!> arguments.
!>
!> A compiler must assume that the two arguments are
!> not the same variable, and the test will not be optimized away.
!> Interprocedural or whole-program optimization may delete this
!> test.  The ISNAN functions will be replaced by the correct
!> Fortran 03 intrinsic once the intrinsic is widely available.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] DIN1
!> \verbatim
!>          DIN1 is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] DIN2
!> \verbatim
!>          DIN2 is DOUBLE PRECISION
!>          Two numbers to compare for inequality.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup laisnan
!
!  =====================================================================

      LOGICAL FUNCTION dlaisnan( DIN1, DIN2 )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION, INTENT(IN) :: din1, din2
!     ..
!
!  =====================================================================
!
!  .. Executable Statements ..
      dlaisnan = (din1.NE.din2)
      RETURN

      END
!> \brief \b DLAQR1 sets a scalar multiple of the first column of the product of 2-by-2 or 3-by-3 matrix H and specified shifts.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download DLAQR1 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaqr1.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaqr1.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaqr1.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAQR1( N, H, LDH, SR1, SI1, SR2, SI2, V )
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION   SI1, SI2, SR1, SR2
!       INTEGER            LDH, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   H( LDH, * ), V( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>      Given a 2-by-2 or 3-by-3 matrix H, DLAQR1 sets v to a
!>      scalar multiple of the first column of the product
!>
!>      (*)  K = (H - (sr1 + i*si1)*I)*(H - (sr2 + i*si2)*I)
!>
!>      scaling to avoid overflows and most underflows. It
!>      is assumed that either
!>
!>              1) sr1 = sr2 and si1 = -si2
!>          or
!>              2) si1 = si2 = 0.
!>
!>      This is useful for starting double implicit shift bulges
!>      in the QR algorithm.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>              Order of the matrix H. N must be either 2 or 3.
!> \endverbatim
!>
!> \param[in] H
!> \verbatim
!>          H is DOUBLE PRECISION array, dimension (LDH,N)
!>              The 2-by-2 or 3-by-3 matrix H in (*).
!> \endverbatim
!>
!> \param[in] LDH
!> \verbatim
!>          LDH is INTEGER
!>              The leading dimension of H as declared in
!>              the calling procedure.  LDH >= N
!> \endverbatim
!>
!> \param[in] SR1
!> \verbatim
!>          SR1 is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] SI1
!> \verbatim
!>          SI1 is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] SR2
!> \verbatim
!>          SR2 is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] SI2
!> \verbatim
!>          SI2 is DOUBLE PRECISION
!>              The shifts in (*).
!> \endverbatim
!>
!> \param[out] V
!> \verbatim
!>          V is DOUBLE PRECISION array, dimension (N)
!>              A scalar multiple of the first column of the
!>              matrix K in (*).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup laqr1
!
!> \par Contributors:
!  ==================
!>
!>       Karen Braman and Ralph Byers, Department of Mathematics,
!>       University of Kansas, USA
!>
!  =====================================================================

      SUBROUTINE dlaqr1( N, H, LDH, SR1, SI1, SR2, SI2, V )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   SI1, SI2, SR1, SR2
      INTEGER            LDH, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   H( LDH, * ), V( * )
!     ..
!
!  ================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      parameter( zero = 0.0d0 )
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION   H21S, H31S, S
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          abs
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF( n.NE.2 .AND. n.NE.3 ) THEN
         RETURN
      END IF
!
      IF( n.EQ.2 ) THEN
         s = abs( h( 1, 1 )-sr2 ) + abs( si2 ) + abs( h( 2, 1 ) )
         IF( s.EQ.zero ) THEN
            v( 1 ) = zero
            v( 2 ) = zero
         ELSE
            h21s = h( 2, 1 ) / s
            v( 1 ) = h21s*h( 1, 2 ) + ( h( 1, 1 )-sr1 )* &
     &               ( ( h( 1, 1 )-sr2 ) / s ) - si1*( si2 / s )
            v( 2 ) = h21s*( h( 1, 1 )+h( 2, 2 )-sr1-sr2 )
         END IF
      ELSE
         s = abs( h( 1, 1 )-sr2 ) + abs( si2 ) + abs( h( 2, 1 ) ) + &
     &       abs( h( 3, 1 ) )
         IF( s.EQ.zero ) THEN
            v( 1 ) = zero
            v( 2 ) = zero
            v( 3 ) = zero
         ELSE
            h21s = h( 2, 1 ) / s
            h31s = h( 3, 1 ) / s
            v( 1 ) = ( h( 1, 1 )-sr1 )*( ( h( 1, 1 )-sr2 ) / s ) - &
     &               si1*( si2 / s ) + h( 1, 2 )*h21s + h( 1, 3 )*h31s
            v( 2 ) = h21s*( h( 1, 1 )+h( 2, 2 )-sr1-sr2 ) + &
     &               h( 2, 3 )*h31s
            v( 3 ) = h31s*( h( 1, 1 )+h( 3, 3 )-sr1-sr2 ) + &
     &               h21s*h( 3, 2 )
         END IF
      END IF

      END
!> \brief \b ILADLR scans a matrix for its last non-zero row.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download ILADLR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/iladlr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/iladlr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/iladlr.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       INTEGER FUNCTION ILADLR( M, N, A, LDA )
!
!       .. Scalar Arguments ..
!       INTEGER            M, N, LDA
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ILADLR scans A for its last non-zero row.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          The m by n matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >= max(1,M).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup ilalr
!
!  =====================================================================

      INTEGER FUNCTION iladlr( M, N, A, LDA )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            m, n, lda
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   a( lda, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION zero
      parameter( zero = 0.0d+0 )
!     ..
!     .. Local Scalars ..
      INTEGER i, j
!     ..
!     .. Executable Statements ..
!
!     Quick test for the common case where one corner is non-zero.
      IF( m.EQ.0 ) THEN
         iladlr = m
      ELSE IF( a(m, 1).NE.zero .OR. a(m, n).NE.zero ) THEN
         iladlr = m
      ELSE
!     Scan up each column tracking the last zero row seen.
         iladlr = 0
         DO j = 1, n
            i=m
            DO WHILE((a(max(i,1),j).EQ.zero).AND.(i.GE.1))
               i=i-1
            ENDDO
            iladlr = max( iladlr, i )
         END DO
      END IF
      RETURN

      END
!> \brief \b DTRMV
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,LDA,N
!       CHARACTER DIAG,TRANS,UPLO
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION A(LDA,*),X(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DTRMV  performs one of the matrix-vector operations
!>
!>    x := A*x,   or   x := A**T*x,
!>
!> where x is an n element vector and  A is an n by n unit, or non-unit,
!> upper or lower triangular matrix.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>           On entry, UPLO specifies whether the matrix is an upper or
!>           lower triangular matrix as follows:
!>
!>              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!>
!>              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>           On entry, TRANS specifies the operation to be performed as
!>           follows:
!>
!>              TRANS = 'N' or 'n'   x := A*x.
!>
!>              TRANS = 'T' or 't'   x := A**T*x.
!>
!>              TRANS = 'C' or 'c'   x := A**T*x.
!> \endverbatim
!>
!> \param[in] DIAG
!> \verbatim
!>          DIAG is CHARACTER*1
!>           On entry, DIAG specifies whether or not A is unit
!>           triangular as follows:
!>
!>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!>
!>              DIAG = 'N' or 'n'   A is not assumed to be unit
!>                                  triangular.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry, N specifies the order of the matrix A.
!>           N must be at least zero.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension ( LDA, N )
!>           Before entry with  UPLO = 'U' or 'u', the leading n by n
!>           upper triangular part of the array A must contain the upper
!>           triangular matrix and the strictly lower triangular part of
!>           A is not referenced.
!>           Before entry with UPLO = 'L' or 'l', the leading n by n
!>           lower triangular part of the array A must contain the lower
!>           triangular matrix and the strictly upper triangular part of
!>           A is not referenced.
!>           Note that when  DIAG = 'U' or 'u', the diagonal elements of
!>           A are not referenced either, but are assumed to be unity.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in the calling (sub) program. LDA must be at least
!>           max( 1, n ).
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension at least
!>           ( 1 + ( n - 1 )*abs( INCX ) ).
!>           Before entry, the incremented array X must contain the n
!>           element vector x. On exit, X is overwritten with the
!>           transformed vector x.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>           On entry, INCX specifies the increment for the elements of
!>           X. INCX must not be zero.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup trmv
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Level 2 Blas routine.
!>  The vector and matrix arguments are not referenced when N = 0, or M = 0
!>
!>  -- Written on 22-October-1986.
!>     Jack Dongarra, Argonne National Lab.
!>     Jeremy Du Croz, Nag Central Office.
!>     Sven Hammarling, Nag Central Office.
!>     Richard Hanson, Sandia National Labs.
!> \endverbatim
!>
!  =====================================================================

      SUBROUTINE dtrmv(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
!
!  -- Reference BLAS level2 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER INCX,LDA,N
      CHARACTER DIAG,TRANS,UPLO
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),X(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO
      parameter(zero=0.0d+0)
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,J,JX,KX
      LOGICAL NOUNIT
!     ..
!     .. External Functions ..
!      LOGICAL LSAME
!      EXTERNAL lsame
!     ..
!     .. External Subroutines ..
!      EXTERNAL xerbla
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC max
!     ..
!
!     Test the input parameters.
!
      info = 0
      IF (.NOT.lsame(uplo,'U') .AND. .NOT.lsame(uplo,'L')) THEN
          info = 1
      ELSE IF (.NOT.lsame(trans,'N') .AND. &
     &         .NOT.lsame(trans,'T') .AND. &
     &         .NOT.lsame(trans,'C')) THEN
          info = 2
      ELSE IF (.NOT.lsame(diag,'U') .AND. &
     &         .NOT.lsame(diag,'N')) THEN
          info = 3
      ELSE IF (n.LT.0) THEN
          info = 4
      ELSE IF (lda.LT.max(1,n)) THEN
          info = 6
      ELSE IF (incx.EQ.0) THEN
          info = 8
      END IF
      IF (info.NE.0) THEN
          CALL xerbla('DTRMV ',info)
          RETURN
      END IF
!
!     Quick return if possible.
!
      IF (n.EQ.0) RETURN
!
      nounit = lsame(diag,'N')
!
!     Set up the start point in X if the increment is not unity. This
!     will be  ( N - 1 )*INCX  too small for descending loops.
!
      IF (incx.LE.0) THEN
          kx = 1 - (n-1)*incx
      ELSE IF (incx.NE.1) THEN
          kx = 1
      END IF
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
      IF (lsame(trans,'N')) THEN
!
!        Form  x := A*x.
!
          IF (lsame(uplo,'U')) THEN
              IF (incx.EQ.1) THEN
                  DO 20 j = 1,n
                      IF (x(j).NE.zero) THEN
                          temp = x(j)
                          DO 10 i = 1,j - 1
                              x(i) = x(i) + temp*a(i,j)
10                        CONTINUE
                          IF (nounit) x(j) = x(j)*a(j,j)
                      END IF
20                CONTINUE
              ELSE
                  jx = kx
                  DO 40 j = 1,n
                      IF (x(jx).NE.zero) THEN
                          temp = x(jx)
                          ix = kx
                          DO 30 i = 1,j - 1
                              x(ix) = x(ix) + temp*a(i,j)
                              ix = ix + incx
30                        CONTINUE
                          IF (nounit) x(jx) = x(jx)*a(j,j)
                      END IF
                      jx = jx + incx
40                CONTINUE
              END IF
          ELSE
              IF (incx.EQ.1) THEN
                  DO 60 j = n,1,-1
                      IF (x(j).NE.zero) THEN
                          temp = x(j)
                          DO 50 i = n,j + 1,-1
                              x(i) = x(i) + temp*a(i,j)
50                        CONTINUE
                          IF (nounit) x(j) = x(j)*a(j,j)
                      END IF
60                CONTINUE
              ELSE
                  kx = kx + (n-1)*incx
                  jx = kx
                  DO 80 j = n,1,-1
                      IF (x(jx).NE.zero) THEN
                          temp = x(jx)
                          ix = kx
                          DO 70 i = n,j + 1,-1
                              x(ix) = x(ix) + temp*a(i,j)
                              ix = ix - incx
70                        CONTINUE
                          IF (nounit) x(jx) = x(jx)*a(j,j)
                      END IF
                      jx = jx - incx
80                CONTINUE
              END IF
          END IF
      ELSE
!
!        Form  x := A**T*x.
!
          IF (lsame(uplo,'U')) THEN
              IF (incx.EQ.1) THEN
                  DO 100 j = n,1,-1
                      temp = x(j)
                      IF (nounit) temp = temp*a(j,j)
                      DO 90 i = j - 1,1,-1
                          temp = temp + a(i,j)*x(i)
90                    CONTINUE
                      x(j) = temp
100               CONTINUE
              ELSE
                  jx = kx + (n-1)*incx
                  DO 120 j = n,1,-1
                      temp = x(jx)
                      ix = jx
                      IF (nounit) temp = temp*a(j,j)
                      DO 110 i = j - 1,1,-1
                          ix = ix - incx
                          temp = temp + a(i,j)*x(ix)
110                   CONTINUE
                      x(jx) = temp
                      jx = jx - incx
120               CONTINUE
              END IF
          ELSE
              IF (incx.EQ.1) THEN
                  DO 140 j = 1,n
                      temp = x(j)
                      IF (nounit) temp = temp*a(j,j)
                      DO 130 i = j + 1,n
                          temp = temp + a(i,j)*x(i)
130                   CONTINUE
                      x(j) = temp
140               CONTINUE
              ELSE
                  jx = kx
                  DO 160 j = 1,n
                      temp = x(jx)
                      ix = jx
                      IF (nounit) temp = temp*a(j,j)
                      DO 150 i = j + 1,n
                          ix = ix + incx
                          temp = temp + a(i,j)*x(ix)
150                   CONTINUE
                      x(jx) = temp
                      jx = jx + incx
160               CONTINUE
              END IF
          END IF
      END IF
!
      RETURN
!
!     End of DTRMV
!

      END

!> \brief \b DGEHD2 reduces a general square matrix to upper Hessenberg form using an unblocked algorithm.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download DGEHD2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgehd2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgehd2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgehd2.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGEHD2( N, ILO, IHI, A, LDA, TAU, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            IHI, ILO, INFO, LDA, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGEHD2 reduces a real general matrix A to upper Hessenberg form H by
!> an orthogonal similarity transformation:  Q**T * A * Q = H .
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] ILO
!> \verbatim
!>          ILO is INTEGER
!> \endverbatim
!>
!> \param[in] IHI
!> \verbatim
!>          IHI is INTEGER
!>
!>          It is assumed that A is already upper triangular in rows
!>          and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally
!>          set by a previous call to DGEBAL; otherwise they should be
!>          set to 1 and N respectively. See Further Details.
!>          1 <= ILO <= IHI <= max(1,N).
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the n by n general matrix to be reduced.
!>          On exit, the upper triangle and the first subdiagonal of A
!>          are overwritten with the upper Hessenberg matrix H, and the
!>          elements below the first subdiagonal, with the array TAU,
!>          represent the orthogonal matrix Q as a product of elementary
!>          reflectors. See Further Details.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION array, dimension (N-1)
!>          The scalar factors of the elementary reflectors (see Further
!>          Details).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup gehd2
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The matrix Q is represented as a product of (ihi-ilo) elementary
!>  reflectors
!>
!>     Q = H(ilo) H(ilo+1) . . . H(ihi-1).
!>
!>  Each H(i) has the form
!>
!>     H(i) = I - tau * v * v**T
!>
!>  where tau is a real scalar, and v is a real vector with
!>  v(1:i) = 0, v(i+1) = 1 and v(ihi+1:n) = 0; v(i+2:ihi) is stored on
!>  exit in A(i+2:ihi,i), and tau in TAU(i).
!>
!>  The contents of A are illustrated by the following example, with
!>  n = 7, ilo = 2 and ihi = 6:
!>
!>  on entry,                        on exit,
!>
!>  ( a   a   a   a   a   a   a )    (  a   a   h   h   h   h   a )
!>  (     a   a   a   a   a   a )    (      a   h   h   h   h   a )
!>  (     a   a   a   a   a   a )    (      h   h   h   h   h   h )
!>  (     a   a   a   a   a   a )    (      v2  h   h   h   h   h )
!>  (     a   a   a   a   a   a )    (      v2  v3  h   h   h   h )
!>  (     a   a   a   a   a   a )    (      v2  v3  v4  h   h   h )
!>  (                         a )    (                          a )
!>
!>  where a denotes an element of the original matrix A, h denotes a
!>  modified element of the upper Hessenberg matrix H, and vi denotes an
!>  element of the vector defining H(i).
!> \endverbatim
!>
!  =====================================================================

      SUBROUTINE dgehd2( N, ILO, IHI, A, LDA, TAU, WORK, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            IHI, ILO, INFO, LDA, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      parameter( one = 1.0d+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I
!     ..
!     .. External Subroutines ..
!      EXTERNAL           dlarf1f, dlarfg, xerbla
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          max, min
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      info = 0
      IF( n.LT.0 ) THEN
         info = -1
      ELSE IF( ilo.LT.1 .OR. ilo.GT.max( 1, n ) ) THEN
         info = -2
      ELSE IF( ihi.LT.min( ilo, n ) .OR. ihi.GT.n ) THEN
         info = -3
      ELSE IF( lda.LT.max( 1, n ) ) THEN
         info = -5
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DGEHD2', -info )
         RETURN
      END IF
!
      DO 10 i = ilo, ihi - 1
!
!        Compute elementary reflector H(i) to annihilate A(i+2:ihi,i)
!
         CALL dlarfg( ihi-i, a( i+1, i ), a( min( i+2, n ), i ), 1, &
     &                tau( i ) )
!
!        Apply H(i) to A(1:ihi,i+1:ihi) from the right
!
         CALL dlarf1f( 'Right', ihi, ihi-i, a( i+1, i ), 1, tau( i ), &
     &               a( 1, i+1 ), lda, work )
!
!        Apply H(i) to A(i+1:ihi,i+1:n) from the left
!
         CALL dlarf1f( 'Left', ihi-i, n-i, a( i+1, i ), 1, tau( i ), &
     &               a( i+1, i+1 ), lda, work )
!
10    CONTINUE
!
      RETURN
!
!     End of DGEHD2
!

      END
!> \brief \b DSWAP
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSWAP(N,DX,INCX,DY,INCY)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,INCY,N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION DX(*),DY(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    DSWAP interchanges two vectors.
!>    uses unrolled loops for increments equal to 1.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         number of elements in input vector(s)
!> \endverbatim
!>
!> \param[in,out] DX
!> \verbatim
!>          DX is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>         storage spacing between elements of DX
!> \endverbatim
!>
!> \param[in,out] DY
!> \verbatim
!>          DY is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>         storage spacing between elements of DY
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup swap
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, linpack, 3/11/78.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================

      SUBROUTINE dswap(N,DX,INCX,DY,INCY)
!
!  -- Reference BLAS level1 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER I,IX,IY,M,MP1
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC mod
!     ..
      IF (n.LE.0) RETURN
      IF (incx.EQ.1 .AND. incy.EQ.1) THEN
!
!       code for both increments equal to 1
!
!
!       clean-up loop
!
         m = mod(n,3)
         IF (m.NE.0) THEN
            DO i = 1,m
               dtemp = dx(i)
               dx(i) = dy(i)
               dy(i) = dtemp
            END DO
            IF (n.LT.3) RETURN
         END IF
         mp1 = m + 1
         DO i = mp1,n,3
            dtemp = dx(i)
            dx(i) = dy(i)
            dy(i) = dtemp
            dtemp = dx(i+1)
            dx(i+1) = dy(i+1)
            dy(i+1) = dtemp
            dtemp = dx(i+2)
            dx(i+2) = dy(i+2)
            dy(i+2) = dtemp
         END DO
      ELSE
!
!       code for unequal increments or equal increments not equal
!         to 1
!
         ix = 1
         iy = 1
         IF (incx.LT.0) ix = (-n+1)*incx + 1
         IF (incy.LT.0) iy = (-n+1)*incy + 1
         DO i = 1,n
            dtemp = dx(ix)
            dx(ix) = dy(iy)
            dy(iy) = dtemp
            ix = ix + incx
            iy = iy + incy
         END DO
      END IF
      RETURN
!
!     End of DSWAP
!

      END
!> \brief \b DLANV2 computes the Schur factorization of a real 2-by-2 nonsymmetric matrix in standard form.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download DLANV2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlanv2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlanv2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlanv2.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLANV2( A, B, C, D, RT1R, RT1I, RT2R, RT2I, CS, SN )
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION   A, B, C, CS, D, RT1I, RT1R, RT2I, RT2R, SN
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLANV2 computes the Schur factorization of a real 2-by-2 nonsymmetric
!> matrix in standard form:
!>
!>      [ A  B ] = [ CS -SN ] [ AA  BB ] [ CS  SN ]
!>      [ C  D ]   [ SN  CS ] [ CC  DD ] [-SN  CS ]
!>
!> where either
!> 1) CC = 0 so that AA and DD are real eigenvalues of the matrix, or
!> 2) AA = DD and BB*CC < 0, so that AA + or - sqrt(BB*CC) are complex
!> conjugate eigenvalues.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is DOUBLE PRECISION
!>          On entry, the elements of the input matrix.
!>          On exit, they are overwritten by the elements of the
!>          standardised Schur form.
!> \endverbatim
!>
!> \param[out] RT1R
!> \verbatim
!>          RT1R is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[out] RT1I
!> \verbatim
!>          RT1I is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[out] RT2R
!> \verbatim
!>          RT2R is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[out] RT2I
!> \verbatim
!>          RT2I is DOUBLE PRECISION
!>          The real and imaginary parts of the eigenvalues. If the
!>          eigenvalues are a complex conjugate pair, RT1I > 0.
!> \endverbatim
!>
!> \param[out] CS
!> \verbatim
!>          CS is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[out] SN
!> \verbatim
!>          SN is DOUBLE PRECISION
!>          Parameters of the rotation matrix.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup lanv2
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Modified by V. Sima, Research Institute for Informatics, Bucharest,
!>  Romania, to reduce the risk of cancellation errors,
!>  when computing real eigenvalues, and to ensure, if possible, that
!>  abs(RT1R) >= abs(RT2R).
!> \endverbatim
!>
!  =====================================================================

      SUBROUTINE dlanv2( A, B, C, D, RT1R, RT1I, RT2R, RT2I, CS, SN )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, CS, D, RT1I, RT1R, RT2I, RT2R, SN
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, HALF, ONE, TWO
      parameter( zero = 0.0d+0, half = 0.5d+0, one = 1.0d+0, &
     &                     two = 2.0d0 )
      DOUBLE PRECISION   MULTPL
      parameter( multpl = 4.0d+0 )
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION   AA, BB, BCMAX, BCMIS, CC, CS1, DD, EPS, P, SAB, &
     &                   SAC, SCALE, SIGMA, SN1, TAU, TEMP, Z, SAFMIN, &
     &                   SAFMN2, SAFMX2
      INTEGER            COUNT
!     ..
!     .. External Functions ..
!      DOUBLE PRECISION   DLAMCH, DLAPY2
!      EXTERNAL           dlamch, dlapy2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          abs, max, min, sign, sqrt
!     ..
!     .. Executable Statements ..
!
      safmin = dlamch( 'S' )
      eps = dlamch( 'P' )
      safmn2 = dlamch( 'B' )**int( log( safmin / eps ) / &
     &            log( dlamch( 'B' ) ) / two )
      safmx2 = one / safmn2
      IF( c.EQ.zero ) THEN
         cs = one
         sn = zero
!
      ELSE IF( b.EQ.zero ) THEN
!
!        Swap rows and columns
!
         cs = zero
         sn = one
         temp = d
         d = a
         a = temp
         b = -c
         c = zero
!
      ELSE IF( ( a-d ).EQ.zero .AND. sign( one, b ).NE.sign( one, c ) ) &
     &          THEN
         cs = one
         sn = zero
!
      ELSE
!
         temp = a - d
         p = half*temp
         bcmax = max( abs( b ), abs( c ) )
         bcmis = min( abs( b ), abs( c ) )*sign( one, b )*sign( one, c )
         scale = max( abs( p ), bcmax )
         z = ( p / scale )*p + ( bcmax / scale )*bcmis
!
!        If Z is of the order of the machine accuracy, postpone the
!        decision on the nature of eigenvalues
!
         IF( z.GE.multpl*eps ) THEN
!
!           Real eigenvalues. Compute A and D.
!
            z = p + sign( sqrt( scale )*sqrt( z ), p )
            a = d + z
            d = d - ( bcmax / z )*bcmis
!
!           Compute B and the rotation matrix
!
            tau = dlapy2( c, z )
            cs = z / tau
            sn = c / tau
            b = b - c
            c = zero
!
         ELSE
!
!           Complex eigenvalues, or real (almost) equal eigenvalues.
!           Make diagonal elements equal.
!
            count = 0
            sigma = b + c
10          CONTINUE
            count = count + 1
            scale = max( abs(temp), abs(sigma) )
            IF( scale.GE.safmx2 ) THEN
               sigma = sigma * safmn2
               temp = temp * safmn2
               IF (count .LE. 20) &
     &            GOTO 10
            END IF
            IF( scale.LE.safmn2 ) THEN
               sigma = sigma * safmx2
               temp = temp * safmx2
               IF (count .LE. 20) &
     &            GOTO 10
            END IF
            p = half*temp
            tau = dlapy2( sigma, temp )
            cs = sqrt( half*( one+abs( sigma ) / tau ) )
            sn = -( p / ( tau*cs ) )*sign( one, sigma )
!
!           Compute [ AA  BB ] = [ A  B ] [ CS -SN ]
!                   [ CC  DD ]   [ C  D ] [ SN  CS ]
!
            aa = a*cs + b*sn
            bb = -a*sn + b*cs
            cc = c*cs + d*sn
            dd = -c*sn + d*cs
!
!           Compute [ A  B ] = [ CS  SN ] [ AA  BB ]
!                   [ C  D ]   [-SN  CS ] [ CC  DD ]
!
!           Note: Some of the multiplications are wrapped in parentheses to
!                 prevent compilers from using FMA instructions. See
!                 https://github.com/Reference-LAPACK/lapack/issues/1031.
!
            a = aa*cs + cc*sn
            b = ( bb*cs ) + ( dd*sn )
            c = -( aa*sn ) + ( cc*cs )
            d = -bb*sn + dd*cs
!
            temp = half*( a+d )
            a = temp
            d = temp
!
            IF( c.NE.zero ) THEN
               IF( b.NE.zero ) THEN
                  IF( sign( one, b ).EQ.sign( one, c ) ) THEN
!
!                    Real eigenvalues: reduce to upper triangular form
!
                     sab = sqrt( abs( b ) )
                     sac = sqrt( abs( c ) )
                     p = sign( sab*sac, c )
                     tau = one / sqrt( abs( b+c ) )
                     a = temp + p
                     d = temp - p
                     b = b - c
                     c = zero
                     cs1 = sab*tau
                     sn1 = sac*tau
                     temp = cs*cs1 - sn*sn1
                     sn = cs*sn1 + sn*cs1
                     cs = temp
                  END IF
               ELSE
                  b = -c
                  c = zero
                  temp = cs
                  cs = -sn
                  sn = temp
               END IF
            END IF
         END IF
!
      END IF
!
!     Store eigenvalues in (RT1R,RT1I) and (RT2R,RT2I).
!
      rt1r = a
      rt2r = d
      IF( c.EQ.zero ) THEN
         rt1i = zero
         rt2i = zero
      ELSE
         rt1i = sqrt( abs( b ) )*sqrt( abs( c ) )
         rt2i = -rt1i
      END IF
      RETURN
!
!     End of DLANV2
!

      END
!> \brief \b DLADIV performs complex division in real arithmetic, avoiding unnecessary overflow.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download DLADIV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dladiv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dladiv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dladiv.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLADIV( A, B, C, D, P, Q )
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION   A, B, C, D, P, Q
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLADIV performs complex division in  real arithmetic
!>
!>                       a + i*b
!>            p + i*q = ---------
!>                       c + i*d
!>
!> The algorithm is due to Michael Baudin and Robert L. Smith
!> and can be found in the paper
!> "A Robust Complex Division in Scilab"
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is DOUBLE PRECISION
!>          The scalars a, b, c, and d in the above expression.
!> \endverbatim
!>
!> \param[out] P
!> \verbatim
!>          P is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is DOUBLE PRECISION
!>          The scalars p and q in the above expression.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup ladiv
!
!  =====================================================================

      SUBROUTINE dladiv( A, B, C, D, P, Q )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, D, P, Q
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   BS
      parameter( bs = 2.0d0 )
      DOUBLE PRECISION   HALF
      parameter( half = 0.5d0 )
      DOUBLE PRECISION   TWO
      parameter( two = 2.0d0 )
!
!     .. Local Scalars ..
      DOUBLE PRECISION   AA, BB, CC, DD, AB, CD, S, OV, UN, BE, EPS
!     ..
!     .. External Functions ..
!      DOUBLE PRECISION   DLAMCH
!      EXTERNAL           dlamch
!     ..
!     .. External Subroutines ..
!      EXTERNAL           dladiv1
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          abs, max
!     ..
!     .. Executable Statements ..
!
      aa = a
      bb = b
      cc = c
      dd = d
      ab = max( abs(a), abs(b) )
      cd = max( abs(c), abs(d) )
      s = 1.0d0

      ov = dlamch( 'Overflow threshold' )
      un = dlamch( 'Safe minimum' )
      eps = dlamch( 'Epsilon' )
      be = bs / (eps*eps)

      IF( ab >= half*ov ) THEN
         aa = half * aa
         bb = half * bb
         s  = two * s
      END IF
      IF( cd >= half*ov ) THEN
         cc = half * cc
         dd = half * dd
         s  = half * s
      END IF
      IF( ab <= un*bs/eps ) THEN
         aa = aa * be
         bb = bb * be
         s  = s / be
      END IF
      IF( cd <= un*bs/eps ) THEN
         cc = cc * be
         dd = dd * be
         s  = s * be
      END IF
      IF( abs( d ).LE.abs( c ) ) THEN
         CALL dladiv1(aa, bb, cc, dd, p, q)
      ELSE
         CALL dladiv1(bb, aa, dd, cc, p, q)
         q = -q
      END IF
      p = p * s
      q = q * s
!
      RETURN
!
!     End of DLADIV
!

      END

!> \ingroup ladiv

      SUBROUTINE dladiv1( A, B, C, D, P, Q )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, D, P, Q
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      parameter( one = 1.0d0 )
!
!     .. Local Scalars ..
      DOUBLE PRECISION   R, T
!     ..
!     .. External Functions ..
!      DOUBLE PRECISION   DLADIV2
!      EXTERNAL           dladiv2
!     ..
!     .. Executable Statements ..
!
      r = d / c
      t = one / (c + d * r)
      p = dladiv2(a, b, c, d, r, t)
      a = -a
      q = dladiv2(b, a, c, d, r, t)
!
      RETURN
!
!     End of DLADIV1
!

      END

!> \ingroup ladiv

      DOUBLE PRECISION FUNCTION dladiv2( A, B, C, D, R, T )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   a, b, c, d, r, t
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   zero
      parameter( zero = 0.0d0 )
!
!     .. Local Scalars ..
      DOUBLE PRECISION   br
!     ..
!     .. Executable Statements ..
!
      IF( r.NE.zero ) THEN
         br = b * r
         IF( br.NE.zero ) THEN
            dladiv2 = (a + br) * t
         ELSE
            dladiv2 = a * t + (b * t) * r
         END IF
      ELSE
         dladiv2 = (a + d * (b / c)) * t
      END IF
!
      RETURN
!
!     End of DLADIV2
!

      END
!> \brief \b DLARF1L applies an elementary reflector to a general rectangular
!              matrix assuming v(lastv) = 1 where lastv is the last non-zero
!              element
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download DLARF1L + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarf1l.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarf1l.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarf1l.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLARF1L( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          SIDE
!       INTEGER            INCV, LDC, M, N
!       DOUBLE PRECISION   TAU
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   C( LDC, * ), V( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLARF1L applies a real elementary reflector H to a real m by n matrix
!> C, from either the left or the right. H is represented in the form
!>
!>       H = I - tau * v * v**T
!>
!> where tau is a real scalar and v is a real vector.
!>
!> If tau = 0, then H is taken to be the unit matrix.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'L': form  H * C
!>          = 'R': form  C * H
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix C.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix C.
!> \endverbatim
!>
!> \param[in] V
!> \verbatim
!>          V is DOUBLE PRECISION array, dimension
!>                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
!>                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
!>          The vector v in the representation of H. V is not used if
!>          TAU = 0.
!> \endverbatim
!>
!> \param[in] INCV
!> \verbatim
!>          INCV is INTEGER
!>          The increment between elements of v. INCV <> 0.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION
!>          The value tau in the representation of H.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension (LDC,N)
!>          On entry, the m by n matrix C.
!>          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
!>          or C * H if SIDE = 'R'.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C. LDC >= max(1,M).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension
!>                         (N) if SIDE = 'L'
!>                      or (M) if SIDE = 'R'
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup larf
!
!  =====================================================================

      SUBROUTINE dlarf1l( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            INCV, LDC, M, N
      DOUBLE PRECISION   TAU
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   C( LDC, * ), V( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      parameter( one = 1.0d+0, zero = 0.0d+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            APPLYLEFT
      INTEGER            I, FIRSTV, LASTV, LASTC
!     ..
!     .. External Subroutines ..
!      EXTERNAL           dgemv, dger
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      INTEGER            ILADLR, ILADLC
!      EXTERNAL           lsame, iladlr, iladlc
!     ..
!     .. Executable Statements ..
!
      applyleft = lsame( side, 'L' )
      firstv = 1
      lastc = 0
      IF( tau.NE.zero ) THEN
!     Set up variables for scanning V.  LASTV begins pointing to the end
!     of V.
         IF( applyleft ) THEN
            lastv = m
         ELSE
            lastv = n
         END IF
         i = 1
!     Look for the last non-zero row in V.
         DO WHILE( lastv.GT.firstv .AND. v( i ).EQ.zero )
            firstv = firstv + 1
            i = i + incv
         END DO
         IF( applyleft ) THEN
!     Scan for the last non-zero column in C(1:lastv,:).
            lastc = iladlc(lastv, n, c, ldc)
         ELSE
!     Scan for the last non-zero row in C(:,1:lastv).
            lastc = iladlr(m, lastv, c, ldc)
         END IF
      END IF
      IF( lastc.EQ.0 ) THEN
         RETURN
      END IF
      IF( applyleft ) THEN
!
!        Form  H * C
!
         IF( lastv.GT.0 ) THEN
            ! Check if m = 1. This means v = 1, So we just need to compu
            ! C := HC = (1-\tau)C.
            IF( lastv.EQ.firstv ) THEN
               CALL dscal(lastc, one - tau, c( firstv, 1), ldc)
            ELSE
!
!              w(1:lastc,1) := C(1:lastv,1:lastc)**T * v(1:lastv,1)
!
               ! w(1:lastc,1) := C(1:lastv-1,1:lastc)**T * v(1:lastv-1,1
               CALL dgemv( 'Transpose', lastv-firstv, lastc, one, &
     &                     c(firstv,1), ldc, v(i), incv, zero, &
     &                     work, 1)
               ! w(1:lastc,1) += C(lastv,1:lastc)**T * v(lastv,1) = C(la
               CALL daxpy(lastc, one, c(lastv,1), ldc, work, 1)
!
!              C(1:lastv,1:lastc) := C(...) - tau * v(1:lastv,1) * w(1:lastc,1)**T
!
               ! C(lastv, 1:lastc)   := C(...) - tau * v(lastv,1) * w(1:
               !                      = C(...) - tau * w(1:lastc,1)**T
               CALL daxpy(lastc, -tau, work, 1, c(lastv,1), ldc)
               ! C(1:lastv-1,1:lastc) := C(...) - tau * v(1:lastv-1,1)*w
               CALL dger(lastv-firstv, lastc, -tau, v(i), incv, &
     &                   work, 1, c(firstv,1), ldc)
            END IF
         END IF
      ELSE
!
!        Form  C * H
!
         IF( lastv.GT.0 ) THEN
            ! Check if n = 1. This means v = 1, so we just need to compu
            ! C := CH = C(1-\tau).
            IF( lastv.EQ.firstv ) THEN
               CALL dscal(lastc, one - tau, c, 1)
            ELSE
!
!              w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1)
!
               ! w(1:lastc,1) := C(1:lastc,1:lastv-1) * v(1:lastv-1,1)
               CALL dgemv( 'No transpose', lastc, lastv-firstv, &
     &            one, c(1,firstv), ldc, v(i), incv, zero, work, 1 )
               ! w(1:lastc,1) += C(1:lastc,lastv) * v(lastv,1) = C(1:las
               CALL daxpy(lastc, one, c(1,lastv), 1, work, 1)
!
!              C(1:lastc,1:lastv) := C(...) - tau * w(1:lastc,1) * v(1:lastv,1)**T
!
               ! C(1:lastc,lastv)     := C(...) - tau * w(1:lastc,1) * v
               !                       = C(...) - tau * w(1:lastc,1)
               CALL daxpy(lastc, -tau, work, 1, c(1,lastv), 1)
               ! C(1:lastc,1:lastv-1) := C(...) - tau * w(1:lastc,1) * v
               CALL dger( lastc, lastv-firstv, -tau, work, 1, v(i), &
     &                     incv, c(1,firstv), ldc )
            END IF
         END IF
      END IF
      RETURN
!
!     End of DLARF1L
!

      END
!> \brief \b DLAMCH
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!      DOUBLE PRECISION FUNCTION DLAMCH( CMACH )
!
!     .. Scalar Arguments ..
!     CHARACTER          CMACH
!     ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAMCH determines double precision machine parameters.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] CMACH
!> \verbatim
!>          CMACH is CHARACTER*1
!>          Specifies the value to be returned by DLAMCH:
!>          = 'E' or 'e',   DLAMCH := eps
!>          = 'S' or 's ,   DLAMCH := sfmin
!>          = 'B' or 'b',   DLAMCH := base
!>          = 'P' or 'p',   DLAMCH := eps*base
!>          = 'N' or 'n',   DLAMCH := t
!>          = 'R' or 'r',   DLAMCH := rnd
!>          = 'M' or 'm',   DLAMCH := emin
!>          = 'U' or 'u',   DLAMCH := rmin
!>          = 'L' or 'l',   DLAMCH := emax
!>          = 'O' or 'o',   DLAMCH := rmax
!>          where
!>          eps   = relative machine precision
!>          sfmin = safe minimum, such that 1/sfmin does not overflow
!>          base  = base of the machine
!>          prec  = eps*base
!>          t     = number of (base) digits in the mantissa
!>          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
!>          emin  = minimum exponent before (gradual) underflow
!>          rmin  = underflow threshold - base**(emin-1)
!>          emax  = largest exponent before overflow
!>          rmax  = overflow threshold  - (base**emax)*(1-eps)
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!

!> \ingroup lamch
!
!  =====================================================================

      DOUBLE PRECISION FUNCTION dlamch( CMACH )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          cmach
!     ..
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   one, zero
      parameter( one = 1.0d+0, zero = 0.0d+0 )
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION   rnd, eps, sfmin, small, rmach
!     ..
!     .. External Functions ..
!      LOGICAL            lsame
!      EXTERNAL           lsame
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          digits, epsilon, huge, maxexponent, &
     &                   minexponent, radix, tiny
!     ..
!     .. Executable Statements ..
!
!
!     Assume rounding, not chopping. Always.
!
      rnd = one
!
      IF( one.EQ.rnd ) THEN
         eps = epsilon(zero) * 0.5
      ELSE
         eps = epsilon(zero)
      END IF
!
      IF( lsame( cmach, 'E' ) ) THEN
         rmach = eps
      ELSE IF( lsame( cmach, 'S' ) ) THEN
         sfmin = tiny(zero)
         small = one / huge(zero)
         IF( small.GE.sfmin ) THEN
!
!           Use SMALL plus a bit, to avoid the possibility of rounding
!           causing overflow when computing  1/sfmin.
!
            sfmin = small*( one+eps )
         END IF
         rmach = sfmin
      ELSE IF( lsame( cmach, 'B' ) ) THEN
         rmach = radix(zero)
      ELSE IF( lsame( cmach, 'P' ) ) THEN
         rmach = eps * radix(zero)
      ELSE IF( lsame( cmach, 'N' ) ) THEN
         rmach = digits(zero)
      ELSE IF( lsame( cmach, 'R' ) ) THEN
         rmach = rnd
      ELSE IF( lsame( cmach, 'M' ) ) THEN
         rmach = minexponent(zero)
      ELSE IF( lsame( cmach, 'U' ) ) THEN
         rmach = tiny(zero)
      ELSE IF( lsame( cmach, 'L' ) ) THEN
         rmach = maxexponent(zero)
      ELSE IF( lsame( cmach, 'O' ) ) THEN
         rmach = huge(zero)
      ELSE
         rmach = zero
      END IF
!
      dlamch = rmach
      RETURN
!
!     End of DLAMCH
!

      END

!***********************************************************************
!> \brief \b DLAMC3
!> \details
!> \b Purpose:
!> \verbatim
!> DLAMC3  is intended to force  A  and  B  to be stored prior to doing
!> the addition of  A  and  B ,  for use in situations where optimizers
!> might hold one of these in a register.
!> \endverbatim
!> \author LAPACK is a software package provided by Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..
!> \param[in] A
!> \verbatim
!>          A is a DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is a DOUBLE PRECISION
!>          The values A and B.
!> \endverbatim
!>
!> \ingroup lamc3
!>

      DOUBLE PRECISION FUNCTION dlamc3( A, B )
!
!  -- LAPACK auxiliary routine --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   a, b
!     ..
! =====================================================================
!
!     .. Executable Statements ..
!
      dlamc3 = a + b
!
      RETURN
!
!     End of DLAMC3
!

      END
!
!***********************************************************************

!> \brief \b DLAQR4 computes the eigenvalues of a Hessenberg matrix, and optionally the matrices from the Schur decomposition.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download DLAQR4 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaqr4.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaqr4.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaqr4.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAQR4( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI,
!                          ILOZ, IHIZ, Z, LDZ, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, LWORK, N
!       LOGICAL            WANTT, WANTZ
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   H( LDH, * ), WI( * ), WORK( * ), WR( * ),
!      $                   Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    DLAQR4 implements one level of recursion for DLAQR0.
!>    It is a complete implementation of the small bulge multi-shift
!>    QR algorithm.  It may be called by DLAQR0 and, for large enough
!>    deflation window size, it may be called by DLAQR3.  This
!>    subroutine is identical to DLAQR0 except that it calls DLAQR2
!>    instead of DLAQR3.
!>
!>    DLAQR4 computes the eigenvalues of a Hessenberg matrix H
!>    and, optionally, the matrices T and Z from the Schur decomposition
!>    H = Z T Z**T, where T is an upper quasi-triangular matrix (the
!>    Schur form), and Z is the orthogonal matrix of Schur vectors.
!>
!>    Optionally Z may be postmultiplied into an input orthogonal
!>    matrix Q so that this routine can give the Schur factorization
!>    of a matrix A which has been reduced to the Hessenberg form H
!>    by the orthogonal matrix Q:  A = Q*H*Q**T = (QZ)*T*(QZ)**T.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] WANTT
!> \verbatim
!>          WANTT is LOGICAL
!>          = .TRUE. : the full Schur form T is required;
!>          = .FALSE.: only eigenvalues are required.
!> \endverbatim
!>
!> \param[in] WANTZ
!> \verbatim
!>          WANTZ is LOGICAL
!>          = .TRUE. : the matrix of Schur vectors Z is required;
!>          = .FALSE.: Schur vectors are not required.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           The order of the matrix H.  N >= 0.
!> \endverbatim
!>
!> \param[in] ILO
!> \verbatim
!>          ILO is INTEGER
!> \endverbatim
!>
!> \param[in] IHI
!> \verbatim
!>          IHI is INTEGER
!>           It is assumed that H is already upper triangular in rows
!>           and columns 1:ILO-1 and IHI+1:N and, if ILO > 1,
!>           H(ILO,ILO-1) is zero. ILO and IHI are normally set by a
!>           previous call to DGEBAL, and then passed to DGEHRD when the
!>           matrix output by DGEBAL is reduced to Hessenberg form.
!>           Otherwise, ILO and IHI should be set to 1 and N,
!>           respectively.  If N > 0, then 1 <= ILO <= IHI <= N.
!>           If N = 0, then ILO = 1 and IHI = 0.
!> \endverbatim
!>
!> \param[in,out] H
!> \verbatim
!>          H is DOUBLE PRECISION array, dimension (LDH,N)
!>           On entry, the upper Hessenberg matrix H.
!>           On exit, if INFO = 0 and WANTT is .TRUE., then H contains
!>           the upper quasi-triangular matrix T from the Schur
!>           decomposition (the Schur form); 2-by-2 diagonal blocks
!>           (corresponding to complex conjugate pairs of eigenvalues)
!>           are returned in standard form, with H(i,i) = H(i+1,i+1)
!>           and H(i+1,i)*H(i,i+1) < 0. If INFO = 0 and WANTT is
!>           .FALSE., then the contents of H are unspecified on exit.
!>           (The output value of H when INFO > 0 is given under the
!>           description of INFO below.)
!>
!>           This subroutine may explicitly set H(i,j) = 0 for i > j and
!>           j = 1, 2, ... ILO-1 or j = IHI+1, IHI+2, ... N.
!> \endverbatim
!>
!> \param[in] LDH
!> \verbatim
!>          LDH is INTEGER
!>           The leading dimension of the array H. LDH >= max(1,N).
!> \endverbatim
!>
!> \param[out] WR
!> \verbatim
!>          WR is DOUBLE PRECISION array, dimension (IHI)
!> \endverbatim
!>
!> \param[out] WI
!> \verbatim
!>          WI is DOUBLE PRECISION array, dimension (IHI)
!>           The real and imaginary parts, respectively, of the computed
!>           eigenvalues of H(ILO:IHI,ILO:IHI) are stored in WR(ILO:IHI)
!>           and WI(ILO:IHI). If two eigenvalues are computed as a
!>           complex conjugate pair, they are stored in consecutive
!>           elements of WR and WI, say the i-th and (i+1)th, with
!>           WI(i) > 0 and WI(i+1) < 0. If WANTT is .TRUE., then
!>           the eigenvalues are stored in the same order as on the
!>           diagonal of the Schur form returned in H, with
!>           WR(i) = H(i,i) and, if H(i:i+1,i:i+1) is a 2-by-2 diagonal
!>           block, WI(i) = sqrt(-H(i+1,i)*H(i,i+1)) and
!>           WI(i+1) = -WI(i).
!> \endverbatim
!>
!> \param[in] ILOZ
!> \verbatim
!>          ILOZ is INTEGER
!> \endverbatim
!>
!> \param[in] IHIZ
!> \verbatim
!>          IHIZ is INTEGER
!>           Specify the rows of Z to which transformations must be
!>           applied if WANTZ is .TRUE..
!>           1 <= ILOZ <= ILO; IHI <= IHIZ <= N.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is DOUBLE PRECISION array, dimension (LDZ,IHI)
!>           If WANTZ is .FALSE., then Z is not referenced.
!>           If WANTZ is .TRUE., then Z(ILO:IHI,ILOZ:IHIZ) is
!>           replaced by Z(ILO:IHI,ILOZ:IHIZ)*U where U is the
!>           orthogonal Schur factor of H(ILO:IHI,ILO:IHI).
!>           (The output value of Z when INFO > 0 is given under
!>           the description of INFO below.)
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>           The leading dimension of the array Z.  if WANTZ is .TRUE.
!>           then LDZ >= MAX(1,IHIZ).  Otherwise, LDZ >= 1.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension LWORK
!>           On exit, if LWORK = -1, WORK(1) returns an estimate of
!>           the optimal value for LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>           The dimension of the array WORK.  LWORK >= max(1,N)
!>           is sufficient, but LWORK typically as large as 6*N may
!>           be required for optimal performance.  A workspace query
!>           to determine the optimal workspace size is recommended.
!>
!>           If LWORK = -1, then DLAQR4 does a workspace query.
!>           In this case, DLAQR4 checks the input parameters and
!>           estimates the optimal workspace size for the given
!>           values of N, ILO and IHI.  The estimate is returned
!>           in WORK(1).  No error message related to LWORK is
!>           issued by XERBLA.  Neither H nor Z are accessed.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>             = 0:  successful exit
!>             > 0:  if INFO = i, DLAQR4 failed to compute all of
!>                the eigenvalues.  Elements 1:ilo-1 and i+1:n of WR
!>                and WI contain those eigenvalues which have been
!>                successfully computed.  (Failures are rare.)
!>
!>                If INFO > 0 and WANT is .FALSE., then on exit,
!>                the remaining unconverged eigenvalues are the eigen-
!>                values of the upper Hessenberg matrix rows and
!>                columns ILO through INFO of the final, output
!>                value of H.
!>
!>                If INFO > 0 and WANTT is .TRUE., then on exit
!>
!>           (*)  (initial value of H)*U  = U*(final value of H)
!>
!>                where U is a orthogonal matrix.  The final
!>                value of  H is upper Hessenberg and triangular in
!>                rows and columns INFO+1 through IHI.
!>
!>                If INFO > 0 and WANTZ is .TRUE., then on exit
!>
!>                  (final value of Z(ILO:IHI,ILOZ:IHIZ)
!>                   =  (initial value of Z(ILO:IHI,ILOZ:IHIZ)*U
!>
!>                where U is the orthogonal matrix in (*) (regard-
!>                less of the value of WANTT.)
!>
!>                If INFO > 0 and WANTZ is .FALSE., then Z is not
!>                accessed.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup laqr4
!
!> \par Contributors:
!  ==================
!>
!>       Karen Braman and Ralph Byers, Department of Mathematics,
!>       University of Kansas, USA
!
!> \par References:
!  ================
!>
!>       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
!>       Algorithm Part I: Maintaining Well Focused Shifts, and Level 3
!>       Performance, SIAM Journal of Matrix Analysis, volume 23, pages
!>       929--947, 2002.
!> \n
!>       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
!>       Algorithm Part II: Aggressive Early Deflation, SIAM Journal
!>       of Matrix Analysis, volume 23, pages 948--973, 2002.
!>
!  =====================================================================

      SUBROUTINE dlaqr4( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, &
     &                   ILOZ, IHIZ, Z, LDZ, WORK, LWORK, INFO )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, LWORK, N
      LOGICAL            WANTT, WANTZ
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   H( LDH, * ), WI( * ), WORK( * ), WR( * ), &
     &                   z( ldz, * )
!     ..
!
!  ================================================================
!     .. Parameters ..
!
!     ==== Matrices of order NTINY or smaller must be processed by
!     .    DLAHQR because of insufficient subdiagonal scratch space.
!     .    (This is a hard limit.) ====
      INTEGER            NTINY
      parameter( ntiny = 15 )
!
!     ==== Exceptional deflation windows:  try to cure rare
!     .    slow convergence by varying the size of the
!     .    deflation window after KEXNW iterations. ====
      INTEGER            KEXNW
      parameter( kexnw = 5 )
!
!     ==== Exceptional shifts: try to cure rare slow convergence
!     .    with ad-hoc exceptional shifts every KEXSH iterations.
!     .    ====
      INTEGER            KEXSH
      parameter( kexsh = 6 )
!
!     ==== The constants WILK1 and WILK2 are used to form the
!     .    exceptional shifts. ====
      DOUBLE PRECISION   WILK1, WILK2
      parameter( wilk1 = 0.75d0, wilk2 = -0.4375d0 )
      DOUBLE PRECISION   ZERO, ONE
      parameter( zero = 0.0d0, one = 1.0d0 )
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION   AA, BB, CC, CS, DD, SN, SS, SWAP
      INTEGER            I, INF, IT, ITMAX, K, KACC22, KBOT, KDU, KS, &
     &                   kt, ktop, ku, kv, kwh, kwtop, kwv, ld, ls, &
     &                   lwkopt, ndec, ndfl, nh, nho, nibble, nmin, ns, &
     &                   nsmax, nsr, nve, nw, nwmax, nwr, nwupbd
      LOGICAL            SORTED
      CHARACTER          JBCMPZ*2
!     ..
!     .. External Functions ..
!      INTEGER            ILAENV
!      EXTERNAL           ilaenv
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION   ZDUM( 1, 1 )
!     ..
!     .. External Subroutines ..
!      EXTERNAL           dlacpy, dlahqr, dlanv2, dlaqr2, &
!     &                   dlaqr5
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          abs, dble, int, max, min, mod
!     ..
!     .. Executable Statements ..
      info = 0
!
!     ==== Quick return for N = 0: nothing to do. ====
!
      IF( n.EQ.0 ) THEN
         work( 1 ) = one
         RETURN
      END IF
!
      IF( n.LE.ntiny ) THEN
!
!        ==== Tiny matrices must use DLAHQR. ====
!
         lwkopt = 1
         IF( lwork.NE.-1 ) &
     &      CALL dlahqr( wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, &
     &                   iloz, ihiz, z, ldz, info )
      ELSE
!
!        ==== Use small bulge multi-shift QR with aggressive early
!        .    deflation on larger-than-tiny matrices. ====
!
!        ==== Hope for the best. ====
!
         info = 0
!
!        ==== Set up job flags for ILAENV. ====
!
         IF( wantt ) THEN
            jbcmpz( 1: 1 ) = 'S'
         ELSE
            jbcmpz( 1: 1 ) = 'E'
         END IF
         IF( wantz ) THEN
            jbcmpz( 2: 2 ) = 'V'
         ELSE
            jbcmpz( 2: 2 ) = 'N'
         END IF
!
!        ==== NWR = recommended deflation window size.  At this
!        .    point,  N .GT. NTINY = 15, so there is enough
!        .    subdiagonal workspace for NWR.GE.2 as required.
!        .    (In fact, there is enough subdiagonal space for
!        .    NWR.GE.4.) ====
!
         nwr = ilaenv( 13, 'DLAQR4', jbcmpz, n, ilo, ihi, lwork )
         nwr = max( 2, nwr )
         nwr = min( ihi-ilo+1, ( n-1 ) / 3, nwr )
!
!        ==== NSR = recommended number of simultaneous shifts.
!        .    At this point N .GT. NTINY = 15, so there is at
!        .    enough subdiagonal workspace for NSR to be even
!        .    and greater than or equal to two as required. ====
!
         nsr = ilaenv( 15, 'DLAQR4', jbcmpz, n, ilo, ihi, lwork )
         nsr = min( nsr, ( n-3 ) / 6, ihi-ilo )
         nsr = max( 2, nsr-mod( nsr, 2 ) )
!
!        ==== Estimate optimal workspace ====
!
!        ==== Workspace query call to DLAQR2 ====
!
         CALL dlaqr2( wantt, wantz, n, ilo, ihi, nwr+1, h, ldh, iloz, &
     &                ihiz, z, ldz, ls, ld, wr, wi, h, ldh, n, h, ldh, &
     &                n, h, ldh, work, -1 )
!
!        ==== Optimal workspace = MAX(DLAQR5, DLAQR2) ====
!
         lwkopt = max( 3*nsr / 2, int( work( 1 ) ) )
!
!        ==== Quick return in case of workspace query. ====
!
         IF( lwork.EQ.-1 ) THEN
            work( 1 ) = dble( lwkopt )
            RETURN
         END IF
!
!        ==== DLAHQR/DLAQR0 crossover point ====
!
         nmin = ilaenv( 12, 'DLAQR4', jbcmpz, n, ilo, ihi, lwork )
         nmin = max( ntiny, nmin )
!
!        ==== Nibble crossover point ====
!
         nibble = ilaenv( 14, 'DLAQR4', jbcmpz, n, ilo, ihi, lwork )
         nibble = max( 0, nibble )
!
!        ==== Accumulate reflections during ttswp?  Use block
!        .    2-by-2 structure during matrix-matrix multiply? ====
!
         kacc22 = ilaenv( 16, 'DLAQR4', jbcmpz, n, ilo, ihi, lwork )
         kacc22 = max( 0, kacc22 )
         kacc22 = min( 2, kacc22 )
!
!        ==== NWMAX = the largest possible deflation window for
!        .    which there is sufficient workspace. ====
!
         nwmax = min( ( n-1 ) / 3, lwork / 2 )
         nw = nwmax
!
!        ==== NSMAX = the Largest number of simultaneous shifts
!        .    for which there is sufficient workspace. ====
!
         nsmax = min( ( n-3 ) / 6, 2*lwork / 3 )
         nsmax = nsmax - mod( nsmax, 2 )
!
!        ==== NDFL: an iteration count restarted at deflation. ====
!
         ndfl = 1
!
!        ==== ITMAX = iteration limit ====
!
         itmax = max( 30, 2*kexsh )*max( 10, ( ihi-ilo+1 ) )
!
!        ==== Last row and column in the active block ====
!
         kbot = ihi
!
!        ==== Main Loop ====
!
         DO 80 it = 1, itmax
!
!           ==== Done when KBOT falls below ILO ====
!
            IF( kbot.LT.ilo ) &
     &         GO TO 90
!
!           ==== Locate active block ====
!
            DO 10 k = kbot, ilo + 1, -1
               IF( h( k, k-1 ).EQ.zero ) &
     &            GO TO 20
10          CONTINUE
            k = ilo
20          CONTINUE
            ktop = k
!
!           ==== Select deflation window size:
!           .    Typical Case:
!           .      If possible and advisable, nibble the entire
!           .      active block.  If not, use size MIN(NWR,NWMAX)
!           .      or MIN(NWR+1,NWMAX) depending upon which has
!           .      the smaller corresponding subdiagonal entry
!           .      (a heuristic).
!           .
!           .    Exceptional Case:
!           .      If there have been no deflations in KEXNW or
!           .      more iterations, then vary the deflation window
!           .      size.   At first, because, larger windows are,
!           .      in general, more powerful than smaller ones,
!           .      rapidly increase the window to the maximum possible.
!           .      Then, gradually reduce the window size. ====
!
            nh = kbot - ktop + 1
            nwupbd = min( nh, nwmax )
            IF( ndfl.LT.kexnw ) THEN
               nw = min( nwupbd, nwr )
            ELSE
               nw = min( nwupbd, 2*nw )
            END IF
            IF( nw.LT.nwmax ) THEN
               IF( nw.GE.nh-1 ) THEN
                  nw = nh
               ELSE
                  kwtop = kbot - nw + 1
                  IF( abs( h( kwtop, kwtop-1 ) ).GT. &
     &                abs( h( kwtop-1, kwtop-2 ) ) )nw = nw + 1
               END IF
            END IF
            IF( ndfl.LT.kexnw ) THEN
               ndec = -1
            ELSE IF( ndec.GE.0 .OR. nw.GE.nwupbd ) THEN
               ndec = ndec + 1
               IF( nw-ndec.LT.2 ) &
     &            ndec = 0
               nw = nw - ndec
            END IF
!
!           ==== Aggressive early deflation:
!           .    split workspace under the subdiagonal into
!           .      - an nw-by-nw work array V in the lower
!           .        left-hand-corner,
!           .      - an NW-by-at-least-NW-but-more-is-better
!           .        (NW-by-NHO) horizontal work array along
!           .        the bottom edge,
!           .      - an at-least-NW-but-more-is-better (NHV-by-NW)
!           .        vertical work array along the left-hand-edge.
!           .        ====
!
            kv = n - nw + 1
            kt = nw + 1
            nho = ( n-nw-1 ) - kt + 1
            kwv = nw + 2
            nve = ( n-nw ) - kwv + 1
!
!           ==== Aggressive early deflation ====
!
            CALL dlaqr2( wantt, wantz, n, ktop, kbot, nw, h, ldh, &
     &                   iloz, &
     &                   ihiz, z, ldz, ls, ld, wr, wi, h( kv, 1 ), ldh, &
     &                   nho, h( kv, kt ), ldh, nve, h( kwv, 1 ), ldh, &
     &                   work, lwork )
!
!           ==== Adjust KBOT accounting for new deflations. ====
!
            kbot = kbot - ld
!
!           ==== KS points to the shifts. ====
!
            ks = kbot - ls + 1
!
!           ==== Skip an expensive QR sweep if there is a (partly
!           .    heuristic) reason to expect that many eigenvalues
!           .    will deflate without it.  Here, the QR sweep is
!           .    skipped if many eigenvalues have just been deflated
!           .    or if the remaining active block is small.
!
            IF( ( ld.EQ.0 ) .OR. ( ( 100*ld.LE.nw*nibble ) .AND. ( kbot- &
     &          ktop+1.GT.min( nmin, nwmax ) ) ) ) THEN
!
!              ==== NS = nominal number of simultaneous shifts.
!              .    This may be lowered (slightly) if DLAQR2
!              .    did not provide that many shifts. ====
!
               ns = min( nsmax, nsr, max( 2, kbot-ktop ) )
               ns = ns - mod( ns, 2 )
!
!              ==== If there have been no deflations
!              .    in a multiple of KEXSH iterations,
!              .    then try exceptional shifts.
!              .    Otherwise use shifts provided by
!              .    DLAQR2 above or from the eigenvalues
!              .    of a trailing principal submatrix. ====
!
               IF( mod( ndfl, kexsh ).EQ.0 ) THEN
                  ks = kbot - ns + 1
                  DO 30 i = kbot, max( ks+1, ktop+2 ), -2
                     ss = abs( h( i, i-1 ) ) + abs( h( i-1, i-2 ) )
                     aa = wilk1*ss + h( i, i )
                     bb = ss
                     cc = wilk2*ss
                     dd = aa
                     CALL dlanv2( aa, bb, cc, dd, wr( i-1 ), &
     &                            wi( i-1 ), &
     &                            wr( i ), wi( i ), cs, sn )
30                CONTINUE
                  IF( ks.EQ.ktop ) THEN
                     wr( ks+1 ) = h( ks+1, ks+1 )
                     wi( ks+1 ) = zero
                     wr( ks ) = wr( ks+1 )
                     wi( ks ) = wi( ks+1 )
                  END IF
               ELSE
!
!                 ==== Got NS/2 or fewer shifts? Use DLAHQR
!                 .    on a trailing principal submatrix to
!                 .    get more. (Since NS.LE.NSMAX.LE.(N-3)/6,
!                 .    there is enough space below the subdiagonal
!                 .    to fit an NS-by-NS scratch array.) ====
!
                  IF( kbot-ks+1.LE.ns / 2 ) THEN
                     ks = kbot - ns + 1
                     kt = n - ns + 1
                     CALL dlacpy( 'A', ns, ns, h( ks, ks ), ldh, &
     &                            h( kt, 1 ), ldh )
                     CALL dlahqr( .false., .false., ns, 1, ns, &
     &                            h( kt, 1 ), ldh, wr( ks ), wi( ks ), &
     &                            1, 1, zdum, 1, inf )
                     ks = ks + inf
!
!                    ==== In case of a rare QR failure use
!                    .    eigenvalues of the trailing 2-by-2
!                    .    principal submatrix.  ====
!
                     IF( ks.GE.kbot ) THEN
                        aa = h( kbot-1, kbot-1 )
                        cc = h( kbot, kbot-1 )
                        bb = h( kbot-1, kbot )
                        dd = h( kbot, kbot )
                        CALL dlanv2( aa, bb, cc, dd, wr( kbot-1 ), &
     &                               wi( kbot-1 ), wr( kbot ), &
     &                               wi( kbot ), cs, sn )
                        ks = kbot - 1
                     END IF
                  END IF
!
                  IF( kbot-ks+1.GT.ns ) THEN
!
!                    ==== Sort the shifts (Helps a little)
!                    .    Bubble sort keeps complex conjugate
!                    .    pairs together. ====
!
                     sorted = .false.
                     DO 50 k = kbot, ks + 1, -1
                        IF( sorted ) &
     &                     GO TO 60
                        sorted = .true.
                        DO 40 i = ks, k - 1
                           IF( abs( wr( i ) )+abs( wi( i ) ).LT. &
     &                         abs( wr( i+1 ) )+abs( wi( i+1 ) ) ) THEN
                              sorted = .false.
!
                              swap = wr( i )
                              wr( i ) = wr( i+1 )
                              wr( i+1 ) = swap
!
                              swap = wi( i )
                              wi( i ) = wi( i+1 )
                              wi( i+1 ) = swap
                           END IF
40                      CONTINUE
50                   CONTINUE
60                   CONTINUE
                  END IF
!
!                 ==== Shuffle shifts into pairs of real shifts
!                 .    and pairs of complex conjugate shifts
!                 .    assuming complex conjugate shifts are
!                 .    already adjacent to one another. (Yes,
!                 .    they are.)  ====
!
                  DO 70 i = kbot, ks + 2, -2
                     IF( wi( i ).NE.-wi( i-1 ) ) THEN
!
                        swap = wr( i )
                        wr( i ) = wr( i-1 )
                        wr( i-1 ) = wr( i-2 )
                        wr( i-2 ) = swap
!
                        swap = wi( i )
                        wi( i ) = wi( i-1 )
                        wi( i-1 ) = wi( i-2 )
                        wi( i-2 ) = swap
                     END IF
70                CONTINUE
               END IF
!
!              ==== If there are only two shifts and both are
!              .    real, then use only one.  ====
!
               IF( kbot-ks+1.EQ.2 ) THEN
                  IF( wi( kbot ).EQ.zero ) THEN
                     IF( abs( wr( kbot )-h( kbot, kbot ) ).LT. &
     &                   abs( wr( kbot-1 )-h( kbot, kbot ) ) ) THEN
                        wr( kbot-1 ) = wr( kbot )
                     ELSE
                        wr( kbot ) = wr( kbot-1 )
                     END IF
                  END IF
               END IF
!
!              ==== Use up to NS of the the smallest magnitude
!              .    shifts.  If there aren't NS shifts available,
!              .    then use them all, possibly dropping one to
!              .    make the number of shifts even. ====
!
               ns = min( ns, kbot-ks+1 )
               ns = ns - mod( ns, 2 )
               ks = kbot - ns + 1
!
!              ==== Small-bulge multi-shift QR sweep:
!              .    split workspace under the subdiagonal into
!              .    - a KDU-by-KDU work array U in the lower
!              .      left-hand-corner,
!              .    - a KDU-by-at-least-KDU-but-more-is-better
!              .      (KDU-by-NHo) horizontal work array WH along
!              .      the bottom edge,
!              .    - and an at-least-KDU-but-more-is-better-by-KDU
!              .      (NVE-by-KDU) vertical work WV arrow along
!              .      the left-hand-edge. ====
!
               kdu = 2*ns
               ku = n - kdu + 1
               kwh = kdu + 1
               nho = ( n-kdu+1-4 ) - ( kdu+1 ) + 1
               kwv = kdu + 4
               nve = n - kdu - kwv + 1
!
!              ==== Small-bulge multi-shift QR sweep ====
!
               CALL dlaqr5( wantt, wantz, kacc22, n, ktop, kbot, ns, &
     &                      wr( ks ), wi( ks ), h, ldh, iloz, ihiz, z, &
     &                      ldz, work, 3, h( ku, 1 ), ldh, nve, &
     &                      h( kwv, 1 ), ldh, nho, h( ku, kwh ), ldh )
            END IF
!
!           ==== Note progress (or the lack of it). ====
!
            IF( ld.GT.0 ) THEN
               ndfl = 1
            ELSE
               ndfl = ndfl + 1
            END IF
!
!           ==== End of main loop ====
80       CONTINUE
!
!        ==== Iteration limit exceeded.  Set INFO to show where
!        .    the problem occurred and exit. ====
!
         info = kbot
90       CONTINUE
      END IF
!
!     ==== Return the optimal value of LWORK. ====
!
      work( 1 ) = dble( lwkopt )
!
!     ==== End of DLAQR4 ====
!

      END
!> \brief \b DLAQR3 performs the orthogonal similarity transformation of a Hessenberg matrix to detect and deflate fully converged eigenvalues from a trailing principal submatrix (aggressive early deflation).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download DLAQR3 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaqr3.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaqr3.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaqr3.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAQR3( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ,
!                          IHIZ, Z, LDZ, NS, ND, SR, SI, V, LDV, NH, T,
!                          LDT, NV, WV, LDWV, WORK, LWORK )
!
!       .. Scalar Arguments ..
!       INTEGER            IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV,
!      $                   LDZ, LWORK, N, ND, NH, NS, NV, NW
!       LOGICAL            WANTT, WANTZ
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   H( LDH, * ), SI( * ), SR( * ), T( LDT, * ),
!      $                   V( LDV, * ), WORK( * ), WV( LDWV, * ),
!      $                   Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    Aggressive early deflation:
!>
!>    DLAQR3 accepts as input an upper Hessenberg matrix
!>    H and performs an orthogonal similarity transformation
!>    designed to detect and deflate fully converged eigenvalues from
!>    a trailing principal submatrix.  On output H has been over-
!>    written by a new Hessenberg matrix that is a perturbation of
!>    an orthogonal similarity transformation of H.  It is to be
!>    hoped that the final version of H has many zero subdiagonal
!>    entries.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] WANTT
!> \verbatim
!>          WANTT is LOGICAL
!>          If .TRUE., then the Hessenberg matrix H is fully updated
!>          so that the quasi-triangular Schur factor may be
!>          computed (in cooperation with the calling subroutine).
!>          If .FALSE., then only enough of H is updated to preserve
!>          the eigenvalues.
!> \endverbatim
!>
!> \param[in] WANTZ
!> \verbatim
!>          WANTZ is LOGICAL
!>          If .TRUE., then the orthogonal matrix Z is updated so
!>          so that the orthogonal Schur factor may be computed
!>          (in cooperation with the calling subroutine).
!>          If .FALSE., then Z is not referenced.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix H and (if WANTZ is .TRUE.) the
!>          order of the orthogonal matrix Z.
!> \endverbatim
!>
!> \param[in] KTOP
!> \verbatim
!>          KTOP is INTEGER
!>          It is assumed that either KTOP = 1 or H(KTOP,KTOP-1)=0.
!>          KBOT and KTOP together determine an isolated block
!>          along the diagonal of the Hessenberg matrix.
!> \endverbatim
!>
!> \param[in] KBOT
!> \verbatim
!>          KBOT is INTEGER
!>          It is assumed without a check that either
!>          KBOT = N or H(KBOT+1,KBOT)=0.  KBOT and KTOP together
!>          determine an isolated block along the diagonal of the
!>          Hessenberg matrix.
!> \endverbatim
!>
!> \param[in] NW
!> \verbatim
!>          NW is INTEGER
!>          Deflation window size.  1 <= NW <= (KBOT-KTOP+1).
!> \endverbatim
!>
!> \param[in,out] H
!> \verbatim
!>          H is DOUBLE PRECISION array, dimension (LDH,N)
!>          On input the initial N-by-N section of H stores the
!>          Hessenberg matrix undergoing aggressive early deflation.
!>          On output H has been transformed by an orthogonal
!>          similarity transformation, perturbed, and the returned
!>          to Hessenberg form that (it is to be hoped) has some
!>          zero subdiagonal entries.
!> \endverbatim
!>
!> \param[in] LDH
!> \verbatim
!>          LDH is INTEGER
!>          Leading dimension of H just as declared in the calling
!>          subroutine.  N <= LDH
!> \endverbatim
!>
!> \param[in] ILOZ
!> \verbatim
!>          ILOZ is INTEGER
!> \endverbatim
!>
!> \param[in] IHIZ
!> \verbatim
!>          IHIZ is INTEGER
!>          Specify the rows of Z to which transformations must be
!>          applied if WANTZ is .TRUE.. 1 <= ILOZ <= IHIZ <= N.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is DOUBLE PRECISION array, dimension (LDZ,N)
!>          IF WANTZ is .TRUE., then on output, the orthogonal
!>          similarity transformation mentioned above has been
!>          accumulated into Z(ILOZ:IHIZ,ILOZ:IHIZ) from the right.
!>          If WANTZ is .FALSE., then Z is unreferenced.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of Z just as declared in the
!>          calling subroutine.  1 <= LDZ.
!> \endverbatim
!>
!> \param[out] NS
!> \verbatim
!>          NS is INTEGER
!>          The number of unconverged (ie approximate) eigenvalues
!>          returned in SR and SI that may be used as shifts by the
!>          calling subroutine.
!> \endverbatim
!>
!> \param[out] ND
!> \verbatim
!>          ND is INTEGER
!>          The number of converged eigenvalues uncovered by this
!>          subroutine.
!> \endverbatim
!>
!> \param[out] SR
!> \verbatim
!>          SR is DOUBLE PRECISION array, dimension (KBOT)
!> \endverbatim
!>
!> \param[out] SI
!> \verbatim
!>          SI is DOUBLE PRECISION array, dimension (KBOT)
!>          On output, the real and imaginary parts of approximate
!>          eigenvalues that may be used for shifts are stored in
!>          SR(KBOT-ND-NS+1) through SR(KBOT-ND) and
!>          SI(KBOT-ND-NS+1) through SI(KBOT-ND), respectively.
!>          The real and imaginary parts of converged eigenvalues
!>          are stored in SR(KBOT-ND+1) through SR(KBOT) and
!>          SI(KBOT-ND+1) through SI(KBOT), respectively.
!> \endverbatim
!>
!> \param[out] V
!> \verbatim
!>          V is DOUBLE PRECISION array, dimension (LDV,NW)
!>          An NW-by-NW work array.
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is INTEGER
!>          The leading dimension of V just as declared in the
!>          calling subroutine.  NW <= LDV
!> \endverbatim
!>
!> \param[in] NH
!> \verbatim
!>          NH is INTEGER
!>          The number of columns of T.  NH >= NW.
!> \endverbatim
!>
!> \param[out] T
!> \verbatim
!>          T is DOUBLE PRECISION array, dimension (LDT,NW)
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of T just as declared in the
!>          calling subroutine.  NW <= LDT
!> \endverbatim
!>
!> \param[in] NV
!> \verbatim
!>          NV is INTEGER
!>          The number of rows of work array WV available for
!>          workspace.  NV >= NW.
!> \endverbatim
!>
!> \param[out] WV
!> \verbatim
!>          WV is DOUBLE PRECISION array, dimension (LDWV,NW)
!> \endverbatim
!>
!> \param[in] LDWV
!> \verbatim
!>          LDWV is INTEGER
!>          The leading dimension of W just as declared in the
!>          calling subroutine.  NW <= LDV
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (LWORK)
!>          On exit, WORK(1) is set to an estimate of the optimal value
!>          of LWORK for the given values of N, NW, KTOP and KBOT.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the work array WORK.  LWORK = 2*NW
!>          suffices, but greater efficiency may result from larger
!>          values of LWORK.
!>
!>          If LWORK = -1, then a workspace query is assumed; DLAQR3
!>          only estimates the optimal workspace size for the given
!>          values of N, NW, KTOP and KBOT.  The estimate is returned
!>          in WORK(1).  No error message related to LWORK is issued
!>          by XERBLA.  Neither H nor Z are accessed.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup laqr3
!
!> \par Contributors:
!  ==================
!>
!>       Karen Braman and Ralph Byers, Department of Mathematics,
!>       University of Kansas, USA
!>
!  =====================================================================

      SUBROUTINE dlaqr3( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, &
     &                   ILOZ, &
     &                   IHIZ, Z, LDZ, NS, ND, SR, SI, V, LDV, NH, T, &
     &                   LDT, NV, WV, LDWV, WORK, LWORK )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV, &
     &                   LDZ, LWORK, N, ND, NH, NS, NV, NW
      LOGICAL            WANTT, WANTZ
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   H( LDH, * ), SI( * ), SR( * ), T( LDT, * ), &
     &                   V( LDV, * ), WORK( * ), WV( LDWV, * ), &
     &                   Z( LDZ, * )
!     ..
!
!  ================================================================
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0d0, one = 1.0d0 )
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION   AA, BB, BETA, CC, CS, DD, EVI, EVK, FOO, S, &
     &                   SAFMAX, SAFMIN, SMLNUM, SN, TAU, ULP
      INTEGER            I, IFST, ILST, INFO, INFQR, J, JW, K, KCOL, &
     &                   kend, kln, krow, kwtop, ltop, lwk1, lwk2, lwk3, &
     &                   lwkopt, nmin
      LOGICAL            BULGE, SORTED
!     ..
!     .. External Functions ..
!      DOUBLE PRECISION   DLAMCH
!      INTEGER            ILAENV
!      EXTERNAL           DLAMCH, ILAENV
!     ..
!     .. External Subroutines ..
!      EXTERNAL           dcopy, dgehrd, dgemm, dlacpy, dlahqr, &
!     &                   dlanv2, &
!     &                   dlaqr4, dlarf1f, dlarfg, dlaset, dormhr, dtrexc
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          abs, dble, int, max, min, sqrt
!     ..
!     .. Executable Statements ..
!
!     ==== Estimate optimal workspace. ====
!
      jw = min( nw, kbot-ktop+1 )
      IF( jw.LE.2 ) THEN
         lwkopt = 1
      ELSE
!
!        ==== Workspace query call to DGEHRD ====
!
         CALL dgehrd( jw, 1, jw-1, t, ldt, work, work, -1, info )
         lwk1 = int( work( 1 ) )
!
!        ==== Workspace query call to DORMHR ====
!
         CALL dormhr( 'R', 'N', jw, jw, 1, jw-1, t, ldt, work, v, &
     &                ldv, &
     &                work, -1, info )
         lwk2 = int( work( 1 ) )
!
!        ==== Workspace query call to DLAQR4 ====
!
         CALL dlaqr4( .true., .true., jw, 1, jw, t, ldt, sr, si, 1, &
     &                jw, &
     &                v, ldv, work, -1, infqr )
         lwk3 = int( work( 1 ) )
!
!        ==== Optimal workspace ====
!
         lwkopt = max( jw+max( lwk1, lwk2 ), lwk3 )
      END IF
!
!     ==== Quick return in case of workspace query. ====
!
      IF( lwork.EQ.-1 ) THEN
         work( 1 ) = dble( lwkopt )
         RETURN
      END IF
!
!     ==== Nothing to do ...
!     ... for an empty active block ... ====
      ns = 0
      nd = 0
      work( 1 ) = one
      IF( ktop.GT.kbot ) &
     &   RETURN
!     ... nor for an empty deflation window. ====
      IF( nw.LT.1 ) &
     &   RETURN
!
!     ==== Machine constants ====
!
      safmin = dlamch( 'SAFE MINIMUM' )
      safmax = one / safmin
      ulp = dlamch( 'PRECISION' )
      smlnum = safmin*( dble( n ) / ulp )
!
!     ==== Setup deflation window ====
!
      jw = min( nw, kbot-ktop+1 )
      kwtop = kbot - jw + 1
      IF( kwtop.EQ.ktop ) THEN
         s = zero
      ELSE
         s = h( kwtop, kwtop-1 )
      END IF
!
      IF( kbot.EQ.kwtop ) THEN
!
!        ==== 1-by-1 deflation window: not much to do ====
!
         sr( kwtop ) = h( kwtop, kwtop )
         si( kwtop ) = zero
         ns = 1
         nd = 0
         IF( abs( s ).LE.max( smlnum, ulp*abs( h( kwtop, kwtop ) ) ) ) &
     &        THEN
            ns = 0
            nd = 1
            IF( kwtop.GT.ktop ) &
     &         h( kwtop, kwtop-1 ) = zero
         END IF
         work( 1 ) = one
         RETURN
      END IF
!
!     ==== Convert to spike-triangular form.  (In case of a
!     .    rare QR failure, this routine continues to do
!     .    aggressive early deflation using that part of
!     .    the deflation window that converged using INFQR
!     .    here and there to keep track.) ====
!
      CALL dlacpy( 'U', jw, jw, h( kwtop, kwtop ), ldh, t, ldt )
      CALL dcopy( jw-1, h( kwtop+1, kwtop ), ldh+1, t( 2, 1 ), &
     &            ldt+1 )
!
      CALL dlaset( 'A', jw, jw, zero, one, v, ldv )
      nmin = ilaenv( 12, 'DLAQR3', 'SV', jw, 1, jw, lwork )
      IF( jw.GT.nmin ) THEN
         CALL dlaqr4( .true., .true., jw, 1, jw, t, ldt, sr( kwtop ), &
     &                si( kwtop ), 1, jw, v, ldv, work, lwork, infqr )
      ELSE
         CALL dlahqr( .true., .true., jw, 1, jw, t, ldt, sr( kwtop ), &
     &                si( kwtop ), 1, jw, v, ldv, infqr )
      END IF
!
!     ==== DTREXC needs a clean margin near the diagonal ====
!
      DO 10 j = 1, jw - 3
         t( j+2, j ) = zero
         t( j+3, j ) = zero
10    CONTINUE
      IF( jw.GT.2 ) &
     &   t( jw, jw-2 ) = zero
!
!     ==== Deflation detection loop ====
!
      ns = jw
      ilst = infqr + 1
20    CONTINUE
      IF( ilst.LE.ns ) THEN
         IF( ns.EQ.1 ) THEN
            bulge = .false.
         ELSE
            bulge = t( ns, ns-1 ).NE.zero
         END IF
!
!        ==== Small spike tip test for deflation ====
!
         IF( .NOT. bulge ) THEN
!
!           ==== Real eigenvalue ====
!
            foo = abs( t( ns, ns ) )
            IF( foo.EQ.zero ) &
     &         foo = abs( s )
            IF( abs( s*v( 1, ns ) ).LE.max( smlnum, ulp*foo ) ) THEN
!
!              ==== Deflatable ====
!
               ns = ns - 1
            ELSE
!
!              ==== Undeflatable.   Move it up out of the way.
!              .    (DTREXC can not fail in this case.) ====
!
               ifst = ns
               CALL dtrexc( 'V', jw, t, ldt, v, ldv, ifst, ilst, &
     &                      work, &
     &                      info )
               ilst = ilst + 1
            END IF
         ELSE
!
!           ==== Complex conjugate pair ====
!
            foo = abs( t( ns, ns ) ) + sqrt( abs( t( ns, ns-1 ) ) )* &
     &            sqrt( abs( t( ns-1, ns ) ) )
            IF( foo.EQ.zero ) &
     &         foo = abs( s )
            IF( max( abs( s*v( 1, ns ) ), abs( s*v( 1, ns-1 ) ) ).LE. &
     &          max( smlnum, ulp*foo ) ) THEN
!
!              ==== Deflatable ====
!
               ns = ns - 2
            ELSE
!
!              ==== Undeflatable. Move them up out of the way.
!              .    Fortunately, DTREXC does the right thing with
!              .    ILST in case of a rare exchange failure. ====
!
               ifst = ns
               CALL dtrexc( 'V', jw, t, ldt, v, ldv, ifst, ilst, &
     &                      work, &
     &                      info )
               ilst = ilst + 2
            END IF
         END IF
!
!        ==== End deflation detection loop ====
!
         GO TO 20
      END IF
!
!        ==== Return to Hessenberg form ====
!
      IF( ns.EQ.0 ) &
     &   s = zero
!
      IF( ns.LT.jw ) THEN
!
!        ==== sorting diagonal blocks of T improves accuracy for
!        .    graded matrices.  Bubble sort deals well with
!        .    exchange failures. ====
!
         sorted = .false.
         i = ns + 1
30       CONTINUE
         IF( sorted ) &
     &      GO TO 50
         sorted = .true.
!
         kend = i - 1
         i = infqr + 1
         IF( i.EQ.ns ) THEN
            k = i + 1
         ELSE IF( t( i+1, i ).EQ.zero ) THEN
            k = i + 1
         ELSE
            k = i + 2
         END IF
40       CONTINUE
         IF( k.LE.kend ) THEN
            IF( k.EQ.i+1 ) THEN
               evi = abs( t( i, i ) )
            ELSE
               evi = abs( t( i, i ) ) + sqrt( abs( t( i+1, i ) ) )* &
     &               sqrt( abs( t( i, i+1 ) ) )
            END IF
!
            IF( k.EQ.kend ) THEN
               evk = abs( t( k, k ) )
            ELSE IF( t( k+1, k ).EQ.zero ) THEN
               evk = abs( t( k, k ) )
            ELSE
               evk = abs( t( k, k ) ) + sqrt( abs( t( k+1, k ) ) )* &
     &               sqrt( abs( t( k, k+1 ) ) )
            END IF
!
            IF( evi.GE.evk ) THEN
               i = k
            ELSE
               sorted = .false.
               ifst = i
               ilst = k
               CALL dtrexc( 'V', jw, t, ldt, v, ldv, ifst, ilst, &
     &                      work, &
     &                      info )
               IF( info.EQ.0 ) THEN
                  i = ilst
               ELSE
                  i = k
               END IF
            END IF
            IF( i.EQ.kend ) THEN
               k = i + 1
            ELSE IF( t( i+1, i ).EQ.zero ) THEN
               k = i + 1
            ELSE
               k = i + 2
            END IF
            GO TO 40
         END IF
         GO TO 30
50       CONTINUE
      END IF
!
!     ==== Restore shift/eigenvalue array from T ====
!
      i = jw
60    CONTINUE
      IF( i.GE.infqr+1 ) THEN
         IF( i.EQ.infqr+1 ) THEN
            sr( kwtop+i-1 ) = t( i, i )
            si( kwtop+i-1 ) = zero
            i = i - 1
         ELSE IF( t( i, i-1 ).EQ.zero ) THEN
            sr( kwtop+i-1 ) = t( i, i )
            si( kwtop+i-1 ) = zero
            i = i - 1
         ELSE
            aa = t( i-1, i-1 )
            cc = t( i, i-1 )
            bb = t( i-1, i )
            dd = t( i, i )
            CALL dlanv2( aa, bb, cc, dd, sr( kwtop+i-2 ), &
     &                   si( kwtop+i-2 ), sr( kwtop+i-1 ), &
     &                   si( kwtop+i-1 ), cs, sn )
            i = i - 2
         END IF
         GO TO 60
      END IF
!
      IF( ns.LT.jw .OR. s.EQ.zero ) THEN
         IF( ns.GT.1 .AND. s.NE.zero ) THEN
!
!           ==== Reflect spike back into lower triangle ====
!
            CALL dcopy( ns, v, ldv, work, 1 )
            beta = work( 1 )
            CALL dlarfg( ns, beta, work( 2 ), 1, tau )
!
            CALL dlaset( 'L', jw-2, jw-2, zero, zero, t( 3, 1 ), &
     &                   ldt )
!
            CALL dlarf1f( 'L', ns, jw, work, 1, tau, t, ldt, &
     &                  work( jw+1 ) )
            CALL dlarf1f( 'R', ns, ns, work, 1, tau, t, ldt, &
     &                  work( jw+1 ) )
            CALL dlarf1f( 'R', jw, ns, work, 1, tau, v, ldv, &
     &                  work( jw+1 ) )
!
            CALL dgehrd( jw, 1, ns, t, ldt, work, work( jw+1 ), &
     &                   lwork-jw, info )
         END IF
!
!        ==== Copy updated reduced window into place ====
!
         IF( kwtop.GT.1 ) &
     &      h( kwtop, kwtop-1 ) = s*v( 1, 1 )
         CALL dlacpy( 'U', jw, jw, t, ldt, h( kwtop, kwtop ), ldh )
         CALL dcopy( jw-1, t( 2, 1 ), ldt+1, h( kwtop+1, kwtop ), &
     &               ldh+1 )
!
!        ==== Accumulate orthogonal matrix in order update
!        .    H and Z, if requested.  ====
!
         IF( ns.GT.1 .AND. s.NE.zero ) &
     &      CALL dormhr( 'R', 'N', jw, ns, 1, ns, t, ldt, work, v, &
     &                   ldv, &
     &                   work( jw+1 ), lwork-jw, info )
!
!        ==== Update vertical slab in H ====
!
         IF( wantt ) THEN
            ltop = 1
         ELSE
            ltop = ktop
         END IF
         DO 70 krow = ltop, kwtop - 1, nv
            kln = min( nv, kwtop-krow )
            CALL dgemm( 'N', 'N', kln, jw, jw, one, h( krow, kwtop ), &
     &                  ldh, v, ldv, zero, wv, ldwv )
            CALL dlacpy( 'A', kln, jw, wv, ldwv, h( krow, kwtop ), &
     &                   ldh )
70       CONTINUE
!
!        ==== Update horizontal slab in H ====
!
         IF( wantt ) THEN
            DO 80 kcol = kbot + 1, n, nh
               kln = min( nh, n-kcol+1 )
               CALL dgemm( 'C', 'N', jw, kln, jw, one, v, ldv, &
     &                     h( kwtop, kcol ), ldh, zero, t, ldt )
               CALL dlacpy( 'A', jw, kln, t, ldt, h( kwtop, kcol ), &
     &                      ldh )
80          CONTINUE
         END IF
!
!        ==== Update vertical slab in Z ====
!
         IF( wantz ) THEN
            DO 90 krow = iloz, ihiz, nv
               kln = min( nv, ihiz-krow+1 )
               CALL dgemm( 'N', 'N', kln, jw, jw, one, z( krow, &
     &                     kwtop ), &
     &                     ldz, v, ldv, zero, wv, ldwv )
               CALL dlacpy( 'A', kln, jw, wv, ldwv, z( krow, kwtop ), &
     &                      ldz )
90          CONTINUE
         END IF
      END IF
!
!     ==== Return the number of deflations ... ====
!
      nd = jw - ns
!
!     ==== ... and the number of shifts. (Subtracting
!     .    INFQR from the spike length takes care
!     .    of the case of a rare QR failure while
!     .    calculating eigenvalues of the deflation
!     .    window.)  ====
!
      ns = ns - infqr
!
!      ==== Return optimal workspace. ====
!
      work( 1 ) = dble( lwkopt )
!
!     ==== End of DLAQR3 ====
!

      END
!> \brief \b DLARFT forms the triangular factor T of a block reflector H = I - vtvH
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download DLARFT + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarft.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarft.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarft.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       RECURSIVE SUBROUTINE DLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIRECT, STOREV
!       INTEGER            K, LDT, LDV, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   T( LDT, * ), TAU( * ), V( LDV, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLARFT forms the triangular factor T of a real block reflector H
!> of order n, which is defined as a product of k elementary reflectors.
!>
!> If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
!>
!> If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
!>
!> If STOREV = 'C', the vector which defines the elementary reflector
!> H(i) is stored in the i-th column of the array V, and
!>
!>    H  =  I - V * T * V**T
!>
!> If STOREV = 'R', the vector which defines the elementary reflector
!> H(i) is stored in the i-th row of the array V, and
!>
!>    H  =  I - V**T * T * V
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] DIRECT
!> \verbatim
!>          DIRECT is CHARACTER*1
!>          Specifies the order in which the elementary reflectors are
!>          multiplied to form the block reflector:
!>          = 'F': H = H(1) H(2) . . . H(k) (Forward)
!>          = 'B': H = H(k) . . . H(2) H(1) (Backward)
!> \endverbatim
!>
!> \param[in] STOREV
!> \verbatim
!>          STOREV is CHARACTER*1
!>          Specifies how the vectors which define the elementary
!>          reflectors are stored (see also Further Details):
!>          = 'C': columnwise
!>          = 'R': rowwise
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the block reflector H. N >= 0.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The order of the triangular factor T (= the number of
!>          elementary reflectors). K >= 1.
!> \endverbatim
!>
!> \param[in] V
!> \verbatim
!>          V is DOUBLE PRECISION array, dimension
!>                               (LDV,K) if STOREV = 'C'
!>                               (LDV,N) if STOREV = 'R'
!>          The matrix V. See further details.
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is INTEGER
!>          The leading dimension of the array V.
!>          If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i).
!> \endverbatim
!>
!> \param[out] T
!> \verbatim
!>          T is DOUBLE PRECISION array, dimension (LDT,K)
!>          The k by k triangular factor T of the block reflector.
!>          If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is
!>          lower triangular. The rest of the array is not used.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of the array T. LDT >= K.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup larft
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The shape of the matrix V and the storage of the vectors which define
!>  the H(i) is best illustrated by the following example with n = 5 and
!>  k = 3. The elements equal to 1 are not stored.
!>
!>  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':
!>
!>               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
!>                   ( v1  1    )                     (     1 v2 v2 v2 )
!>                   ( v1 v2  1 )                     (        1 v3 v3 )
!>                   ( v1 v2 v3 )
!>                   ( v1 v2 v3 )
!>
!>  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':
!>
!>               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
!>                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
!>                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
!>                   (     1 v3 )
!>                   (        1 )
!> \endverbatim
!>
!  =====================================================================

      RECURSIVE SUBROUTINE dlarft( DIRECT, STOREV, N, K, V, LDV, &
     &                             TAU, T, LDT )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!        .. Scalar Arguments
!
      CHARACTER          direct, storev
      INTEGER            k, ldt, ldv, n
!     ..
!     .. Array Arguments ..
!
      DOUBLE PRECISION   t( ldt, * ), tau( * ), v( ldv, * )
!     ..
!
!     .. Parameters ..
!
      DOUBLE PRECISION one, neg_one, zero
      parameter(one=1.0d+0, zero = 0.0d+0, neg_one=-1.0d+0)
!
!     .. Local Scalars ..
!
      INTEGER           i,j,l
      LOGICAL           qr,lq,ql,dirf,colv
!
!     .. External Subroutines ..
!
!      EXTERNAL          dtrmm,dgemm,dlacpy
!
!     .. External Functions..
!
!      LOGICAL           lsame
!      EXTERNAL          lsame
!     
!     The general scheme used is inspired by the approach inside DGEQRT3
!     which was (at the time of writing this code):
!     Based on the algorithm of Elmroth and Gustavson,
!     IBM J. Res. Develop. Vol 44 No. 4 July 2000.
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF(n.EQ.0.OR.k.EQ.0) THEN
         RETURN
      END IF
!
!     Base case
!
      IF(n.EQ.1.OR.k.EQ.1) THEN
         t(1,1) = tau(1)
         RETURN
      END IF
!
!     Beginning of executable statements
!
      l = k / 2
!
!     Determine what kind of Q we need to compute
!     We assume that if the user doesn't provide 'F' for DIRECT,
!     then they meant to provide 'B' and if they don't provide
!     'C' for STOREV, then they meant to provide 'R'
!
      dirf = lsame(direct,'F')
      colv = lsame(storev,'C')
!
!     QR happens when we have forward direction in column storage
!
      qr = dirf.AND.colv
!
!     LQ happens when we have forward direction in row storage
!
      lq = dirf.AND.(.NOT.colv)
!
!     QL happens when we have backward direction in column storage
!
      ql = (.NOT.dirf).AND.colv
!
!     The last case is RQ. Due to how we structured this, if the
!     above 3 are false, then RQ must be true, so we never store 
!     this
!     RQ happens when we have backward direction in row storage
!     RQ = (.NOT.DIRF).AND.(.NOT.COLV)
!
      IF(qr) THEN
!
!        Break V apart into 6 components
!
!        V = |---------------|
!            |V_{1,1} 0      |
!            |V_{2,1} V_{2,2}|
!            |V_{3,1} V_{3,2}|
!            |---------------|
!
!        V_{1,1}\in\R^{l,l}      unit lower triangular
!        V_{2,1}\in\R^{k-l,l}    rectangular
!        V_{3,1}\in\R^{n-k,l}    rectangular
!        
!        V_{2,2}\in\R^{k-l,k-l}  unit lower triangular
!        V_{3,2}\in\R^{n-k,k-l}  rectangular
!
!        We will construct the T matrix 
!        T = |---------------|
!            |T_{1,1} T_{1,2}|
!            |0       T_{2,2}|
!            |---------------|
!
!        T is the triangular factor obtained from block reflectors. 
!        To motivate the structure, assume we have already computed T_{1,1}
!        and T_{2,2}. Then collect the associated reflectors in V_1 and V_2
!
!        T_{1,1}\in\R^{l, l}     upper triangular
!        T_{2,2}\in\R^{k-l, k-l} upper triangular
!        T_{1,2}\in\R^{l, k-l}   rectangular
!
!        Where l = floor(k/2)
!
!        Then, consider the product:
!        
!        (I - V_1*T_{1,1}*V_1')*(I - V_2*T_{2,2}*V_2')
!        = I - V_1*T_{1,1}*V_1' - V_2*T_{2,2}*V_2' + V_1*T_{1,1}*V_1'*V_2*T_{2,2}*V_2'
!        
!        Define T_{1,2} = -T_{1,1}*V_1'*V_2*T_{2,2}
!        
!        Then, we can define the matrix V as 
!        V = |-------|
!            |V_1 V_2|
!            |-------|
!        
!        So, our product is equivalent to the matrix product
!        I - V*T*V'
!        This means, we can compute T_{1,1} and T_{2,2}, then use this information
!        to compute T_{1,2}
!
!        Compute T_{1,1} recursively
!
         CALL dlarft(direct, storev, n, l, v, ldv, tau, t, ldt)
!
!        Compute T_{2,2} recursively
!
         CALL dlarft(direct, storev, n-l, k-l, v(l+1, l+1), ldv, &
     &               tau(l+1), t(l+1, l+1), ldt)
!
!        Compute T_{1,2} 
!        T_{1,2} = V_{2,1}'
!
         DO j = 1, l
            DO i = 1, k-l
               t(j, l+i) = v(l+i, j)
            END DO
         END DO
!
!        T_{1,2} = T_{1,2}*V_{2,2}
!
         CALL dtrmm('Right', 'Lower', 'No transpose', 'Unit', l, &
     &               k-l, one, v(l+1, l+1), ldv, t(1, l+1), ldt)

!
!        T_{1,2} = V_{3,1}'*V_{3,2} + T_{1,2}
!        Note: We assume K <= N, and GEMM will do nothing if N=K
!
         CALL dgemm('Transpose', 'No transpose', l, k-l, n-k, one, &
     &               v(k+1, 1), ldv, v(k+1, l+1), ldv, one, &
     &               t(1, l+1), ldt)
!
!        At this point, we have that T_{1,2} = V_1'*V_2
!        All that is left is to pre and post multiply by -T_{1,1} and T_{2,2}
!        respectively.
!
!        T_{1,2} = -T_{1,1}*T_{1,2}
!
         CALL dtrmm('Left', 'Upper', 'No transpose', 'Non-unit', l, &
     &               k-l, neg_one, t, ldt, t(1, l+1), ldt)
!
!        T_{1,2} = T_{1,2}*T_{2,2}
!
         CALL dtrmm('Right', 'Upper', 'No transpose', 'Non-unit', l, &
     &               k-l, one, t(l+1, l+1), ldt, t(1, l+1), ldt)

      ELSE IF(lq) THEN
!
!        Break V apart into 6 components
!
!        V = |----------------------|
!            |V_{1,1} V_{1,2} V{1,3}|
!            |0       V_{2,2} V{2,3}|
!            |----------------------|
!
!        V_{1,1}\in\R^{l,l}      unit upper triangular
!        V_{1,2}\in\R^{l,k-l}    rectangular
!        V_{1,3}\in\R^{l,n-k}    rectangular
!        
!        V_{2,2}\in\R^{k-l,k-l}  unit upper triangular
!        V_{2,3}\in\R^{k-l,n-k}  rectangular
!
!        Where l = floor(k/2)
!
!        We will construct the T matrix 
!        T = |---------------|
!            |T_{1,1} T_{1,2}|
!            |0       T_{2,2}|
!            |---------------|
!
!        T is the triangular factor obtained from block reflectors. 
!        To motivate the structure, assume we have already computed T_{1,1}
!        and T_{2,2}. Then collect the associated reflectors in V_1 and V_2
!
!        T_{1,1}\in\R^{l, l}     upper triangular
!        T_{2,2}\in\R^{k-l, k-l} upper triangular
!        T_{1,2}\in\R^{l, k-l}   rectangular
!
!        Then, consider the product:
!        
!        (I - V_1'*T_{1,1}*V_1)*(I - V_2'*T_{2,2}*V_2)
!        = I - V_1'*T_{1,1}*V_1 - V_2'*T_{2,2}*V_2 + V_1'*T_{1,1}*V_1*V_2'*T_{2,2}*V_2
!        
!        Define T_{1,2} = -T_{1,1}*V_1*V_2'*T_{2,2}
!        
!        Then, we can define the matrix V as 
!        V = |---|
!            |V_1|
!            |V_2|
!            |---|
!        
!        So, our product is equivalent to the matrix product
!        I - V'*T*V
!        This means, we can compute T_{1,1} and T_{2,2}, then use this information
!        to compute T_{1,2}
!
!        Compute T_{1,1} recursively
!
         CALL dlarft(direct, storev, n, l, v, ldv, tau, t, ldt)
!
!        Compute T_{2,2} recursively
!
         CALL dlarft(direct, storev, n-l, k-l, v(l+1, l+1), ldv, &
     &      tau(l+1), t(l+1, l+1), ldt)

!
!        Compute T_{1,2}
!        T_{1,2} = V_{1,2}
!
         CALL dlacpy('All', l, k-l, v(1, l+1), ldv, t(1, l+1), ldt)
!
!        T_{1,2} = T_{1,2}*V_{2,2}'
!
         CALL dtrmm('Right', 'Upper', 'Transpose', 'Unit', l, k-l, &
     &               one, v(l+1, l+1), ldv, t(1, l+1), ldt)

!
!        T_{1,2} = V_{1,3}*V_{2,3}' + T_{1,2}
!        Note: We assume K <= N, and GEMM will do nothing if N=K
!
         CALL dgemm('No transpose', 'Transpose', l, k-l, n-k, one, &
     &               v(1, k+1), ldv, v(l+1, k+1), ldv, one, &
     &               t(1, l+1), ldt)
!
!        At this point, we have that T_{1,2} = V_1*V_2'
!        All that is left is to pre and post multiply by -T_{1,1} and T_{2,2}
!        respectively.
!
!        T_{1,2} = -T_{1,1}*T_{1,2}
!
         CALL dtrmm('Left', 'Upper', 'No transpose', 'Non-unit', l, &
     &               k-l, neg_one, t, ldt, t(1, l+1), ldt)

!
!        T_{1,2} = T_{1,2}*T_{2,2}
!
         CALL dtrmm('Right', 'Upper', 'No transpose', 'Non-unit', l, &
     &               k-l, one, t(l+1, l+1), ldt, t(1, l+1), ldt)
      ELSE IF(ql) THEN
!
!        Break V apart into 6 components
!
!        V = |---------------|
!            |V_{1,1} V_{1,2}|
!            |V_{2,1} V_{2,2}|
!            |0       V_{3,2}|
!            |---------------|
!
!        V_{1,1}\in\R^{n-k,k-l}  rectangular
!        V_{2,1}\in\R^{k-l,k-l}  unit upper triangular
!        
!        V_{1,2}\in\R^{n-k,l}    rectangular
!        V_{2,2}\in\R^{k-l,l}    rectangular
!        V_{3,2}\in\R^{l,l}      unit upper triangular
!
!        We will construct the T matrix 
!        T = |---------------|
!            |T_{1,1} 0      |
!            |T_{2,1} T_{2,2}|
!            |---------------|
!
!        T is the triangular factor obtained from block reflectors. 
!        To motivate the structure, assume we have already computed T_{1,1}
!        and T_{2,2}. Then collect the associated reflectors in V_1 and V_2
!
!        T_{1,1}\in\R^{k-l, k-l} non-unit lower triangular
!        T_{2,2}\in\R^{l, l}     non-unit lower triangular
!        T_{2,1}\in\R^{k-l, l}   rectangular
!
!        Where l = floor(k/2)
!
!        Then, consider the product:
!        
!        (I - V_2*T_{2,2}*V_2')*(I - V_1*T_{1,1}*V_1')
!        = I - V_2*T_{2,2}*V_2' - V_1*T_{1,1}*V_1' + V_2*T_{2,2}*V_2'*V_1*T_{1,1}*V_1'
!        
!        Define T_{2,1} = -T_{2,2}*V_2'*V_1*T_{1,1}
!        
!        Then, we can define the matrix V as 
!        V = |-------|
!            |V_1 V_2|
!            |-------|
!        
!        So, our product is equivalent to the matrix product
!        I - V*T*V'
!        This means, we can compute T_{1,1} and T_{2,2}, then use this information
!        to compute T_{2,1}
!
!        Compute T_{1,1} recursively
!
         CALL dlarft(direct, storev, n-l, k-l, v, ldv, tau, t, ldt)
!
!        Compute T_{2,2} recursively
!
         CALL dlarft(direct, storev, n, l, v(1, k-l+1), ldv, &
     &               tau(k-l+1), t(k-l+1, k-l+1), ldt)
!
!        Compute T_{2,1}
!        T_{2,1} = V_{2,2}'
!
         DO j = 1, k-l
            DO i = 1, l
               t(k-l+i, j) = v(n-k+j, k-l+i)
            END DO
         END DO
!
!        T_{2,1} = T_{2,1}*V_{2,1}
!
         CALL dtrmm('Right', 'Upper', 'No transpose', 'Unit', l, &
     &               k-l, one, v(n-k+1, 1), ldv, t(k-l+1, 1), ldt)

!
!        T_{2,1} = V_{2,2}'*V_{2,1} + T_{2,1}
!        Note: We assume K <= N, and GEMM will do nothing if N=K
!
         CALL dgemm('Transpose', 'No transpose', l, k-l, n-k, one, &
     &               v(1, k-l+1), ldv, v, ldv, one, t(k-l+1, 1), &
     &               ldt)
!
!        At this point, we have that T_{2,1} = V_2'*V_1
!        All that is left is to pre and post multiply by -T_{2,2} and T_{1,1}
!        respectively.
!
!        T_{2,1} = -T_{2,2}*T_{2,1}
!
         CALL dtrmm('Left', 'Lower', 'No transpose', 'Non-unit', l, &
     &               k-l, neg_one, t(k-l+1, k-l+1), ldt, &
     &               t(k-l+1, 1), ldt)
!
!        T_{2,1} = T_{2,1}*T_{1,1}
!
         CALL dtrmm('Right', 'Lower', 'No transpose', 'Non-unit', l, &
     &               k-l, one, t, ldt, t(k-l+1, 1), ldt)
      ELSE
!
!        Else means RQ case
!
!        Break V apart into 6 components
!
!        V = |-----------------------|
!            |V_{1,1} V_{1,2} 0      |
!            |V_{2,1} V_{2,2} V_{2,3}|
!            |-----------------------|
!
!        V_{1,1}\in\R^{k-l,n-k}  rectangular
!        V_{1,2}\in\R^{k-l,k-l}  unit lower triangular
!
!        V_{2,1}\in\R^{l,n-k}    rectangular
!        V_{2,2}\in\R^{l,k-l}    rectangular
!        V_{2,3}\in\R^{l,l}      unit lower triangular
!
!        We will construct the T matrix 
!        T = |---------------|
!            |T_{1,1} 0      |
!            |T_{2,1} T_{2,2}|
!            |---------------|
!
!        T is the triangular factor obtained from block reflectors. 
!        To motivate the structure, assume we have already computed T_{1,1}
!        and T_{2,2}. Then collect the associated reflectors in V_1 and V_2
!
!        T_{1,1}\in\R^{k-l, k-l} non-unit lower triangular
!        T_{2,2}\in\R^{l, l}     non-unit lower triangular
!        T_{2,1}\in\R^{k-l, l}   rectangular
!
!        Where l = floor(k/2)
!
!        Then, consider the product:
!        
!        (I - V_2'*T_{2,2}*V_2)*(I - V_1'*T_{1,1}*V_1)
!        = I - V_2'*T_{2,2}*V_2 - V_1'*T_{1,1}*V_1 + V_2'*T_{2,2}*V_2*V_1'*T_{1,1}*V_1
!        
!        Define T_{2,1} = -T_{2,2}*V_2*V_1'*T_{1,1}
!        
!        Then, we can define the matrix V as 
!        V = |---|
!            |V_1|
!            |V_2|
!            |---|
!        
!        So, our product is equivalent to the matrix product
!        I - V'*T*V
!        This means, we can compute T_{1,1} and T_{2,2}, then use this information
!        to compute T_{2,1}
!
!        Compute T_{1,1} recursively
!
         CALL dlarft(direct, storev, n-l, k-l, v, ldv, tau, t, ldt)
!
!        Compute T_{2,2} recursively
!
         CALL dlarft(direct, storev, n, l, v(k-l+1, 1), ldv, &
     &               tau(k-l+1), t(k-l+1, k-l+1), ldt)
!
!        Compute T_{2,1}
!        T_{2,1} = V_{2,2}
!
         CALL dlacpy('All', l, k-l, v(k-l+1, n-k+1), ldv, &
     &               t(k-l+1, 1), ldt)

!
!        T_{2,1} = T_{2,1}*V_{1,2}'
!
         CALL dtrmm('Right', 'Lower', 'Transpose', 'Unit', l, k-l, &
     &               one, v(1, n-k+1), ldv, t(k-l+1, 1), ldt)

!
!        T_{2,1} = V_{2,1}*V_{1,1}' + T_{2,1} 
!        Note: We assume K <= N, and GEMM will do nothing if N=K
!
         CALL dgemm('No transpose', 'Transpose', l, k-l, n-k, one, &
     &               v(k-l+1, 1), ldv, v, ldv, one, t(k-l+1, 1), &
     &               ldt)

!
!        At this point, we have that T_{2,1} = V_2*V_1'
!        All that is left is to pre and post multiply by -T_{2,2} and T_{1,1}
!        respectively.
!
!        T_{2,1} = -T_{2,2}*T_{2,1}
!
         CALL dtrmm('Left', 'Lower', 'No tranpose', 'Non-unit', l, &
     &               k-l, neg_one, t(k-l+1, k-l+1), ldt, &
     &               t(k-l+1, 1), ldt)

!
!        T_{2,1} = T_{2,1}*T_{1,1}
!
         CALL dtrmm('Right', 'Lower', 'No tranpose', 'Non-unit', l, &
     &               k-l, one, t, ldt, t(k-l+1, 1), ldt)
      END IF

      END SUBROUTINE
!> \brief \b DGEMM
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION ALPHA,BETA
!       INTEGER K,LDA,LDB,LDC,M,N
!       CHARACTER TRANSA,TRANSB
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGEMM  performs one of the matrix-matrix operations
!>
!>    C := alpha*op( A )*op( B ) + beta*C,
!>
!> where  op( X ) is one of
!>
!>    op( X ) = X   or   op( X ) = X**T,
!>
!> alpha and beta are scalars, and A, B and C are matrices, with op( A )
!> an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANSA
!> \verbatim
!>          TRANSA is CHARACTER*1
!>           On entry, TRANSA specifies the form of op( A ) to be used in
!>           the matrix multiplication as follows:
!>
!>              TRANSA = 'N' or 'n',  op( A ) = A.
!>
!>              TRANSA = 'T' or 't',  op( A ) = A**T.
!>
!>              TRANSA = 'C' or 'c',  op( A ) = A**T.
!> \endverbatim
!>
!> \param[in] TRANSB
!> \verbatim
!>          TRANSB is CHARACTER*1
!>           On entry, TRANSB specifies the form of op( B ) to be used in
!>           the matrix multiplication as follows:
!>
!>              TRANSB = 'N' or 'n',  op( B ) = B.
!>
!>              TRANSB = 'T' or 't',  op( B ) = B**T.
!>
!>              TRANSB = 'C' or 'c',  op( B ) = B**T.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>           On entry,  M  specifies  the number  of rows  of the  matrix
!>           op( A )  and of the  matrix  C.  M  must  be at least  zero.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry,  N  specifies the number  of columns of the matrix
!>           op( B ) and the number of columns of the matrix C. N must be
!>           at least zero.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>           On entry,  K  specifies  the number of columns of the matrix
!>           op( A ) and the number of rows of the matrix op( B ). K must
!>           be at least  zero.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is DOUBLE PRECISION.
!>           On entry, ALPHA specifies the scalar alpha.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension ( LDA, ka ), where ka is
!>           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!>           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!>           part of the array  A  must contain the matrix  A,  otherwise
!>           the leading  k by m  part of the array  A  must contain  the
!>           matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!>           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!>           least  max( 1, k ).
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension ( LDB, kb ), where kb is
!>           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!>           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!>           part of the array  B  must contain the matrix  B,  otherwise
!>           the leading  n by k  part of the array  B  must contain  the
!>           matrix B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>           On entry, LDB specifies the first dimension of B as declared
!>           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!>           LDB must be at least  max( 1, k ), otherwise  LDB must be at
!>           least  max( 1, n ).
!> \endverbatim
!>
!> \param[in] BETA
!> \verbatim
!>          BETA is DOUBLE PRECISION.
!>           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!>           supplied as zero then C need not be set on input.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension ( LDC, N )
!>           Before entry, the leading  m by n  part of the array  C must
!>           contain the matrix  C,  except when  beta  is zero, in which
!>           case C need not be set on entry.
!>           On exit, the array  C  is overwritten by the  m by n  matrix
!>           ( alpha*op( A )*op( B ) + beta*C ).
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>           On entry, LDC specifies the first dimension of C as declared
!>           in  the  calling  (sub)  program.   LDC  must  be  at  least
!>           max( 1, m ).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup gemm
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Level 3 Blas routine.
!>
!>  -- Written on 8-February-1989.
!>     Jack Dongarra, Argonne National Laboratory.
!>     Iain Duff, AERE Harwell.
!>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!>     Sven Hammarling, Numerical Algorithms Group Ltd.
!> \endverbatim
!>
!  =====================================================================

      SUBROUTINE dgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB, &
     &         BETA,C,LDC)
!
!  -- Reference BLAS level3 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BETA
      INTEGER K,LDA,LDB,LDC,M,N
      CHARACTER TRANSA,TRANSB
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
!     ..
!
!  =====================================================================
!
!     .. External Functions ..
!      LOGICAL LSAME
!      EXTERNAL lsame
!     ..
!     .. External Subroutines ..
!      EXTERNAL xerbla
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC max
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,J,L,NROWA,NROWB
      LOGICAL NOTA,NOTB
!     ..
!     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      parameter(one=1.0d+0,zero=0.0d+0)
!     ..
!
!     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!     transposed and set  NROWA and NROWB  as the number of rows of  A
!     and  B  respectively.
!
      nota = lsame(transa,'N')
      notb = lsame(transb,'N')
      IF (nota) THEN
          nrowa = m
      ELSE
          nrowa = k
      END IF
      IF (notb) THEN
          nrowb = k
      ELSE
          nrowb = n
      END IF
!
!     Test the input parameters.
!
      info = 0
      IF ((.NOT.nota) .AND. (.NOT.lsame(transa,'C')) .AND. &
     &    (.NOT.lsame(transa,'T'))) THEN
          info = 1
      ELSE IF ((.NOT.notb) .AND. (.NOT.lsame(transb,'C')) .AND. &
     &         (.NOT.lsame(transb,'T'))) THEN
          info = 2
      ELSE IF (m.LT.0) THEN
          info = 3
      ELSE IF (n.LT.0) THEN
          info = 4
      ELSE IF (k.LT.0) THEN
          info = 5
      ELSE IF (lda.LT.max(1,nrowa)) THEN
          info = 8
      ELSE IF (ldb.LT.max(1,nrowb)) THEN
          info = 10
      ELSE IF (ldc.LT.max(1,m)) THEN
          info = 13
      END IF
      IF (info.NE.0) THEN
          CALL xerbla('DGEMM ',info)
          RETURN
      END IF
!
!     Quick return if possible.
!
      IF ((m.EQ.0) .OR. (n.EQ.0) .OR. &
     &    (((alpha.EQ.zero).OR. (k.EQ.0)).AND. (beta.EQ.one))) RETURN
!
!     And if  alpha.eq.zero.
!
      IF (alpha.EQ.zero) THEN
          IF (beta.EQ.zero) THEN
              DO 20 j = 1,n
                  DO 10 i = 1,m
                      c(i,j) = zero
10                CONTINUE
20            CONTINUE
          ELSE
              DO 40 j = 1,n
                  DO 30 i = 1,m
                      c(i,j) = beta*c(i,j)
30                CONTINUE
40            CONTINUE
          END IF
          RETURN
      END IF
!
!     Start the operations.
!
      IF (notb) THEN
          IF (nota) THEN
!
!           Form  C := alpha*A*B + beta*C.
!
              DO 90 j = 1,n
                  IF (beta.EQ.zero) THEN
                      DO 50 i = 1,m
                          c(i,j) = zero
50                    CONTINUE
                  ELSE IF (beta.NE.one) THEN
                      DO 60 i = 1,m
                          c(i,j) = beta*c(i,j)
60                    CONTINUE
                  END IF
                  DO 80 l = 1,k
                      temp = alpha*b(l,j)
                      DO 70 i = 1,m
                          c(i,j) = c(i,j) + temp*a(i,l)
70                    CONTINUE
80                CONTINUE
90            CONTINUE
          ELSE
!
!           Form  C := alpha*A**T*B + beta*C
!
              DO 120 j = 1,n
                  DO 110 i = 1,m
                      temp = zero
                      DO 100 l = 1,k
                          temp = temp + a(l,i)*b(l,j)
100                   CONTINUE
                      IF (beta.EQ.zero) THEN
                          c(i,j) = alpha*temp
                      ELSE
                          c(i,j) = alpha*temp + beta*c(i,j)
                      END IF
110               CONTINUE
120           CONTINUE
          END IF
      ELSE
          IF (nota) THEN
!
!           Form  C := alpha*A*B**T + beta*C
!
              DO 170 j = 1,n
                  IF (beta.EQ.zero) THEN
                      DO 130 i = 1,m
                          c(i,j) = zero
130                   CONTINUE
                  ELSE IF (beta.NE.one) THEN
                      DO 140 i = 1,m
                          c(i,j) = beta*c(i,j)
140                   CONTINUE
                  END IF
                  DO 160 l = 1,k
                      temp = alpha*b(j,l)
                      DO 150 i = 1,m
                          c(i,j) = c(i,j) + temp*a(i,l)
150                   CONTINUE
160               CONTINUE
170           CONTINUE
          ELSE
!
!           Form  C := alpha*A**T*B**T + beta*C
!
              DO 200 j = 1,n
                  DO 190 i = 1,m
                      temp = zero
                      DO 180 l = 1,k
                          temp = temp + a(l,i)*b(j,l)
180                   CONTINUE
                      IF (beta.EQ.zero) THEN
                          c(i,j) = alpha*temp
                      ELSE
                          c(i,j) = alpha*temp + beta*c(i,j)
                      END IF
190               CONTINUE
200           CONTINUE
          END IF
      END IF
!
      RETURN
!
!     End of DGEMM
!

      END
!> \brief \b DLACPY copies all or part of one two-dimensional array to another.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download DLACPY + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlacpy.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlacpy.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlacpy.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLACPY( UPLO, M, N, A, LDA, B, LDB )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            LDA, LDB, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLACPY copies all or part of a two-dimensional matrix A to another
!> matrix B.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies the part of the matrix A to be copied to B.
!>          = 'U':      Upper triangular part
!>          = 'L':      Lower triangular part
!>          Otherwise:  All of the matrix A
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          The m by n matrix A.  If UPLO = 'U', only the upper triangle
!>          or trapezoid is accessed; if UPLO = 'L', only the lower
!>          triangle or trapezoid is accessed.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (LDB,N)
!>          On exit, B = A in the locations specified by UPLO.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,M).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup lacpy
!
!  =====================================================================

      SUBROUTINE dlacpy( UPLO, M, N, A, LDA, B, LDB )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDB, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, J
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      EXTERNAL           lsame
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          min
!     ..
!     .. Executable Statements ..
!
      IF( lsame( uplo, 'U' ) ) THEN
         DO 20 j = 1, n
            DO 10 i = 1, min( j, m )
               b( i, j ) = a( i, j )
10          CONTINUE
20       CONTINUE
      ELSE IF( lsame( uplo, 'L' ) ) THEN
         DO 40 j = 1, n
            DO 30 i = j, m
               b( i, j ) = a( i, j )
30          CONTINUE
40       CONTINUE
      ELSE
         DO 60 j = 1, n
            DO 50 i = 1, m
               b( i, j ) = a( i, j )
50          CONTINUE
60       CONTINUE
      END IF
      RETURN
!
!     End of DLACPY
!

      END
!> \brief \b DORMHR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download DORMHR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dormhr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dormhr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dormhr.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       SUBROUTINE DORMHR( SIDE, TRANS, M, N, ILO, IHI, A, LDA, TAU, C,
!                          LDC, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          SIDE, TRANS
!       INTEGER            IHI, ILO, INFO, LDA, LDC, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DORMHR overwrites the general real M-by-N matrix C with
!>
!>                 SIDE = 'L'     SIDE = 'R'
!> TRANS = 'N':      Q * C          C * Q
!> TRANS = 'T':      Q**T * C       C * Q**T
!>
!> where Q is a real orthogonal matrix of order nq, with nq = m if
!> SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of
!> IHI-ILO elementary reflectors, as returned by DGEHRD:
!>
!> Q = H(ilo) H(ilo+1) . . . H(ihi-1).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'L': apply Q or Q**T from the Left;
!>          = 'R': apply Q or Q**T from the Right.
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          = 'N':  No transpose, apply Q;
!>          = 'T':  Transpose, apply Q**T.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix C. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix C. N >= 0.
!> \endverbatim
!>
!> \param[in] ILO
!> \verbatim
!>          ILO is INTEGER
!> \endverbatim
!>
!> \param[in] IHI
!> \verbatim
!>          IHI is INTEGER
!>
!>          ILO and IHI must have the same values as in the previous call
!>          of DGEHRD. Q is equal to the unit matrix except in the
!>          submatrix Q(ilo+1:ihi,ilo+1:ihi).
!>          If SIDE = 'L', then 1 <= ILO <= IHI <= M, if M > 0, and
!>          ILO = 1 and IHI = 0, if M = 0;
!>          if SIDE = 'R', then 1 <= ILO <= IHI <= N, if N > 0, and
!>          ILO = 1 and IHI = 0, if N = 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension
!>                               (LDA,M) if SIDE = 'L'
!>                               (LDA,N) if SIDE = 'R'
!>          The vectors which define the elementary reflectors, as
!>          returned by DGEHRD.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.
!>          LDA >= max(1,M) if SIDE = 'L'; LDA >= max(1,N) if SIDE = 'R'.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION array, dimension
!>                               (M-1) if SIDE = 'L'
!>                               (N-1) if SIDE = 'R'
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by DGEHRD.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension (LDC,N)
!>          On entry, the M-by-N matrix C.
!>          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C. LDC >= max(1,M).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.
!>          If SIDE = 'L', LWORK >= max(1,N);
!>          if SIDE = 'R', LWORK >= max(1,M).
!>          For optimum performance LWORK >= N*NB if SIDE = 'L', and
!>          LWORK >= M*NB if SIDE = 'R', where NB is the optimal
!>          blocksize.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup unmhr
!
!  =====================================================================

      SUBROUTINE dormhr( SIDE, TRANS, M, N, ILO, IHI, A, LDA, TAU, C, &
     &                   LDC, WORK, LWORK, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            IHI, ILO, INFO, LDA, LDC, LWORK, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            LEFT, LQUERY
      INTEGER            I1, I2, IINFO, LWKOPT, MI, NB, NH, NI, NQ, NW
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      INTEGER            ILAENV
!      EXTERNAL           lsame, ilaenv
!     ..
!     .. External Subroutines ..
!      EXTERNAL           dormqr, xerbla
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          max, min
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      info = 0
      nh = ihi - ilo
      left = lsame( side, 'L' )
      lquery = ( lwork.EQ.-1 )
!
!     NQ is the order of Q and NW is the minimum dimension of WORK
!
      IF( left ) THEN
         nq = m
         nw = max( 1, n )
      ELSE
         nq = n
         nw = max( 1, m )
      END IF
      IF( .NOT.left .AND. .NOT.lsame( side, 'R' ) ) THEN
         info = -1
      ELSE IF( .NOT.lsame( trans, 'N' ) .AND. &
     &         .NOT.lsame( trans, 'T' ) ) &
     &          THEN
         info = -2
      ELSE IF( m.LT.0 ) THEN
         info = -3
      ELSE IF( n.LT.0 ) THEN
         info = -4
      ELSE IF( ilo.LT.1 .OR. ilo.GT.max( 1, nq ) ) THEN
         info = -5
      ELSE IF( ihi.LT.min( ilo, nq ) .OR. ihi.GT.nq ) THEN
         info = -6
      ELSE IF( lda.LT.max( 1, nq ) ) THEN
         info = -8
      ELSE IF( ldc.LT.max( 1, m ) ) THEN
         info = -11
      ELSE IF( lwork.LT.nw .AND. .NOT.lquery ) THEN
         info = -13
      END IF
!
      IF( info.EQ.0 ) THEN
         IF( left ) THEN
            nb = ilaenv( 1, 'DORMQR', side // trans, nh, n, nh, -1 )
         ELSE
            nb = ilaenv( 1, 'DORMQR', side // trans, m, nh, nh, -1 )
         END IF
         lwkopt = nw*nb
         work( 1 ) = lwkopt
      END IF
!
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DORMHR', -info )
         RETURN
      ELSE IF( lquery ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( m.EQ.0 .OR. n.EQ.0 .OR. nh.EQ.0 ) THEN
         work( 1 ) = 1
         RETURN
      END IF
!
      IF( left ) THEN
         mi = nh
         ni = n
         i1 = ilo + 1
         i2 = 1
      ELSE
         mi = m
         ni = nh
         i1 = 1
         i2 = ilo + 1
      END IF
!
      CALL dormqr( side, trans, mi, ni, nh, a( ilo+1, ilo ), lda, &
     &             tau( ilo ), c( i1, i2 ), ldc, work, lwork, iinfo )
!
      work( 1 ) = lwkopt
      RETURN
!
!     End of DORMHR
!

      END
!> \brief \b DTREXC
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download DTREXC + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtrexc.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtrexc.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtrexc.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       SUBROUTINE DTREXC( COMPQ, N, T, LDT, Q, LDQ, IFST, ILST, WORK,
!                          INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          COMPQ
!       INTEGER            IFST, ILST, INFO, LDQ, LDT, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   Q( LDQ, * ), T( LDT, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DTREXC reorders the real Schur factorization of a real matrix
!> A = Q*T*Q**T, so that the diagonal block of T with row index IFST is
!> moved to row ILST.
!>
!> The real Schur form T is reordered by an orthogonal similarity
!> transformation Z**T*T*Z, and optionally the matrix Q of Schur vectors
!> is updated by postmultiplying it with Z.
!>
!> T must be in Schur canonical form (as returned by DHSEQR), that is,
!> block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each
!> 2-by-2 diagonal block has its diagonal elements equal and its
!> off-diagonal elements of opposite sign.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] COMPQ
!> \verbatim
!>          COMPQ is CHARACTER*1
!>          = 'V':  update the matrix Q of Schur vectors;
!>          = 'N':  do not update Q.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix T. N >= 0.
!>          If N == 0 arguments ILST and IFST may be any value.
!> \endverbatim
!>
!> \param[in,out] T
!> \verbatim
!>          T is DOUBLE PRECISION array, dimension (LDT,N)
!>          On entry, the upper quasi-triangular matrix T, in Schur
!>          Schur canonical form.
!>          On exit, the reordered upper quasi-triangular matrix, again
!>          in Schur canonical form.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of the array T. LDT >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] Q
!> \verbatim
!>          Q is DOUBLE PRECISION array, dimension (LDQ,N)
!>          On entry, if COMPQ = 'V', the matrix Q of Schur vectors.
!>          On exit, if COMPQ = 'V', Q has been postmultiplied by the
!>          orthogonal transformation matrix Z which reorders T.
!>          If COMPQ = 'N', Q is not referenced.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of the array Q.  LDQ >= 1, and if
!>          COMPQ = 'V', LDQ >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] IFST
!> \verbatim
!>          IFST is INTEGER
!> \endverbatim
!>
!> \param[in,out] ILST
!> \verbatim
!>          ILST is INTEGER
!>
!>          Specify the reordering of the diagonal blocks of T.
!>          The block with row index IFST is moved to row ILST, by a
!>          sequence of transpositions between adjacent blocks.
!>          On exit, if IFST pointed on entry to the second row of a
!>          2-by-2 block, it is changed to point to the first row; ILST
!>          always points to the first row of the block in its final
!>          position (which may differ from its input value by +1 or -1).
!>          1 <= IFST <= N; 1 <= ILST <= N.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          = 1:  two adjacent blocks were too close to swap (the problem
!>                is very ill-conditioned); T may have been partially
!>                reordered, and ILST points to the first row of the
!>                current position of the block being moved.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup trexc
!
!  =====================================================================

      SUBROUTINE dtrexc( COMPQ, N, T, LDT, Q, LDQ, IFST, ILST, WORK, &
     &                   INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          COMPQ
      INTEGER            IFST, ILST, INFO, LDQ, LDT, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   Q( LDQ, * ), T( LDT, * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      parameter( zero = 0.0d+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            WANTQ
      INTEGER            HERE, NBF, NBL, NBNEXT
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      EXTERNAL           lsame
!     ..
!     .. External Subroutines ..
!      EXTERNAL           dlaexc, xerbla
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          max
!     ..
!     .. Executable Statements ..
!
!     Decode and test the input arguments.
!
      info = 0
      wantq = lsame( compq, 'V' )
      IF( .NOT.wantq .AND. .NOT.lsame( compq, 'N' ) ) THEN
         info = -1
      ELSE IF( n.LT.0 ) THEN
         info = -2
      ELSE IF( ldt.LT.max( 1, n ) ) THEN
         info = -4
      ELSE IF( ldq.LT.1 .OR. ( wantq .AND. ldq.LT.max( 1, n ) ) ) THEN
         info = -6
      ELSE IF(( ifst.LT.1 .OR. ifst.GT.n ).AND.( n.GT.0 )) THEN
         info = -7
      ELSE IF(( ilst.LT.1 .OR. ilst.GT.n ).AND.( n.GT.0 )) THEN
         info = -8
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DTREXC', -info )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( n.LE.1 ) &
     &   RETURN
!
!     Determine the first row of specified block
!     and find out it is 1 by 1 or 2 by 2.
!
      IF( ifst.GT.1 ) THEN
         IF( t( ifst, ifst-1 ).NE.zero ) &
     &      ifst = ifst - 1
      END IF
      nbf = 1
      IF( ifst.LT.n ) THEN
         IF( t( ifst+1, ifst ).NE.zero ) &
     &      nbf = 2
      END IF
!
!     Determine the first row of the final block
!     and find out it is 1 by 1 or 2 by 2.
!
      IF( ilst.GT.1 ) THEN
         IF( t( ilst, ilst-1 ).NE.zero ) &
     &      ilst = ilst - 1
      END IF
      nbl = 1
      IF( ilst.LT.n ) THEN
         IF( t( ilst+1, ilst ).NE.zero ) &
     &      nbl = 2
      END IF
!
      IF( ifst.EQ.ilst ) &
     &   RETURN
!
      IF( ifst.LT.ilst ) THEN
!
!        Update ILST
!
         IF( nbf.EQ.2 .AND. nbl.EQ.1 ) &
     &      ilst = ilst - 1
         IF( nbf.EQ.1 .AND. nbl.EQ.2 ) &
     &      ilst = ilst + 1
!
         here = ifst
!
10       CONTINUE
!
!        Swap block with next one below
!
         IF( nbf.EQ.1 .OR. nbf.EQ.2 ) THEN
!
!           Current block either 1 by 1 or 2 by 2
!
            nbnext = 1
            IF( here+nbf+1.LE.n ) THEN
               IF( t( here+nbf+1, here+nbf ).NE.zero ) &
     &            nbnext = 2
            END IF
            CALL dlaexc( wantq, n, t, ldt, q, ldq, here, nbf, nbnext, &
     &                   work, info )
            IF( info.NE.0 ) THEN
               ilst = here
               RETURN
            END IF
            here = here + nbnext
!
!           Test if 2 by 2 block breaks into two 1 by 1 blocks
!
            IF( nbf.EQ.2 ) THEN
               IF( t( here+1, here ).EQ.zero ) &
     &            nbf = 3
            END IF
!
         ELSE
!
!           Current block consists of two 1 by 1 blocks each of which
!           must be swapped individually
!
            nbnext = 1
            IF( here+3.LE.n ) THEN
               IF( t( here+3, here+2 ).NE.zero ) &
     &            nbnext = 2
            END IF
            CALL dlaexc( wantq, n, t, ldt, q, ldq, here+1, 1, nbnext, &
     &                   work, info )
            IF( info.NE.0 ) THEN
               ilst = here
               RETURN
            END IF
            IF( nbnext.EQ.1 ) THEN
!
!              Swap two 1 by 1 blocks, no problems possible
!
               CALL dlaexc( wantq, n, t, ldt, q, ldq, here, 1, &
     &                      nbnext, &
     &                      work, info )
               here = here + 1
            ELSE
!
!              Recompute NBNEXT in case 2 by 2 split
!
               IF( t( here+2, here+1 ).EQ.zero ) &
     &            nbnext = 1
               IF( nbnext.EQ.2 ) THEN
!
!                 2 by 2 Block did not split
!
                  CALL dlaexc( wantq, n, t, ldt, q, ldq, here, 1, &
     &                         nbnext, work, info )
                  IF( info.NE.0 ) THEN
                     ilst = here
                     RETURN
                  END IF
                  here = here + 2
               ELSE
!
!                 2 by 2 Block did split
!
                  CALL dlaexc( wantq, n, t, ldt, q, ldq, here, 1, 1, &
     &                         work, info )
                  CALL dlaexc( wantq, n, t, ldt, q, ldq, here+1, 1, &
     &                         1, &
     &                         work, info )
                  here = here + 2
               END IF
            END IF
         END IF
         IF( here.LT.ilst ) &
     &      GO TO 10
!
      ELSE
!
         here = ifst
20       CONTINUE
!
!        Swap block with next one above
!
         IF( nbf.EQ.1 .OR. nbf.EQ.2 ) THEN
!
!           Current block either 1 by 1 or 2 by 2
!
            nbnext = 1
            IF( here.GE.3 ) THEN
               IF( t( here-1, here-2 ).NE.zero ) &
     &            nbnext = 2
            END IF
            CALL dlaexc( wantq, n, t, ldt, q, ldq, here-nbnext, &
     &                   nbnext, &
     &                   nbf, work, info )
            IF( info.NE.0 ) THEN
               ilst = here
               RETURN
            END IF
            here = here - nbnext
!
!           Test if 2 by 2 block breaks into two 1 by 1 blocks
!
            IF( nbf.EQ.2 ) THEN
               IF( t( here+1, here ).EQ.zero ) &
     &            nbf = 3
            END IF
!
         ELSE
!
!           Current block consists of two 1 by 1 blocks each of which
!           must be swapped individually
!
            nbnext = 1
            IF( here.GE.3 ) THEN
               IF( t( here-1, here-2 ).NE.zero ) &
     &            nbnext = 2
            END IF
            CALL dlaexc( wantq, n, t, ldt, q, ldq, here-nbnext, &
     &                   nbnext, &
     &                   1, work, info )
            IF( info.NE.0 ) THEN
               ilst = here
               RETURN
            END IF
            IF( nbnext.EQ.1 ) THEN
!
!              Swap two 1 by 1 blocks, no problems possible
!
               CALL dlaexc( wantq, n, t, ldt, q, ldq, here, nbnext, &
     &                      1, &
     &                      work, info )
               here = here - 1
            ELSE
!
!              Recompute NBNEXT in case 2 by 2 split
!
               IF( t( here, here-1 ).EQ.zero ) &
     &            nbnext = 1
               IF( nbnext.EQ.2 ) THEN
!
!                 2 by 2 Block did not split
!
                  CALL dlaexc( wantq, n, t, ldt, q, ldq, here-1, 2, &
     &                         1, &
     &                         work, info )
                  IF( info.NE.0 ) THEN
                     ilst = here
                     RETURN
                  END IF
                  here = here - 2
               ELSE
!
!                 2 by 2 Block did split
!
                  CALL dlaexc( wantq, n, t, ldt, q, ldq, here, 1, 1, &
     &                         work, info )
                  CALL dlaexc( wantq, n, t, ldt, q, ldq, here-1, 1, &
     &                         1, &
     &                         work, info )
                  here = here - 2
               END IF
            END IF
         END IF
         IF( here.GT.ilst ) &
     &      GO TO 20
      END IF
      ilst = here
!
      RETURN
!
!     End of DTREXC
!

      END
!> \brief \b DLAQR2 performs the orthogonal similarity transformation of a Hessenberg matrix to detect and deflate fully converged eigenvalues from a trailing principal submatrix (aggressive early deflation).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download DLAQR2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaqr2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaqr2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaqr2.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAQR2( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ,
!                          IHIZ, Z, LDZ, NS, ND, SR, SI, V, LDV, NH, T,
!                          LDT, NV, WV, LDWV, WORK, LWORK )
!
!       .. Scalar Arguments ..
!       INTEGER            IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV,
!      $                   LDZ, LWORK, N, ND, NH, NS, NV, NW
!       LOGICAL            WANTT, WANTZ
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   H( LDH, * ), SI( * ), SR( * ), T( LDT, * ),
!      $                   V( LDV, * ), WORK( * ), WV( LDWV, * ),
!      $                   Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    DLAQR2 is identical to DLAQR3 except that it avoids
!>    recursion by calling DLAHQR instead of DLAQR4.
!>
!>    Aggressive early deflation:
!>
!>    This subroutine accepts as input an upper Hessenberg matrix
!>    H and performs an orthogonal similarity transformation
!>    designed to detect and deflate fully converged eigenvalues from
!>    a trailing principal submatrix.  On output H has been over-
!>    written by a new Hessenberg matrix that is a perturbation of
!>    an orthogonal similarity transformation of H.  It is to be
!>    hoped that the final version of H has many zero subdiagonal
!>    entries.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] WANTT
!> \verbatim
!>          WANTT is LOGICAL
!>          If .TRUE., then the Hessenberg matrix H is fully updated
!>          so that the quasi-triangular Schur factor may be
!>          computed (in cooperation with the calling subroutine).
!>          If .FALSE., then only enough of H is updated to preserve
!>          the eigenvalues.
!> \endverbatim
!>
!> \param[in] WANTZ
!> \verbatim
!>          WANTZ is LOGICAL
!>          If .TRUE., then the orthogonal matrix Z is updated so
!>          so that the orthogonal Schur factor may be computed
!>          (in cooperation with the calling subroutine).
!>          If .FALSE., then Z is not referenced.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix H and (if WANTZ is .TRUE.) the
!>          order of the orthogonal matrix Z.
!> \endverbatim
!>
!> \param[in] KTOP
!> \verbatim
!>          KTOP is INTEGER
!>          It is assumed that either KTOP = 1 or H(KTOP,KTOP-1)=0.
!>          KBOT and KTOP together determine an isolated block
!>          along the diagonal of the Hessenberg matrix.
!> \endverbatim
!>
!> \param[in] KBOT
!> \verbatim
!>          KBOT is INTEGER
!>          It is assumed without a check that either
!>          KBOT = N or H(KBOT+1,KBOT)=0.  KBOT and KTOP together
!>          determine an isolated block along the diagonal of the
!>          Hessenberg matrix.
!> \endverbatim
!>
!> \param[in] NW
!> \verbatim
!>          NW is INTEGER
!>          Deflation window size.  1 <= NW <= (KBOT-KTOP+1).
!> \endverbatim
!>
!> \param[in,out] H
!> \verbatim
!>          H is DOUBLE PRECISION array, dimension (LDH,N)
!>          On input the initial N-by-N section of H stores the
!>          Hessenberg matrix undergoing aggressive early deflation.
!>          On output H has been transformed by an orthogonal
!>          similarity transformation, perturbed, and the returned
!>          to Hessenberg form that (it is to be hoped) has some
!>          zero subdiagonal entries.
!> \endverbatim
!>
!> \param[in] LDH
!> \verbatim
!>          LDH is INTEGER
!>          Leading dimension of H just as declared in the calling
!>          subroutine.  N <= LDH
!> \endverbatim
!>
!> \param[in] ILOZ
!> \verbatim
!>          ILOZ is INTEGER
!> \endverbatim
!>
!> \param[in] IHIZ
!> \verbatim
!>          IHIZ is INTEGER
!>          Specify the rows of Z to which transformations must be
!>          applied if WANTZ is .TRUE.. 1 <= ILOZ <= IHIZ <= N.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is DOUBLE PRECISION array, dimension (LDZ,N)
!>          IF WANTZ is .TRUE., then on output, the orthogonal
!>          similarity transformation mentioned above has been
!>          accumulated into Z(ILOZ:IHIZ,ILOZ:IHIZ) from the right.
!>          If WANTZ is .FALSE., then Z is unreferenced.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of Z just as declared in the
!>          calling subroutine.  1 <= LDZ.
!> \endverbatim
!>
!> \param[out] NS
!> \verbatim
!>          NS is INTEGER
!>          The number of unconverged (ie approximate) eigenvalues
!>          returned in SR and SI that may be used as shifts by the
!>          calling subroutine.
!> \endverbatim
!>
!> \param[out] ND
!> \verbatim
!>          ND is INTEGER
!>          The number of converged eigenvalues uncovered by this
!>          subroutine.
!> \endverbatim
!>
!> \param[out] SR
!> \verbatim
!>          SR is DOUBLE PRECISION array, dimension (KBOT)
!> \endverbatim
!>
!> \param[out] SI
!> \verbatim
!>          SI is DOUBLE PRECISION array, dimension (KBOT)
!>          On output, the real and imaginary parts of approximate
!>          eigenvalues that may be used for shifts are stored in
!>          SR(KBOT-ND-NS+1) through SR(KBOT-ND) and
!>          SI(KBOT-ND-NS+1) through SI(KBOT-ND), respectively.
!>          The real and imaginary parts of converged eigenvalues
!>          are stored in SR(KBOT-ND+1) through SR(KBOT) and
!>          SI(KBOT-ND+1) through SI(KBOT), respectively.
!> \endverbatim
!>
!> \param[out] V
!> \verbatim
!>          V is DOUBLE PRECISION array, dimension (LDV,NW)
!>          An NW-by-NW work array.
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is INTEGER
!>          The leading dimension of V just as declared in the
!>          calling subroutine.  NW <= LDV
!> \endverbatim
!>
!> \param[in] NH
!> \verbatim
!>          NH is INTEGER
!>          The number of columns of T.  NH >= NW.
!> \endverbatim
!>
!> \param[out] T
!> \verbatim
!>          T is DOUBLE PRECISION array, dimension (LDT,NW)
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of T just as declared in the
!>          calling subroutine.  NW <= LDT
!> \endverbatim
!>
!> \param[in] NV
!> \verbatim
!>          NV is INTEGER
!>          The number of rows of work array WV available for
!>          workspace.  NV >= NW.
!> \endverbatim
!>
!> \param[out] WV
!> \verbatim
!>          WV is DOUBLE PRECISION array, dimension (LDWV,NW)
!> \endverbatim
!>
!> \param[in] LDWV
!> \verbatim
!>          LDWV is INTEGER
!>          The leading dimension of W just as declared in the
!>          calling subroutine.  NW <= LDV
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (LWORK)
!>          On exit, WORK(1) is set to an estimate of the optimal value
!>          of LWORK for the given values of N, NW, KTOP and KBOT.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the work array WORK.  LWORK = 2*NW
!>          suffices, but greater efficiency may result from larger
!>          values of LWORK.
!>
!>          If LWORK = -1, then a workspace query is assumed; DLAQR2
!>          only estimates the optimal workspace size for the given
!>          values of N, NW, KTOP and KBOT.  The estimate is returned
!>          in WORK(1).  No error message related to LWORK is issued
!>          by XERBLA.  Neither H nor Z are accessed.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup laqr2
!
!> \par Contributors:
!  ==================
!>
!>       Karen Braman and Ralph Byers, Department of Mathematics,
!>       University of Kansas, USA
!>
!  =====================================================================

      SUBROUTINE dlaqr2( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, &
     &                   ILOZ, &
     &                   IHIZ, Z, LDZ, NS, ND, SR, SI, V, LDV, NH, T, &
     &                   LDT, NV, WV, LDWV, WORK, LWORK )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV, &
     &                   LDZ, LWORK, N, ND, NH, NS, NV, NW
      LOGICAL            WANTT, WANTZ
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   H( LDH, * ), SI( * ), SR( * ), T( LDT, * ), &
     &                   V( LDV, * ), WORK( * ), WV( LDWV, * ), &
     &                   Z( LDZ, * )
!     ..
!
!  ================================================================
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0d0, one = 1.0d0 )
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION   AA, BB, BETA, CC, CS, DD, EVI, EVK, FOO, S, &
     &                   SAFMAX, SAFMIN, SMLNUM, SN, TAU, ULP
      INTEGER            I, IFST, ILST, INFO, INFQR, J, JW, K, KCOL, &
     &                   kend, kln, krow, kwtop, ltop, lwk1, lwk2, &
     &                   lwkopt
      LOGICAL            BULGE, SORTED
!     ..
!     .. External Functions ..
!      DOUBLE PRECISION   DLAMCH
!      EXTERNAL           DLAMCH
!     ..
!     .. External Subroutines ..
!      EXTERNAL           dcopy, dgehrd, dgemm, dlacpy, &
!     &                   dlahqr, &
!     &                   dlanv2, dlarf1f, dlarfg, dlaset, dormhr, dtrexc
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          abs, dble, int, max, min, sqrt
!     ..
!     .. Executable Statements ..
!
!     ==== Estimate optimal workspace. ====
!
      jw = min( nw, kbot-ktop+1 )
      IF( jw.LE.2 ) THEN
         lwkopt = 1
      ELSE
!
!        ==== Workspace query call to DGEHRD ====
!
         CALL dgehrd( jw, 1, jw-1, t, ldt, work, work, -1, info )
         lwk1 = int( work( 1 ) )
!
!        ==== Workspace query call to DORMHR ====
!
         CALL dormhr( 'R', 'N', jw, jw, 1, jw-1, t, ldt, work, v, &
     &                ldv, &
     &                work, -1, info )
         lwk2 = int( work( 1 ) )
!
!        ==== Optimal workspace ====
!
         lwkopt = jw + max( lwk1, lwk2 )
      END IF
!
!     ==== Quick return in case of workspace query. ====
!
      IF( lwork.EQ.-1 ) THEN
         work( 1 ) = dble( lwkopt )
         RETURN
      END IF
!
!     ==== Nothing to do ...
!     ... for an empty active block ... ====
      ns = 0
      nd = 0
      work( 1 ) = one
      IF( ktop.GT.kbot ) &
     &   RETURN
!     ... nor for an empty deflation window. ====
      IF( nw.LT.1 ) &
     &   RETURN
!
!     ==== Machine constants ====
!
      safmin = dlamch( 'SAFE MINIMUM' )
      safmax = one / safmin
      ulp = dlamch( 'PRECISION' )
      smlnum = safmin*( dble( n ) / ulp )
!
!     ==== Setup deflation window ====
!
      jw = min( nw, kbot-ktop+1 )
      kwtop = kbot - jw + 1
      IF( kwtop.EQ.ktop ) THEN
         s = zero
      ELSE
         s = h( kwtop, kwtop-1 )
      END IF
!
      IF( kbot.EQ.kwtop ) THEN
!
!        ==== 1-by-1 deflation window: not much to do ====
!
         sr( kwtop ) = h( kwtop, kwtop )
         si( kwtop ) = zero
         ns = 1
         nd = 0
         IF( abs( s ).LE.max( smlnum, ulp*abs( h( kwtop, kwtop ) ) ) ) &
     &        THEN
            ns = 0
            nd = 1
            IF( kwtop.GT.ktop ) &
     &         h( kwtop, kwtop-1 ) = zero
         END IF
         work( 1 ) = one
         RETURN
      END IF
!
!     ==== Convert to spike-triangular form.  (In case of a
!     .    rare QR failure, this routine continues to do
!     .    aggressive early deflation using that part of
!     .    the deflation window that converged using INFQR
!     .    here and there to keep track.) ====
!
      CALL dlacpy( 'U', jw, jw, h( kwtop, kwtop ), ldh, t, ldt )
      CALL dcopy( jw-1, h( kwtop+1, kwtop ), ldh+1, t( 2, 1 ), &
     &            ldt+1 )
!
      CALL dlaset( 'A', jw, jw, zero, one, v, ldv )
      CALL dlahqr( .true., .true., jw, 1, jw, t, ldt, sr( kwtop ), &
     &             si( kwtop ), 1, jw, v, ldv, infqr )
!
!     ==== DTREXC needs a clean margin near the diagonal ====
!
      DO 10 j = 1, jw - 3
         t( j+2, j ) = zero
         t( j+3, j ) = zero
10    CONTINUE
      IF( jw.GT.2 ) &
     &   t( jw, jw-2 ) = zero
!
!     ==== Deflation detection loop ====
!
      ns = jw
      ilst = infqr + 1
20    CONTINUE
      IF( ilst.LE.ns ) THEN
         IF( ns.EQ.1 ) THEN
            bulge = .false.
         ELSE
            bulge = t( ns, ns-1 ).NE.zero
         END IF
!
!        ==== Small spike tip test for deflation ====
!
         IF( .NOT.bulge ) THEN
!
!           ==== Real eigenvalue ====
!
            foo = abs( t( ns, ns ) )
            IF( foo.EQ.zero ) &
     &         foo = abs( s )
            IF( abs( s*v( 1, ns ) ).LE.max( smlnum, ulp*foo ) ) THEN
!
!              ==== Deflatable ====
!
               ns = ns - 1
            ELSE
!
!              ==== Undeflatable.   Move it up out of the way.
!              .    (DTREXC can not fail in this case.) ====
!
               ifst = ns
               CALL dtrexc( 'V', jw, t, ldt, v, ldv, ifst, ilst, &
     &                      work, &
     &                      info )
               ilst = ilst + 1
            END IF
         ELSE
!
!           ==== Complex conjugate pair ====
!
            foo = abs( t( ns, ns ) ) + sqrt( abs( t( ns, ns-1 ) ) )* &
     &            sqrt( abs( t( ns-1, ns ) ) )
            IF( foo.EQ.zero ) &
     &         foo = abs( s )
            IF( max( abs( s*v( 1, ns ) ), abs( s*v( 1, ns-1 ) ) ).LE. &
     &          max( smlnum, ulp*foo ) ) THEN
!
!              ==== Deflatable ====
!
               ns = ns - 2
            ELSE
!
!              ==== Undeflatable. Move them up out of the way.
!              .    Fortunately, DTREXC does the right thing with
!              .    ILST in case of a rare exchange failure. ====
!
               ifst = ns
               CALL dtrexc( 'V', jw, t, ldt, v, ldv, ifst, ilst, &
     &                      work, &
     &                      info )
               ilst = ilst + 2
            END IF
         END IF
!
!        ==== End deflation detection loop ====
!
         GO TO 20
      END IF
!
!        ==== Return to Hessenberg form ====
!
      IF( ns.EQ.0 ) &
     &   s = zero
!
      IF( ns.LT.jw ) THEN
!
!        ==== sorting diagonal blocks of T improves accuracy for
!        .    graded matrices.  Bubble sort deals well with
!        .    exchange failures. ====
!
         sorted = .false.
         i = ns + 1
30       CONTINUE
         IF( sorted ) &
     &      GO TO 50
         sorted = .true.
!
         kend = i - 1
         i = infqr + 1
         IF( i.EQ.ns ) THEN
            k = i + 1
         ELSE IF( t( i+1, i ).EQ.zero ) THEN
            k = i + 1
         ELSE
            k = i + 2
         END IF
40       CONTINUE
         IF( k.LE.kend ) THEN
            IF( k.EQ.i+1 ) THEN
               evi = abs( t( i, i ) )
            ELSE
               evi = abs( t( i, i ) ) + sqrt( abs( t( i+1, i ) ) )* &
     &               sqrt( abs( t( i, i+1 ) ) )
            END IF
!
            IF( k.EQ.kend ) THEN
               evk = abs( t( k, k ) )
            ELSE IF( t( k+1, k ).EQ.zero ) THEN
               evk = abs( t( k, k ) )
            ELSE
               evk = abs( t( k, k ) ) + sqrt( abs( t( k+1, k ) ) )* &
     &               sqrt( abs( t( k, k+1 ) ) )
            END IF
!
            IF( evi.GE.evk ) THEN
               i = k
            ELSE
               sorted = .false.
               ifst = i
               ilst = k
               CALL dtrexc( 'V', jw, t, ldt, v, ldv, ifst, ilst, &
     &                      work, &
     &                      info )
               IF( info.EQ.0 ) THEN
                  i = ilst
               ELSE
                  i = k
               END IF
            END IF
            IF( i.EQ.kend ) THEN
               k = i + 1
            ELSE IF( t( i+1, i ).EQ.zero ) THEN
               k = i + 1
            ELSE
               k = i + 2
            END IF
            GO TO 40
         END IF
         GO TO 30
50       CONTINUE
      END IF
!
!     ==== Restore shift/eigenvalue array from T ====
!
      i = jw
60    CONTINUE
      IF( i.GE.infqr+1 ) THEN
         IF( i.EQ.infqr+1 ) THEN
            sr( kwtop+i-1 ) = t( i, i )
            si( kwtop+i-1 ) = zero
            i = i - 1
         ELSE IF( t( i, i-1 ).EQ.zero ) THEN
            sr( kwtop+i-1 ) = t( i, i )
            si( kwtop+i-1 ) = zero
            i = i - 1
         ELSE
            aa = t( i-1, i-1 )
            cc = t( i, i-1 )
            bb = t( i-1, i )
            dd = t( i, i )
            CALL dlanv2( aa, bb, cc, dd, sr( kwtop+i-2 ), &
     &                   si( kwtop+i-2 ), sr( kwtop+i-1 ), &
     &                   si( kwtop+i-1 ), cs, sn )
            i = i - 2
         END IF
         GO TO 60
      END IF
!
      IF( ns.LT.jw .OR. s.EQ.zero ) THEN
         IF( ns.GT.1 .AND. s.NE.zero ) THEN
!
!           ==== Reflect spike back into lower triangle ====
!
            CALL dcopy( ns, v, ldv, work, 1 )
            beta = work( 1 )
            CALL dlarfg( ns, beta, work( 2 ), 1, tau )
!
            CALL dlaset( 'L', jw-2, jw-2, zero, zero, t( 3, 1 ), &
     &                   ldt )
!
            CALL dlarf1f( 'L', ns, jw, work, 1, tau, t, ldt, &
     &                  work( jw+1 ) )
            CALL dlarf1f( 'R', ns, ns, work, 1, tau, t, ldt, &
     &                  work( jw+1 ) )
            CALL dlarf1f( 'R', jw, ns, work, 1, tau, v, ldv, &
     &                  work( jw+1 ) )
!
            CALL dgehrd( jw, 1, ns, t, ldt, work, work( jw+1 ), &
     &                   lwork-jw, info )
         END IF
!
!        ==== Copy updated reduced window into place ====
!
         IF( kwtop.GT.1 ) &
     &      h( kwtop, kwtop-1 ) = s*v( 1, 1 )
         CALL dlacpy( 'U', jw, jw, t, ldt, h( kwtop, kwtop ), ldh )
         CALL dcopy( jw-1, t( 2, 1 ), ldt+1, h( kwtop+1, kwtop ), &
     &               ldh+1 )
!
!        ==== Accumulate orthogonal matrix in order update
!        .    H and Z, if requested.  ====
!
         IF( ns.GT.1 .AND. s.NE.zero ) &
     &      CALL dormhr( 'R', 'N', jw, ns, 1, ns, t, ldt, work, v, &
     &                   ldv, &
     &                   work( jw+1 ), lwork-jw, info )
!
!        ==== Update vertical slab in H ====
!
         IF( wantt ) THEN
            ltop = 1
         ELSE
            ltop = ktop
         END IF
         DO 70 krow = ltop, kwtop - 1, nv
            kln = min( nv, kwtop-krow )
            CALL dgemm( 'N', 'N', kln, jw, jw, one, h( krow, kwtop ), &
     &                  ldh, v, ldv, zero, wv, ldwv )
            CALL dlacpy( 'A', kln, jw, wv, ldwv, h( krow, kwtop ), &
     &                   ldh )
70       CONTINUE
!
!        ==== Update horizontal slab in H ====
!
         IF( wantt ) THEN
            DO 80 kcol = kbot + 1, n, nh
               kln = min( nh, n-kcol+1 )
               CALL dgemm( 'C', 'N', jw, kln, jw, one, v, ldv, &
     &                     h( kwtop, kcol ), ldh, zero, t, ldt )
               CALL dlacpy( 'A', jw, kln, t, ldt, h( kwtop, kcol ), &
     &                      ldh )
80          CONTINUE
         END IF
!
!        ==== Update vertical slab in Z ====
!
         IF( wantz ) THEN
            DO 90 krow = iloz, ihiz, nv
               kln = min( nv, ihiz-krow+1 )
               CALL dgemm( 'N', 'N', kln, jw, jw, one, z( krow, &
     &                     kwtop ), &
     &                     ldz, v, ldv, zero, wv, ldwv )
               CALL dlacpy( 'A', kln, jw, wv, ldwv, z( krow, kwtop ), &
     &                      ldz )
90          CONTINUE
         END IF
      END IF
!
!     ==== Return the number of deflations ... ====
!
      nd = jw - ns
!
!     ==== ... and the number of shifts. (Subtracting
!     .    INFQR from the spike length takes care
!     .    of the case of a rare QR failure while
!     .    calculating eigenvalues of the deflation
!     .    window.)  ====
!
      ns = ns - infqr
!
!      ==== Return optimal workspace. ====
!
      work( 1 ) = dble( lwkopt )
!
!     ==== End of DLAQR2 ====
!

      END
!> \brief \b DORMQR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download DORMQR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dormqr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dormqr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dormqr.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       SUBROUTINE DORMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
!                          WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          SIDE, TRANS
!       INTEGER            INFO, K, LDA, LDC, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DORMQR overwrites the general real M-by-N matrix C with
!>
!>                 SIDE = 'L'     SIDE = 'R'
!> TRANS = 'N':      Q * C          C * Q
!> TRANS = 'T':      Q**T * C       C * Q**T
!>
!> where Q is a real orthogonal matrix defined as the product of k
!> elementary reflectors
!>
!>       Q = H(1) H(2) . . . H(k)
!>
!> as returned by DGEQRF. Q is of order M if SIDE = 'L' and of order N
!> if SIDE = 'R'.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'L': apply Q or Q**T from the Left;
!>          = 'R': apply Q or Q**T from the Right.
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          = 'N':  No transpose, apply Q;
!>          = 'T':  Transpose, apply Q**T.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix C. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix C. N >= 0.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of elementary reflectors whose product defines
!>          the matrix Q.
!>          If SIDE = 'L', M >= K >= 0;
!>          if SIDE = 'R', N >= K >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,K)
!>          The i-th column must contain the vector which defines the
!>          elementary reflector H(i), for i = 1,2,...,k, as returned by
!>          DGEQRF in the first k columns of its array argument A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.
!>          If SIDE = 'L', LDA >= max(1,M);
!>          if SIDE = 'R', LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by DGEQRF.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension (LDC,N)
!>          On entry, the M-by-N matrix C.
!>          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C. LDC >= max(1,M).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.
!>          If SIDE = 'L', LWORK >= max(1,N);
!>          if SIDE = 'R', LWORK >= max(1,M).
!>          For good performance, LWORK should generally be larger.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup unmqr
!
!  =====================================================================

      SUBROUTINE dormqr( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &
     &                   WORK, LWORK, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, LWORK, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            NBMAX, LDT, TSIZE
      parameter( nbmax = 64, ldt = nbmax+1, &
     &                     tsize = ldt*nbmax )
!     ..
!     .. Local Scalars ..
      LOGICAL            LEFT, LQUERY, NOTRAN
      INTEGER            I, I1, I2, I3, IB, IC, IINFO, IWT, JC, LDWORK, &
     &                   lwkopt, mi, nb, nbmin, ni, nq, nw
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      INTEGER            ILAENV
!      EXTERNAL           lsame, ilaenv
!     ..
!     .. External Subroutines ..
!      EXTERNAL           dlarfb, dlarft, dorm2r, xerbla
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          max, min
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      info = 0
      left = lsame( side, 'L' )
      notran = lsame( trans, 'N' )
      lquery = ( lwork.EQ.-1 )
!
!     NQ is the order of Q and NW is the minimum dimension of WORK
!
      IF( left ) THEN
         nq = m
         nw = max( 1, n )
      ELSE
         nq = n
         nw = max( 1, m )
      END IF
      IF( .NOT.left .AND. .NOT.lsame( side, 'R' ) ) THEN
         info = -1
      ELSE IF( .NOT.notran .AND. .NOT.lsame( trans, 'T' ) ) THEN
         info = -2
      ELSE IF( m.LT.0 ) THEN
         info = -3
      ELSE IF( n.LT.0 ) THEN
         info = -4
      ELSE IF( k.LT.0 .OR. k.GT.nq ) THEN
         info = -5
      ELSE IF( lda.LT.max( 1, nq ) ) THEN
         info = -7
      ELSE IF( ldc.LT.max( 1, m ) ) THEN
         info = -10
      ELSE IF( lwork.LT.nw .AND. .NOT.lquery ) THEN
         info = -12
      END IF
!
      IF( info.EQ.0 ) THEN
!
!        Compute the workspace requirements
!
         nb = min( nbmax, ilaenv( 1, 'DORMQR', side // trans, m, n, &
     &             k, &
     &        -1 ) )
         lwkopt = nw*nb + tsize
         work( 1 ) = lwkopt
      END IF
!
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DORMQR', -info )
         RETURN
      ELSE IF( lquery ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( m.EQ.0 .OR. n.EQ.0 .OR. k.EQ.0 ) THEN
         work( 1 ) = 1
         RETURN
      END IF
!
      nbmin = 2
      ldwork = nw
      IF( nb.GT.1 .AND. nb.LT.k ) THEN
         IF( lwork.LT.lwkopt ) THEN
            nb = (lwork-tsize) / ldwork
            nbmin = max( 2, ilaenv( 2, 'DORMQR', side // trans, m, n, &
     &                   k, &
     &              -1 ) )
         END IF
      END IF
!
      IF( nb.LT.nbmin .OR. nb.GE.k ) THEN
!
!        Use unblocked code
!
         CALL dorm2r( side, trans, m, n, k, a, lda, tau, c, ldc, &
     &                work, &
     &                iinfo )
      ELSE
!
!        Use blocked code
!
         iwt = 1 + nw*nb
         IF( ( left .AND. .NOT.notran ) .OR. &
     &       ( .NOT.left .AND. notran ) ) THEN
            i1 = 1
            i2 = k
            i3 = nb
         ELSE
            i1 = ( ( k-1 ) / nb )*nb + 1
            i2 = 1
            i3 = -nb
         END IF
!
         IF( left ) THEN
            ni = n
            jc = 1
         ELSE
            mi = m
            ic = 1
         END IF
!
         DO 10 i = i1, i2, i3
            ib = min( nb, k-i+1 )
!
!           Form the triangular factor of the block reflector
!           H = H(i) H(i+1) . . . H(i+ib-1)
!
            CALL dlarft( 'Forward', 'Columnwise', nq-i+1, ib, a( i, &
     &                   i ), &
     &                   lda, tau( i ), work( iwt ), ldt )
            IF( left ) THEN
!
!              H or H**T is applied to C(i:m,1:n)
!
               mi = m - i + 1
               ic = i
            ELSE
!
!              H or H**T is applied to C(1:m,i:n)
!
               ni = n - i + 1
               jc = i
            END IF
!
!           Apply H or H**T
!
            CALL dlarfb( side, trans, 'Forward', 'Columnwise', mi, &
     &                   ni, &
     &                   ib, a( i, i ), lda, work( iwt ), ldt, &
     &                   c( ic, jc ), ldc, work, ldwork )
10       CONTINUE
      END IF
      work( 1 ) = lwkopt
      RETURN
!
!     End of DORMQR
!

      END
!> \brief \b DLAEXC swaps adjacent diagonal blocks of a real upper quasi-triangular matrix in Schur canonical form, by an orthogonal similarity transformation.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download DLAEXC + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaexc.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaexc.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaexc.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAEXC( WANTQ, N, T, LDT, Q, LDQ, J1, N1, N2, WORK,
!                          INFO )
!
!       .. Scalar Arguments ..
!       LOGICAL            WANTQ
!       INTEGER            INFO, J1, LDQ, LDT, N, N1, N2
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   Q( LDQ, * ), T( LDT, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAEXC swaps adjacent diagonal blocks T11 and T22 of order 1 or 2 in
!> an upper quasi-triangular matrix T by an orthogonal similarity
!> transformation.
!>
!> T must be in Schur canonical form, that is, block upper triangular
!> with 1-by-1 and 2-by-2 diagonal blocks; each 2-by-2 diagonal block
!> has its diagonal elements equal and its off-diagonal elements of
!> opposite sign.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] WANTQ
!> \verbatim
!>          WANTQ is LOGICAL
!>          = .TRUE. : accumulate the transformation in the matrix Q;
!>          = .FALSE.: do not accumulate the transformation.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix T. N >= 0.
!> \endverbatim
!>
!> \param[in,out] T
!> \verbatim
!>          T is DOUBLE PRECISION array, dimension (LDT,N)
!>          On entry, the upper quasi-triangular matrix T, in Schur
!>          canonical form.
!>          On exit, the updated matrix T, again in Schur canonical form.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of the array T. LDT >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] Q
!> \verbatim
!>          Q is DOUBLE PRECISION array, dimension (LDQ,N)
!>          On entry, if WANTQ is .TRUE., the orthogonal matrix Q.
!>          On exit, if WANTQ is .TRUE., the updated matrix Q.
!>          If WANTQ is .FALSE., Q is not referenced.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of the array Q.
!>          LDQ >= 1; and if WANTQ is .TRUE., LDQ >= N.
!> \endverbatim
!>
!> \param[in] J1
!> \verbatim
!>          J1 is INTEGER
!>          The index of the first row of the first block T11.
!> \endverbatim
!>
!> \param[in] N1
!> \verbatim
!>          N1 is INTEGER
!>          The order of the first block T11. N1 = 0, 1 or 2.
!> \endverbatim
!>
!> \param[in] N2
!> \verbatim
!>          N2 is INTEGER
!>          The order of the second block T22. N2 = 0, 1 or 2.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          = 1: the transformed matrix T would be too far from Schur
!>               form; the blocks are not swapped and T and Q are
!>               unchanged.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup laexc
!
!  =====================================================================

      SUBROUTINE dlaexc( WANTQ, N, T, LDT, Q, LDQ, J1, N1, N2, WORK, &
     &                   INFO )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      LOGICAL            WANTQ
      INTEGER            INFO, J1, LDQ, LDT, N, N1, N2
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   Q( LDQ, * ), T( LDT, * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      parameter( zero = 0.0d+0, one = 1.0d+0 )
      DOUBLE PRECISION   TEN
      parameter( ten = 1.0d+1 )
      INTEGER            LDD, LDX
      parameter( ldd = 4, ldx = 2 )
!     ..
!     .. Local Scalars ..
      INTEGER            IERR, J2, J3, J4, K, ND
      DOUBLE PRECISION   CS, DNORM, EPS, SCALE, SMLNUM, SN, T11, T22, &
     &                   t33, tau, tau1, tau2, temp, thresh, wi1, wi2, &
     &                   wr1, wr2, xnorm
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION   D( LDD, 4 ), U( 3 ), U1( 3 ), U2( 3 ), &
     &                   x( ldx, 2 )
!     ..
!     .. External Functions ..
!      DOUBLE PRECISION   DLAMCH, DLANGE
!      EXTERNAL           dlamch, dlange
!     ..
!     .. External Subroutines ..
!      EXTERNAL           dlacpy, dlanv2, dlarfg, dlarfx, dlartg, &
!     &                   dlasy2, &
!     &                   drot
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          abs, max
!     ..
!     .. Executable Statements ..
!
      info = 0
!
!     Quick return if possible
!
      IF( n.EQ.0 .OR. n1.EQ.0 .OR. n2.EQ.0 ) &
     &   RETURN
      IF( j1+n1.GT.n ) &
     &   RETURN
!
      j2 = j1 + 1
      j3 = j1 + 2
      j4 = j1 + 3
!
      IF( n1.EQ.1 .AND. n2.EQ.1 ) THEN
!
!        Swap two 1-by-1 blocks.
!
         t11 = t( j1, j1 )
         t22 = t( j2, j2 )
!
!        Determine the transformation to perform the interchange.
!
         CALL dlartg( t( j1, j2 ), t22-t11, cs, sn, temp )
!
!        Apply transformation to the matrix T.
!
         IF( j3.LE.n ) &
     &      CALL drot( n-j1-1, t( j1, j3 ), ldt, t( j2, j3 ), ldt, &
     &                 cs, &
     &                 sn )
         CALL drot( j1-1, t( 1, j1 ), 1, t( 1, j2 ), 1, cs, sn )
!
         t( j1, j1 ) = t22
         t( j2, j2 ) = t11
!
         IF( wantq ) THEN
!
!           Accumulate transformation in the matrix Q.
!
            CALL drot( n, q( 1, j1 ), 1, q( 1, j2 ), 1, cs, sn )
         END IF
!
      ELSE
!
!        Swapping involves at least one 2-by-2 block.
!
!        Copy the diagonal block of order N1+N2 to the local array D
!        and compute its norm.
!
         nd = n1 + n2
         CALL dlacpy( 'Full', nd, nd, t( j1, j1 ), ldt, d, ldd )
         dnorm = dlange( 'Max', nd, nd, d, ldd, work )
!
!        Compute machine-dependent threshold for test for accepting
!        swap.
!
         eps = dlamch( 'P' )
         smlnum = dlamch( 'S' ) / eps
         thresh = max( ten*eps*dnorm, smlnum )
!
!        Solve T11*X - X*T22 = scale*T12 for X.
!
         CALL dlasy2( .false., .false., -1, n1, n2, d, ldd, &
     &                d( n1+1, n1+1 ), ldd, d( 1, n1+1 ), ldd, scale, x, &
     &                ldx, xnorm, ierr )
!
!        Swap the adjacent diagonal blocks.
!
         k = n1 + n1 + n2 - 3
         GO TO ( 10, 20, 30 )k
!
10       CONTINUE
!
!        N1 = 1, N2 = 2: generate elementary reflector H so that:
!
!        ( scale, X11, X12 ) H = ( 0, 0, * )
!
         u( 1 ) = scale
         u( 2 ) = x( 1, 1 )
         u( 3 ) = x( 1, 2 )
         CALL dlarfg( 3, u( 3 ), u, 1, tau )
         u( 3 ) = one
         t11 = t( j1, j1 )
!
!        Perform swap provisionally on diagonal block in D.
!
         CALL dlarfx( 'L', 3, 3, u, tau, d, ldd, work )
         CALL dlarfx( 'R', 3, 3, u, tau, d, ldd, work )
!
!        Test whether to reject swap.
!
         IF( max( abs( d( 3, 1 ) ), abs( d( 3, 2 ) ), abs( d( 3, &
     &       3 )-t11 ) ).GT.thresh )GO TO 50
!
!        Accept swap: apply transformation to the entire matrix T.
!
         CALL dlarfx( 'L', 3, n-j1+1, u, tau, t( j1, j1 ), ldt, &
     &                work )
         CALL dlarfx( 'R', j2, 3, u, tau, t( 1, j1 ), ldt, work )
!
         t( j3, j1 ) = zero
         t( j3, j2 ) = zero
         t( j3, j3 ) = t11
!
         IF( wantq ) THEN
!
!           Accumulate transformation in the matrix Q.
!
            CALL dlarfx( 'R', n, 3, u, tau, q( 1, j1 ), ldq, work )
         END IF
         GO TO 40
!
20       CONTINUE
!
!        N1 = 2, N2 = 1: generate elementary reflector H so that:
!
!        H (  -X11 ) = ( * )
!          (  -X21 ) = ( 0 )
!          ( scale ) = ( 0 )
!
         u( 1 ) = -x( 1, 1 )
         u( 2 ) = -x( 2, 1 )
         u( 3 ) = scale
         CALL dlarfg( 3, u( 1 ), u( 2 ), 1, tau )
         u( 1 ) = one
         t33 = t( j3, j3 )
!
!        Perform swap provisionally on diagonal block in D.
!
         CALL dlarfx( 'L', 3, 3, u, tau, d, ldd, work )
         CALL dlarfx( 'R', 3, 3, u, tau, d, ldd, work )
!
!        Test whether to reject swap.
!
         IF( max( abs( d( 2, 1 ) ), abs( d( 3, 1 ) ), abs( d( 1, &
     &       1 )-t33 ) ).GT.thresh )GO TO 50
!
!        Accept swap: apply transformation to the entire matrix T.
!
         CALL dlarfx( 'R', j3, 3, u, tau, t( 1, j1 ), ldt, work )
         CALL dlarfx( 'L', 3, n-j1, u, tau, t( j1, j2 ), ldt, work )
!
         t( j1, j1 ) = t33
         t( j2, j1 ) = zero
         t( j3, j1 ) = zero
!
         IF( wantq ) THEN
!
!           Accumulate transformation in the matrix Q.
!
            CALL dlarfx( 'R', n, 3, u, tau, q( 1, j1 ), ldq, work )
         END IF
         GO TO 40
!
30       CONTINUE
!
!        N1 = 2, N2 = 2: generate elementary reflectors H(1) and H(2) so
!        that:
!
!        H(2) H(1) (  -X11  -X12 ) = (  *  * )
!                  (  -X21  -X22 )   (  0  * )
!                  ( scale    0  )   (  0  0 )
!                  (    0  scale )   (  0  0 )
!
         u1( 1 ) = -x( 1, 1 )
         u1( 2 ) = -x( 2, 1 )
         u1( 3 ) = scale
         CALL dlarfg( 3, u1( 1 ), u1( 2 ), 1, tau1 )
         u1( 1 ) = one
!
         temp = -tau1*( x( 1, 2 )+u1( 2 )*x( 2, 2 ) )
         u2( 1 ) = -temp*u1( 2 ) - x( 2, 2 )
         u2( 2 ) = -temp*u1( 3 )
         u2( 3 ) = scale
         CALL dlarfg( 3, u2( 1 ), u2( 2 ), 1, tau2 )
         u2( 1 ) = one
!
!        Perform swap provisionally on diagonal block in D.
!
         CALL dlarfx( 'L', 3, 4, u1, tau1, d, ldd, work )
         CALL dlarfx( 'R', 4, 3, u1, tau1, d, ldd, work )
         CALL dlarfx( 'L', 3, 4, u2, tau2, d( 2, 1 ), ldd, work )
         CALL dlarfx( 'R', 4, 3, u2, tau2, d( 1, 2 ), ldd, work )
!
!        Test whether to reject swap.
!
         IF( max( abs( d( 3, 1 ) ), abs( d( 3, 2 ) ), abs( d( 4, 1 ) ), &
     &       abs( d( 4, 2 ) ) ).GT.thresh )GO TO 50
!
!        Accept swap: apply transformation to the entire matrix T.
!
         CALL dlarfx( 'L', 3, n-j1+1, u1, tau1, t( j1, j1 ), ldt, &
     &                work )
         CALL dlarfx( 'R', j4, 3, u1, tau1, t( 1, j1 ), ldt, work )
         CALL dlarfx( 'L', 3, n-j1+1, u2, tau2, t( j2, j1 ), ldt, &
     &                work )
         CALL dlarfx( 'R', j4, 3, u2, tau2, t( 1, j2 ), ldt, work )
!
         t( j3, j1 ) = zero
         t( j3, j2 ) = zero
         t( j4, j1 ) = zero
         t( j4, j2 ) = zero
!
         IF( wantq ) THEN
!
!           Accumulate transformation in the matrix Q.
!
            CALL dlarfx( 'R', n, 3, u1, tau1, q( 1, j1 ), ldq, work )
            CALL dlarfx( 'R', n, 3, u2, tau2, q( 1, j2 ), ldq, work )
         END IF
!
40       CONTINUE
!
         IF( n2.EQ.2 ) THEN
!
!           Standardize new 2-by-2 block T11
!
            CALL dlanv2( t( j1, j1 ), t( j1, j2 ), t( j2, j1 ), &
     &                   t( j2, j2 ), wr1, wi1, wr2, wi2, cs, sn )
            CALL drot( n-j1-1, t( j1, j1+2 ), ldt, t( j2, j1+2 ), &
     &                 ldt, &
     &                 cs, sn )
            CALL drot( j1-1, t( 1, j1 ), 1, t( 1, j2 ), 1, cs, sn )
            IF( wantq ) &
     &         CALL drot( n, q( 1, j1 ), 1, q( 1, j2 ), 1, cs, sn )
         END IF
!
         IF( n1.EQ.2 ) THEN
!
!           Standardize new 2-by-2 block T22
!
            j3 = j1 + n2
            j4 = j3 + 1
            CALL dlanv2( t( j3, j3 ), t( j3, j4 ), t( j4, j3 ), &
     &                   t( j4, j4 ), wr1, wi1, wr2, wi2, cs, sn )
            IF( j3+2.LE.n ) &
     &         CALL drot( n-j3-1, t( j3, j3+2 ), ldt, t( j4, j3+2 ), &
     &                    ldt, cs, sn )
            CALL drot( j3-1, t( 1, j3 ), 1, t( 1, j4 ), 1, cs, sn )
            IF( wantq ) &
     &         CALL drot( n, q( 1, j3 ), 1, q( 1, j4 ), 1, cs, sn )
         END IF
!
      END IF
      RETURN
!
!     Exit with INFO = 1 if swap was rejected.
!
50    CONTINUE
      info = 1
      RETURN
!
!     End of DLAEXC
!

      END
!> \brief \b DLASY2 solves the Sylvester matrix equation where the matrices are of order 1 or 2.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download DLASY2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasy2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasy2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasy2.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLASY2( LTRANL, LTRANR, ISGN, N1, N2, TL, LDTL, TR,
!                          LDTR, B, LDB, SCALE, X, LDX, XNORM, INFO )
!
!       .. Scalar Arguments ..
!       LOGICAL            LTRANL, LTRANR
!       INTEGER            INFO, ISGN, LDB, LDTL, LDTR, LDX, N1, N2
!       DOUBLE PRECISION   SCALE, XNORM
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   B( LDB, * ), TL( LDTL, * ), TR( LDTR, * ),
!      $                   X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLASY2 solves for the N1 by N2 matrix X, 1 <= N1,N2 <= 2, in
!>
!>        op(TL)*X + ISGN*X*op(TR) = SCALE*B,
!>
!> where TL is N1 by N1, TR is N2 by N2, B is N1 by N2, and ISGN = 1 or
!> -1.  op(T) = T or T**T, where T**T denotes the transpose of T.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] LTRANL
!> \verbatim
!>          LTRANL is LOGICAL
!>          On entry, LTRANL specifies the op(TL):
!>             = .FALSE., op(TL) = TL,
!>             = .TRUE., op(TL) = TL**T.
!> \endverbatim
!>
!> \param[in] LTRANR
!> \verbatim
!>          LTRANR is LOGICAL
!>          On entry, LTRANR specifies the op(TR):
!>            = .FALSE., op(TR) = TR,
!>            = .TRUE., op(TR) = TR**T.
!> \endverbatim
!>
!> \param[in] ISGN
!> \verbatim
!>          ISGN is INTEGER
!>          On entry, ISGN specifies the sign of the equation
!>          as described before. ISGN may only be 1 or -1.
!> \endverbatim
!>
!> \param[in] N1
!> \verbatim
!>          N1 is INTEGER
!>          On entry, N1 specifies the order of matrix TL.
!>          N1 may only be 0, 1 or 2.
!> \endverbatim
!>
!> \param[in] N2
!> \verbatim
!>          N2 is INTEGER
!>          On entry, N2 specifies the order of matrix TR.
!>          N2 may only be 0, 1 or 2.
!> \endverbatim
!>
!> \param[in] TL
!> \verbatim
!>          TL is DOUBLE PRECISION array, dimension (LDTL,2)
!>          On entry, TL contains an N1 by N1 matrix.
!> \endverbatim
!>
!> \param[in] LDTL
!> \verbatim
!>          LDTL is INTEGER
!>          The leading dimension of the matrix TL. LDTL >= max(1,N1).
!> \endverbatim
!>
!> \param[in] TR
!> \verbatim
!>          TR is DOUBLE PRECISION array, dimension (LDTR,2)
!>          On entry, TR contains an N2 by N2 matrix.
!> \endverbatim
!>
!> \param[in] LDTR
!> \verbatim
!>          LDTR is INTEGER
!>          The leading dimension of the matrix TR. LDTR >= max(1,N2).
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (LDB,2)
!>          On entry, the N1 by N2 matrix B contains the right-hand
!>          side of the equation.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the matrix B. LDB >= max(1,N1).
!> \endverbatim
!>
!> \param[out] SCALE
!> \verbatim
!>          SCALE is DOUBLE PRECISION
!>          On exit, SCALE contains the scale factor. SCALE is chosen
!>          less than or equal to 1 to prevent the solution overflowing.
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension (LDX,2)
!>          On exit, X contains the N1 by N2 solution.
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of the matrix X. LDX >= max(1,N1).
!> \endverbatim
!>
!> \param[out] XNORM
!> \verbatim
!>          XNORM is DOUBLE PRECISION
!>          On exit, XNORM is the infinity-norm of the solution.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          On exit, INFO is set to
!>             0: successful exit.
!>             1: TL and TR have too close eigenvalues, so TL or
!>                TR is perturbed to get a nonsingular equation.
!>          NOTE: In the interests of speed, this routine does not
!>                check the inputs for errors.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup lasy2
!
!  =====================================================================

      SUBROUTINE dlasy2( LTRANL, LTRANR, ISGN, N1, N2, TL, LDTL, TR, &
     &                   LDTR, B, LDB, SCALE, X, LDX, XNORM, INFO )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      LOGICAL            LTRANL, LTRANR
      INTEGER            INFO, ISGN, LDB, LDTL, LDTR, LDX, N1, N2
      DOUBLE PRECISION   SCALE, XNORM
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   B( LDB, * ), TL( LDTL, * ), TR( LDTR, * ), &
     &                   x( ldx, * )
!     ..
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      parameter( zero = 0.0d+0, one = 1.0d+0 )
      DOUBLE PRECISION   TWO, HALF, EIGHT
      parameter( two = 2.0d+0, half = 0.5d+0, eight = 8.0d+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            BSWAP, XSWAP
      INTEGER            I, IP, IPIV, IPSV, J, JP, JPSV, K
      DOUBLE PRECISION   BET, EPS, GAM, L21, SGN, SMIN, SMLNUM, TAU1, &
     &                   temp, u11, u12, u22, xmax
!     ..
!     .. Local Arrays ..
      LOGICAL            BSWPIV( 4 ), XSWPIV( 4 )
      INTEGER            JPIV( 4 ), LOCL21( 4 ), LOCU12( 4 ), &
     &                   locu22( 4 )
      DOUBLE PRECISION   BTMP( 4 ), T16( 4, 4 ), TMP( 4 ), X2( 2 )
!     ..
!     .. External Functions ..
!      INTEGER            IDAMAX
!      DOUBLE PRECISION   DLAMCH
!      EXTERNAL           idamax, dlamch
!     ..
!     .. External Subroutines ..
!      EXTERNAL           dcopy, dswap
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          abs, max
!     ..
!     .. Data statements ..
      DATA               locu12 / 3, 4, 1, 2 / , locl21 / 2, 1, 4, 3 / , &
     &                   locu22 / 4, 3, 2, 1 /
      DATA               xswpiv / .false., .false., .true., .true. /
      DATA               bswpiv / .false., .true., .false., .true. /
!     ..
!     .. Executable Statements ..
!
!     Do not check the input parameters for errors
!
      info = 0
!
!     Quick return if possible
!
      IF( n1.EQ.0 .OR. n2.EQ.0 ) &
     &   RETURN
!
!     Set constants to control overflow
!
      eps = dlamch( 'P' )
      smlnum = dlamch( 'S' ) / eps
      sgn = isgn
!
      k = n1 + n1 + n2 - 2
      GO TO ( 10, 20, 30, 50 )k
!
!     1 by 1: TL11*X + SGN*X*TR11 = B11
!
10    CONTINUE
      tau1 = tl( 1, 1 ) + sgn*tr( 1, 1 )
      bet = abs( tau1 )
      IF( bet.LE.smlnum ) THEN
         tau1 = smlnum
         bet = smlnum
         info = 1
      END IF
!
      scale = one
      gam = abs( b( 1, 1 ) )
      IF( smlnum*gam.GT.bet ) &
     &   scale = one / gam
!
      x( 1, 1 ) = ( b( 1, 1 )*scale ) / tau1
      xnorm = abs( x( 1, 1 ) )
      RETURN
!
!     1 by 2:
!     TL11*[X11 X12] + ISGN*[X11 X12]*op[TR11 TR12]  = [B11 B12]
!                                       [TR21 TR22]
!
20    CONTINUE
!
      smin = max( eps*max( abs( tl( 1, 1 ) ), abs( tr( 1, 1 ) ), &
     &       abs( tr( 1, 2 ) ), abs( tr( 2, 1 ) ), abs( tr( 2, 2 ) ) ), &
     &       smlnum )
      tmp( 1 ) = tl( 1, 1 ) + sgn*tr( 1, 1 )
      tmp( 4 ) = tl( 1, 1 ) + sgn*tr( 2, 2 )
      IF( ltranr ) THEN
         tmp( 2 ) = sgn*tr( 2, 1 )
         tmp( 3 ) = sgn*tr( 1, 2 )
      ELSE
         tmp( 2 ) = sgn*tr( 1, 2 )
         tmp( 3 ) = sgn*tr( 2, 1 )
      END IF
      btmp( 1 ) = b( 1, 1 )
      btmp( 2 ) = b( 1, 2 )
      GO TO 40
!
!     2 by 1:
!          op[TL11 TL12]*[X11] + ISGN* [X11]*TR11  = [B11]
!            [TL21 TL22] [X21]         [X21]         [B21]
!
30    CONTINUE
      smin = max( eps*max( abs( tr( 1, 1 ) ), abs( tl( 1, 1 ) ), &
     &       abs( tl( 1, 2 ) ), abs( tl( 2, 1 ) ), abs( tl( 2, 2 ) ) ), &
     &       smlnum )
      tmp( 1 ) = tl( 1, 1 ) + sgn*tr( 1, 1 )
      tmp( 4 ) = tl( 2, 2 ) + sgn*tr( 1, 1 )
      IF( ltranl ) THEN
         tmp( 2 ) = tl( 1, 2 )
         tmp( 3 ) = tl( 2, 1 )
      ELSE
         tmp( 2 ) = tl( 2, 1 )
         tmp( 3 ) = tl( 1, 2 )
      END IF
      btmp( 1 ) = b( 1, 1 )
      btmp( 2 ) = b( 2, 1 )
40    CONTINUE
!
!     Solve 2 by 2 system using complete pivoting.
!     Set pivots less than SMIN to SMIN.
!
      ipiv = idamax( 4, tmp, 1 )
      u11 = tmp( ipiv )
      IF( abs( u11 ).LE.smin ) THEN
         info = 1
         u11 = smin
      END IF
      u12 = tmp( locu12( ipiv ) )
      l21 = tmp( locl21( ipiv ) ) / u11
      u22 = tmp( locu22( ipiv ) ) - u12*l21
      xswap = xswpiv( ipiv )
      bswap = bswpiv( ipiv )
      IF( abs( u22 ).LE.smin ) THEN
         info = 1
         u22 = smin
      END IF
      IF( bswap ) THEN
         temp = btmp( 2 )
         btmp( 2 ) = btmp( 1 ) - l21*temp
         btmp( 1 ) = temp
      ELSE
         btmp( 2 ) = btmp( 2 ) - l21*btmp( 1 )
      END IF
      scale = one
      IF( ( two*smlnum )*abs( btmp( 2 ) ).GT.abs( u22 ) .OR. &
     &    ( two*smlnum )*abs( btmp( 1 ) ).GT.abs( u11 ) ) THEN
         scale = half / max( abs( btmp( 1 ) ), abs( btmp( 2 ) ) )
         btmp( 1 ) = btmp( 1 )*scale
         btmp( 2 ) = btmp( 2 )*scale
      END IF
      x2( 2 ) = btmp( 2 ) / u22
      x2( 1 ) = btmp( 1 ) / u11 - ( u12 / u11 )*x2( 2 )
      IF( xswap ) THEN
         temp = x2( 2 )
         x2( 2 ) = x2( 1 )
         x2( 1 ) = temp
      END IF
      x( 1, 1 ) = x2( 1 )
      IF( n1.EQ.1 ) THEN
         x( 1, 2 ) = x2( 2 )
         xnorm = abs( x( 1, 1 ) ) + abs( x( 1, 2 ) )
      ELSE
         x( 2, 1 ) = x2( 2 )
         xnorm = max( abs( x( 1, 1 ) ), abs( x( 2, 1 ) ) )
      END IF
      RETURN
!
!     2 by 2:
!     op[TL11 TL12]*[X11 X12] +ISGN* [X11 X12]*op[TR11 TR12] = [B11 B12]
!       [TL21 TL22] [X21 X22]        [X21 X22]   [TR21 TR22]   [B21 B22]
!
!     Solve equivalent 4 by 4 system using complete pivoting.
!     Set pivots less than SMIN to SMIN.
!
50    CONTINUE
      smin = max( abs( tr( 1, 1 ) ), abs( tr( 1, 2 ) ), &
     &       abs( tr( 2, 1 ) ), abs( tr( 2, 2 ) ) )
      smin = max( smin, abs( tl( 1, 1 ) ), abs( tl( 1, 2 ) ), &
     &       abs( tl( 2, 1 ) ), abs( tl( 2, 2 ) ) )
      smin = max( eps*smin, smlnum )
      btmp( 1 ) = zero
      CALL dcopy( 16, btmp, 0, t16, 1 )
      t16( 1, 1 ) = tl( 1, 1 ) + sgn*tr( 1, 1 )
      t16( 2, 2 ) = tl( 2, 2 ) + sgn*tr( 1, 1 )
      t16( 3, 3 ) = tl( 1, 1 ) + sgn*tr( 2, 2 )
      t16( 4, 4 ) = tl( 2, 2 ) + sgn*tr( 2, 2 )
      IF( ltranl ) THEN
         t16( 1, 2 ) = tl( 2, 1 )
         t16( 2, 1 ) = tl( 1, 2 )
         t16( 3, 4 ) = tl( 2, 1 )
         t16( 4, 3 ) = tl( 1, 2 )
      ELSE
         t16( 1, 2 ) = tl( 1, 2 )
         t16( 2, 1 ) = tl( 2, 1 )
         t16( 3, 4 ) = tl( 1, 2 )
         t16( 4, 3 ) = tl( 2, 1 )
      END IF
      IF( ltranr ) THEN
         t16( 1, 3 ) = sgn*tr( 1, 2 )
         t16( 2, 4 ) = sgn*tr( 1, 2 )
         t16( 3, 1 ) = sgn*tr( 2, 1 )
         t16( 4, 2 ) = sgn*tr( 2, 1 )
      ELSE
         t16( 1, 3 ) = sgn*tr( 2, 1 )
         t16( 2, 4 ) = sgn*tr( 2, 1 )
         t16( 3, 1 ) = sgn*tr( 1, 2 )
         t16( 4, 2 ) = sgn*tr( 1, 2 )
      END IF
      btmp( 1 ) = b( 1, 1 )
      btmp( 2 ) = b( 2, 1 )
      btmp( 3 ) = b( 1, 2 )
      btmp( 4 ) = b( 2, 2 )
!
!     Perform elimination
!
      DO 100 i = 1, 3
         xmax = zero
         DO 70 ip = i, 4
            DO 60 jp = i, 4
               IF( abs( t16( ip, jp ) ).GE.xmax ) THEN
                  xmax = abs( t16( ip, jp ) )
                  ipsv = ip
                  jpsv = jp
               END IF
60          CONTINUE
70       CONTINUE
         IF( ipsv.NE.i ) THEN
            CALL dswap( 4, t16( ipsv, 1 ), 4, t16( i, 1 ), 4 )
            temp = btmp( i )
            btmp( i ) = btmp( ipsv )
            btmp( ipsv ) = temp
         END IF
         IF( jpsv.NE.i ) &
     &      CALL dswap( 4, t16( 1, jpsv ), 1, t16( 1, i ), 1 )
         jpiv( i ) = jpsv
         IF( abs( t16( i, i ) ).LT.smin ) THEN
            info = 1
            t16( i, i ) = smin
         END IF
         DO 90 j = i + 1, 4
            t16( j, i ) = t16( j, i ) / t16( i, i )
            btmp( j ) = btmp( j ) - t16( j, i )*btmp( i )
            DO 80 k = i + 1, 4
               t16( j, k ) = t16( j, k ) - t16( j, i )*t16( i, k )
80          CONTINUE
90       CONTINUE
100   CONTINUE
      IF( abs( t16( 4, 4 ) ).LT.smin ) THEN
         info = 1
         t16( 4, 4 ) = smin
      END IF
      scale = one
      IF( ( eight*smlnum )*abs( btmp( 1 ) ).GT.abs( t16( 1, 1 ) ) .OR. &
     &    ( eight*smlnum )*abs( btmp( 2 ) ).GT.abs( t16( 2, 2 ) ) .OR. &
     &    ( eight*smlnum )*abs( btmp( 3 ) ).GT.abs( t16( 3, 3 ) ) .OR. &
     &    ( eight*smlnum )*abs( btmp( 4 ) ).GT.abs( t16( 4, 4 ) ) ) THEN
         scale = ( one / eight ) / max( abs( btmp( 1 ) ), &
     &           abs( btmp( 2 ) ), abs( btmp( 3 ) ), abs( btmp( 4 ) ) )
         btmp( 1 ) = btmp( 1 )*scale
         btmp( 2 ) = btmp( 2 )*scale
         btmp( 3 ) = btmp( 3 )*scale
         btmp( 4 ) = btmp( 4 )*scale
      END IF
      DO 120 i = 1, 4
         k = 5 - i
         temp = one / t16( k, k )
         tmp( k ) = btmp( k )*temp
         DO 110 j = k + 1, 4
            tmp( k ) = tmp( k ) - ( temp*t16( k, j ) )*tmp( j )
110      CONTINUE
120   CONTINUE
      DO 130 i = 1, 3
         IF( jpiv( 4-i ).NE.4-i ) THEN
            temp = tmp( 4-i )
            tmp( 4-i ) = tmp( jpiv( 4-i ) )
            tmp( jpiv( 4-i ) ) = temp
         END IF
130   CONTINUE
      x( 1, 1 ) = tmp( 1 )
      x( 2, 1 ) = tmp( 2 )
      x( 1, 2 ) = tmp( 3 )
      x( 2, 2 ) = tmp( 4 )
      xnorm = max( abs( tmp( 1 ) )+abs( tmp( 3 ) ), &
     &        abs( tmp( 2 ) )+abs( tmp( 4 ) ) )
      RETURN
!
!     End of DLASY2
!

      END
!> \brief \b DLARFX applies an elementary reflector to a general rectangular matrix, with loop unrolling when the reflector has order ≤ 10.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download DLARFX + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarfx.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarfx.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarfx.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLARFX( SIDE, M, N, V, TAU, C, LDC, WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          SIDE
!       INTEGER            LDC, M, N
!       DOUBLE PRECISION   TAU
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   C( LDC, * ), V( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLARFX applies a real elementary reflector H to a real m by n
!> matrix C, from either the left or the right. H is represented in the
!> form
!>
!>       H = I - tau * v * v**T
!>
!> where tau is a real scalar and v is a real vector.
!>
!> If tau = 0, then H is taken to be the unit matrix
!>
!> This version uses inline code if H has order < 11.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'L': form  H * C
!>          = 'R': form  C * H
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix C.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix C.
!> \endverbatim
!>
!> \param[in] V
!> \verbatim
!>          V is DOUBLE PRECISION array, dimension (M) if SIDE = 'L'
!>                                     or (N) if SIDE = 'R'
!>          The vector v in the representation of H.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION
!>          The value tau in the representation of H.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension (LDC,N)
!>          On entry, the m by n matrix C.
!>          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
!>          or C * H if SIDE = 'R'.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C. LDC >= (1,M).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension
!>                      (N) if SIDE = 'L'
!>                      or (M) if SIDE = 'R'
!>          WORK is not referenced if H has order < 11.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup larfx
!
!  =====================================================================

      SUBROUTINE dlarfx( SIDE, M, N, V, TAU, C, LDC, WORK )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            LDC, M, N
      DOUBLE PRECISION   TAU
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   C( LDC, * ), V( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      parameter( zero = 0.0d+0, one = 1.0d+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            J
      DOUBLE PRECISION   SUM, T1, T10, T2, T3, T4, T5, T6, T7, T8, T9, &
     &                   V1, V10, V2, V3, V4, V5, V6, V7, V8, V9
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      EXTERNAL           lsame
!     ..
!     .. External Subroutines ..
!      EXTERNAL           dlarf
!     ..
!     .. Executable Statements ..
!
      IF( tau.EQ.zero ) &
     &   RETURN
      IF( lsame( side, 'L' ) ) THEN
!
!        Form  H * C, where H has order m.
!
         GO TO ( 10, 30, 50, 70, 90, 110, 130, 150, &
     &           170, 190 )m
!
!        Code for general M
!
         CALL dlarf( side, m, n, v, 1, tau, c, ldc, work )
         GO TO 410
10       CONTINUE
!
!        Special code for 1 x 1 Householder
!
         t1 = one - tau*v( 1 )*v( 1 )
         DO 20 j = 1, n
            c( 1, j ) = t1*c( 1, j )
20       CONTINUE
         GO TO 410
30       CONTINUE
!
!        Special code for 2 x 2 Householder
!
         v1 = v( 1 )
         t1 = tau*v1
         v2 = v( 2 )
         t2 = tau*v2
         DO 40 j = 1, n
            sum = v1*c( 1, j ) + v2*c( 2, j )
            c( 1, j ) = c( 1, j ) - sum*t1
            c( 2, j ) = c( 2, j ) - sum*t2
40       CONTINUE
         GO TO 410
50       CONTINUE
!
!        Special code for 3 x 3 Householder
!
         v1 = v( 1 )
         t1 = tau*v1
         v2 = v( 2 )
         t2 = tau*v2
         v3 = v( 3 )
         t3 = tau*v3
         DO 60 j = 1, n
            sum = v1*c( 1, j ) + v2*c( 2, j ) + v3*c( 3, j )
            c( 1, j ) = c( 1, j ) - sum*t1
            c( 2, j ) = c( 2, j ) - sum*t2
            c( 3, j ) = c( 3, j ) - sum*t3
60       CONTINUE
         GO TO 410
70       CONTINUE
!
!        Special code for 4 x 4 Householder
!
         v1 = v( 1 )
         t1 = tau*v1
         v2 = v( 2 )
         t2 = tau*v2
         v3 = v( 3 )
         t3 = tau*v3
         v4 = v( 4 )
         t4 = tau*v4
         DO 80 j = 1, n
            sum = v1*c( 1, j ) + v2*c( 2, j ) + v3*c( 3, j ) + &
     &            v4*c( 4, j )
            c( 1, j ) = c( 1, j ) - sum*t1
            c( 2, j ) = c( 2, j ) - sum*t2
            c( 3, j ) = c( 3, j ) - sum*t3
            c( 4, j ) = c( 4, j ) - sum*t4
80       CONTINUE
         GO TO 410
90       CONTINUE
!
!        Special code for 5 x 5 Householder
!
         v1 = v( 1 )
         t1 = tau*v1
         v2 = v( 2 )
         t2 = tau*v2
         v3 = v( 3 )
         t3 = tau*v3
         v4 = v( 4 )
         t4 = tau*v4
         v5 = v( 5 )
         t5 = tau*v5
         DO 100 j = 1, n
            sum = v1*c( 1, j ) + v2*c( 2, j ) + v3*c( 3, j ) + &
     &            v4*c( 4, j ) + v5*c( 5, j )
            c( 1, j ) = c( 1, j ) - sum*t1
            c( 2, j ) = c( 2, j ) - sum*t2
            c( 3, j ) = c( 3, j ) - sum*t3
            c( 4, j ) = c( 4, j ) - sum*t4
            c( 5, j ) = c( 5, j ) - sum*t5
100      CONTINUE
         GO TO 410
110      CONTINUE
!
!        Special code for 6 x 6 Householder
!
         v1 = v( 1 )
         t1 = tau*v1
         v2 = v( 2 )
         t2 = tau*v2
         v3 = v( 3 )
         t3 = tau*v3
         v4 = v( 4 )
         t4 = tau*v4
         v5 = v( 5 )
         t5 = tau*v5
         v6 = v( 6 )
         t6 = tau*v6
         DO 120 j = 1, n
            sum = v1*c( 1, j ) + v2*c( 2, j ) + v3*c( 3, j ) + &
     &            v4*c( 4, j ) + v5*c( 5, j ) + v6*c( 6, j )
            c( 1, j ) = c( 1, j ) - sum*t1
            c( 2, j ) = c( 2, j ) - sum*t2
            c( 3, j ) = c( 3, j ) - sum*t3
            c( 4, j ) = c( 4, j ) - sum*t4
            c( 5, j ) = c( 5, j ) - sum*t5
            c( 6, j ) = c( 6, j ) - sum*t6
120      CONTINUE
         GO TO 410
130      CONTINUE
!
!        Special code for 7 x 7 Householder
!
         v1 = v( 1 )
         t1 = tau*v1
         v2 = v( 2 )
         t2 = tau*v2
         v3 = v( 3 )
         t3 = tau*v3
         v4 = v( 4 )
         t4 = tau*v4
         v5 = v( 5 )
         t5 = tau*v5
         v6 = v( 6 )
         t6 = tau*v6
         v7 = v( 7 )
         t7 = tau*v7
         DO 140 j = 1, n
            sum = v1*c( 1, j ) + v2*c( 2, j ) + v3*c( 3, j ) + &
     &            v4*c( 4, j ) + v5*c( 5, j ) + v6*c( 6, j ) + &
     &            v7*c( 7, j )
            c( 1, j ) = c( 1, j ) - sum*t1
            c( 2, j ) = c( 2, j ) - sum*t2
            c( 3, j ) = c( 3, j ) - sum*t3
            c( 4, j ) = c( 4, j ) - sum*t4
            c( 5, j ) = c( 5, j ) - sum*t5
            c( 6, j ) = c( 6, j ) - sum*t6
            c( 7, j ) = c( 7, j ) - sum*t7
140      CONTINUE
         GO TO 410
150      CONTINUE
!
!        Special code for 8 x 8 Householder
!
         v1 = v( 1 )
         t1 = tau*v1
         v2 = v( 2 )
         t2 = tau*v2
         v3 = v( 3 )
         t3 = tau*v3
         v4 = v( 4 )
         t4 = tau*v4
         v5 = v( 5 )
         t5 = tau*v5
         v6 = v( 6 )
         t6 = tau*v6
         v7 = v( 7 )
         t7 = tau*v7
         v8 = v( 8 )
         t8 = tau*v8
         DO 160 j = 1, n
            sum = v1*c( 1, j ) + v2*c( 2, j ) + v3*c( 3, j ) + &
     &            v4*c( 4, j ) + v5*c( 5, j ) + v6*c( 6, j ) + &
     &            v7*c( 7, j ) + v8*c( 8, j )
            c( 1, j ) = c( 1, j ) - sum*t1
            c( 2, j ) = c( 2, j ) - sum*t2
            c( 3, j ) = c( 3, j ) - sum*t3
            c( 4, j ) = c( 4, j ) - sum*t4
            c( 5, j ) = c( 5, j ) - sum*t5
            c( 6, j ) = c( 6, j ) - sum*t6
            c( 7, j ) = c( 7, j ) - sum*t7
            c( 8, j ) = c( 8, j ) - sum*t8
160      CONTINUE
         GO TO 410
170      CONTINUE
!
!        Special code for 9 x 9 Householder
!
         v1 = v( 1 )
         t1 = tau*v1
         v2 = v( 2 )
         t2 = tau*v2
         v3 = v( 3 )
         t3 = tau*v3
         v4 = v( 4 )
         t4 = tau*v4
         v5 = v( 5 )
         t5 = tau*v5
         v6 = v( 6 )
         t6 = tau*v6
         v7 = v( 7 )
         t7 = tau*v7
         v8 = v( 8 )
         t8 = tau*v8
         v9 = v( 9 )
         t9 = tau*v9
         DO 180 j = 1, n
            sum = v1*c( 1, j ) + v2*c( 2, j ) + v3*c( 3, j ) + &
     &            v4*c( 4, j ) + v5*c( 5, j ) + v6*c( 6, j ) + &
     &            v7*c( 7, j ) + v8*c( 8, j ) + v9*c( 9, j )
            c( 1, j ) = c( 1, j ) - sum*t1
            c( 2, j ) = c( 2, j ) - sum*t2
            c( 3, j ) = c( 3, j ) - sum*t3
            c( 4, j ) = c( 4, j ) - sum*t4
            c( 5, j ) = c( 5, j ) - sum*t5
            c( 6, j ) = c( 6, j ) - sum*t6
            c( 7, j ) = c( 7, j ) - sum*t7
            c( 8, j ) = c( 8, j ) - sum*t8
            c( 9, j ) = c( 9, j ) - sum*t9
180      CONTINUE
         GO TO 410
190      CONTINUE
!
!        Special code for 10 x 10 Householder
!
         v1 = v( 1 )
         t1 = tau*v1
         v2 = v( 2 )
         t2 = tau*v2
         v3 = v( 3 )
         t3 = tau*v3
         v4 = v( 4 )
         t4 = tau*v4
         v5 = v( 5 )
         t5 = tau*v5
         v6 = v( 6 )
         t6 = tau*v6
         v7 = v( 7 )
         t7 = tau*v7
         v8 = v( 8 )
         t8 = tau*v8
         v9 = v( 9 )
         t9 = tau*v9
         v10 = v( 10 )
         t10 = tau*v10
         DO 200 j = 1, n
            sum = v1*c( 1, j ) + v2*c( 2, j ) + v3*c( 3, j ) + &
     &            v4*c( 4, j ) + v5*c( 5, j ) + v6*c( 6, j ) + &
     &            v7*c( 7, j ) + v8*c( 8, j ) + v9*c( 9, j ) + &
     &            v10*c( 10, j )
            c( 1, j ) = c( 1, j ) - sum*t1
            c( 2, j ) = c( 2, j ) - sum*t2
            c( 3, j ) = c( 3, j ) - sum*t3
            c( 4, j ) = c( 4, j ) - sum*t4
            c( 5, j ) = c( 5, j ) - sum*t5
            c( 6, j ) = c( 6, j ) - sum*t6
            c( 7, j ) = c( 7, j ) - sum*t7
            c( 8, j ) = c( 8, j ) - sum*t8
            c( 9, j ) = c( 9, j ) - sum*t9
            c( 10, j ) = c( 10, j ) - sum*t10
200      CONTINUE
         GO TO 410
      ELSE
!
!        Form  C * H, where H has order n.
!
         GO TO ( 210, 230, 250, 270, 290, 310, 330, 350, &
     &           370, 390 )n
!
!        Code for general N
!
         CALL dlarf( side, m, n, v, 1, tau, c, ldc, work )
         GO TO 410
210      CONTINUE
!
!        Special code for 1 x 1 Householder
!
         t1 = one - tau*v( 1 )*v( 1 )
         DO 220 j = 1, m
            c( j, 1 ) = t1*c( j, 1 )
220      CONTINUE
         GO TO 410
230      CONTINUE
!
!        Special code for 2 x 2 Householder
!
         v1 = v( 1 )
         t1 = tau*v1
         v2 = v( 2 )
         t2 = tau*v2
         DO 240 j = 1, m
            sum = v1*c( j, 1 ) + v2*c( j, 2 )
            c( j, 1 ) = c( j, 1 ) - sum*t1
            c( j, 2 ) = c( j, 2 ) - sum*t2
240      CONTINUE
         GO TO 410
250      CONTINUE
!
!        Special code for 3 x 3 Householder
!
         v1 = v( 1 )
         t1 = tau*v1
         v2 = v( 2 )
         t2 = tau*v2
         v3 = v( 3 )
         t3 = tau*v3
         DO 260 j = 1, m
            sum = v1*c( j, 1 ) + v2*c( j, 2 ) + v3*c( j, 3 )
            c( j, 1 ) = c( j, 1 ) - sum*t1
            c( j, 2 ) = c( j, 2 ) - sum*t2
            c( j, 3 ) = c( j, 3 ) - sum*t3
260      CONTINUE
         GO TO 410
270      CONTINUE
!
!        Special code for 4 x 4 Householder
!
         v1 = v( 1 )
         t1 = tau*v1
         v2 = v( 2 )
         t2 = tau*v2
         v3 = v( 3 )
         t3 = tau*v3
         v4 = v( 4 )
         t4 = tau*v4
         DO 280 j = 1, m
            sum = v1*c( j, 1 ) + v2*c( j, 2 ) + v3*c( j, 3 ) + &
     &            v4*c( j, 4 )
            c( j, 1 ) = c( j, 1 ) - sum*t1
            c( j, 2 ) = c( j, 2 ) - sum*t2
            c( j, 3 ) = c( j, 3 ) - sum*t3
            c( j, 4 ) = c( j, 4 ) - sum*t4
280      CONTINUE
         GO TO 410
290      CONTINUE
!
!        Special code for 5 x 5 Householder
!
         v1 = v( 1 )
         t1 = tau*v1
         v2 = v( 2 )
         t2 = tau*v2
         v3 = v( 3 )
         t3 = tau*v3
         v4 = v( 4 )
         t4 = tau*v4
         v5 = v( 5 )
         t5 = tau*v5
         DO 300 j = 1, m
            sum = v1*c( j, 1 ) + v2*c( j, 2 ) + v3*c( j, 3 ) + &
     &            v4*c( j, 4 ) + v5*c( j, 5 )
            c( j, 1 ) = c( j, 1 ) - sum*t1
            c( j, 2 ) = c( j, 2 ) - sum*t2
            c( j, 3 ) = c( j, 3 ) - sum*t3
            c( j, 4 ) = c( j, 4 ) - sum*t4
            c( j, 5 ) = c( j, 5 ) - sum*t5
300      CONTINUE
         GO TO 410
310      CONTINUE
!
!        Special code for 6 x 6 Householder
!
         v1 = v( 1 )
         t1 = tau*v1
         v2 = v( 2 )
         t2 = tau*v2
         v3 = v( 3 )
         t3 = tau*v3
         v4 = v( 4 )
         t4 = tau*v4
         v5 = v( 5 )
         t5 = tau*v5
         v6 = v( 6 )
         t6 = tau*v6
         DO 320 j = 1, m
            sum = v1*c( j, 1 ) + v2*c( j, 2 ) + v3*c( j, 3 ) + &
     &            v4*c( j, 4 ) + v5*c( j, 5 ) + v6*c( j, 6 )
            c( j, 1 ) = c( j, 1 ) - sum*t1
            c( j, 2 ) = c( j, 2 ) - sum*t2
            c( j, 3 ) = c( j, 3 ) - sum*t3
            c( j, 4 ) = c( j, 4 ) - sum*t4
            c( j, 5 ) = c( j, 5 ) - sum*t5
            c( j, 6 ) = c( j, 6 ) - sum*t6
320      CONTINUE
         GO TO 410
330      CONTINUE
!
!        Special code for 7 x 7 Householder
!
         v1 = v( 1 )
         t1 = tau*v1
         v2 = v( 2 )
         t2 = tau*v2
         v3 = v( 3 )
         t3 = tau*v3
         v4 = v( 4 )
         t4 = tau*v4
         v5 = v( 5 )
         t5 = tau*v5
         v6 = v( 6 )
         t6 = tau*v6
         v7 = v( 7 )
         t7 = tau*v7
         DO 340 j = 1, m
            sum = v1*c( j, 1 ) + v2*c( j, 2 ) + v3*c( j, 3 ) + &
     &            v4*c( j, 4 ) + v5*c( j, 5 ) + v6*c( j, 6 ) + &
     &            v7*c( j, 7 )
            c( j, 1 ) = c( j, 1 ) - sum*t1
            c( j, 2 ) = c( j, 2 ) - sum*t2
            c( j, 3 ) = c( j, 3 ) - sum*t3
            c( j, 4 ) = c( j, 4 ) - sum*t4
            c( j, 5 ) = c( j, 5 ) - sum*t5
            c( j, 6 ) = c( j, 6 ) - sum*t6
            c( j, 7 ) = c( j, 7 ) - sum*t7
340      CONTINUE
         GO TO 410
350      CONTINUE
!
!        Special code for 8 x 8 Householder
!
         v1 = v( 1 )
         t1 = tau*v1
         v2 = v( 2 )
         t2 = tau*v2
         v3 = v( 3 )
         t3 = tau*v3
         v4 = v( 4 )
         t4 = tau*v4
         v5 = v( 5 )
         t5 = tau*v5
         v6 = v( 6 )
         t6 = tau*v6
         v7 = v( 7 )
         t7 = tau*v7
         v8 = v( 8 )
         t8 = tau*v8
         DO 360 j = 1, m
            sum = v1*c( j, 1 ) + v2*c( j, 2 ) + v3*c( j, 3 ) + &
     &            v4*c( j, 4 ) + v5*c( j, 5 ) + v6*c( j, 6 ) + &
     &            v7*c( j, 7 ) + v8*c( j, 8 )
            c( j, 1 ) = c( j, 1 ) - sum*t1
            c( j, 2 ) = c( j, 2 ) - sum*t2
            c( j, 3 ) = c( j, 3 ) - sum*t3
            c( j, 4 ) = c( j, 4 ) - sum*t4
            c( j, 5 ) = c( j, 5 ) - sum*t5
            c( j, 6 ) = c( j, 6 ) - sum*t6
            c( j, 7 ) = c( j, 7 ) - sum*t7
            c( j, 8 ) = c( j, 8 ) - sum*t8
360      CONTINUE
         GO TO 410
370      CONTINUE
!
!        Special code for 9 x 9 Householder
!
         v1 = v( 1 )
         t1 = tau*v1
         v2 = v( 2 )
         t2 = tau*v2
         v3 = v( 3 )
         t3 = tau*v3
         v4 = v( 4 )
         t4 = tau*v4
         v5 = v( 5 )
         t5 = tau*v5
         v6 = v( 6 )
         t6 = tau*v6
         v7 = v( 7 )
         t7 = tau*v7
         v8 = v( 8 )
         t8 = tau*v8
         v9 = v( 9 )
         t9 = tau*v9
         DO 380 j = 1, m
            sum = v1*c( j, 1 ) + v2*c( j, 2 ) + v3*c( j, 3 ) + &
     &            v4*c( j, 4 ) + v5*c( j, 5 ) + v6*c( j, 6 ) + &
     &            v7*c( j, 7 ) + v8*c( j, 8 ) + v9*c( j, 9 )
            c( j, 1 ) = c( j, 1 ) - sum*t1
            c( j, 2 ) = c( j, 2 ) - sum*t2
            c( j, 3 ) = c( j, 3 ) - sum*t3
            c( j, 4 ) = c( j, 4 ) - sum*t4
            c( j, 5 ) = c( j, 5 ) - sum*t5
            c( j, 6 ) = c( j, 6 ) - sum*t6
            c( j, 7 ) = c( j, 7 ) - sum*t7
            c( j, 8 ) = c( j, 8 ) - sum*t8
            c( j, 9 ) = c( j, 9 ) - sum*t9
380      CONTINUE
         GO TO 410
390      CONTINUE
!
!        Special code for 10 x 10 Householder
!
         v1 = v( 1 )
         t1 = tau*v1
         v2 = v( 2 )
         t2 = tau*v2
         v3 = v( 3 )
         t3 = tau*v3
         v4 = v( 4 )
         t4 = tau*v4
         v5 = v( 5 )
         t5 = tau*v5
         v6 = v( 6 )
         t6 = tau*v6
         v7 = v( 7 )
         t7 = tau*v7
         v8 = v( 8 )
         t8 = tau*v8
         v9 = v( 9 )
         t9 = tau*v9
         v10 = v( 10 )
         t10 = tau*v10
         DO 400 j = 1, m
            sum = v1*c( j, 1 ) + v2*c( j, 2 ) + v3*c( j, 3 ) + &
     &            v4*c( j, 4 ) + v5*c( j, 5 ) + v6*c( j, 6 ) + &
     &            v7*c( j, 7 ) + v8*c( j, 8 ) + v9*c( j, 9 ) + &
     &            v10*c( j, 10 )
            c( j, 1 ) = c( j, 1 ) - sum*t1
            c( j, 2 ) = c( j, 2 ) - sum*t2
            c( j, 3 ) = c( j, 3 ) - sum*t3
            c( j, 4 ) = c( j, 4 ) - sum*t4
            c( j, 5 ) = c( j, 5 ) - sum*t5
            c( j, 6 ) = c( j, 6 ) - sum*t6
            c( j, 7 ) = c( j, 7 ) - sum*t7
            c( j, 8 ) = c( j, 8 ) - sum*t8
            c( j, 9 ) = c( j, 9 ) - sum*t9
            c( j, 10 ) = c( j, 10 ) - sum*t10
400      CONTINUE
         GO TO 410
      END IF
410   CONTINUE
      RETURN
!
!     End of DLARFX
!

      END
!> \brief \b DORM2R multiplies a general matrix by the orthogonal matrix from a QR factorization determined by sgeqrf (unblocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> Download DORM2R + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorm2r.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorm2r.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorm2r.f">
!> [TXT]</a>
!
!  Definition:
!  ===========
!
!       SUBROUTINE DORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
!                          WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          SIDE, TRANS
!       INTEGER            INFO, K, LDA, LDC, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DORM2R overwrites the general real m by n matrix C with
!>
!>       Q * C  if SIDE = 'L' and TRANS = 'N', or
!>
!>       Q**T* C  if SIDE = 'L' and TRANS = 'T', or
!>
!>       C * Q  if SIDE = 'R' and TRANS = 'N', or
!>
!>       C * Q**T if SIDE = 'R' and TRANS = 'T',
!>
!> where Q is a real orthogonal matrix defined as the product of k
!> elementary reflectors
!>
!>       Q = H(1) H(2) . . . H(k)
!>
!> as returned by DGEQRF. Q is of order m if SIDE = 'L' and of order n
!> if SIDE = 'R'.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'L': apply Q or Q**T from the Left
!>          = 'R': apply Q or Q**T from the Right
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          = 'N': apply Q  (No transpose)
!>          = 'T': apply Q**T (Transpose)
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix C. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix C. N >= 0.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of elementary reflectors whose product defines
!>          the matrix Q.
!>          If SIDE = 'L', M >= K >= 0;
!>          if SIDE = 'R', N >= K >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,K)
!>          The i-th column must contain the vector which defines the
!>          elementary reflector H(i), for i = 1,2,...,k, as returned by
!>          DGEQRF in the first k columns of its array argument A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.
!>          If SIDE = 'L', LDA >= max(1,M);
!>          if SIDE = 'R', LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by DGEQRF.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension (LDC,N)
!>          On entry, the m by n matrix C.
!>          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C. LDC >= max(1,M).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension
!>                                   (N) if SIDE = 'L',
!>                                   (M) if SIDE = 'R'
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup unm2r
!
!  =====================================================================

      SUBROUTINE dorm2r( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &
     &                   WORK, INFO )
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      parameter( one = 1.0d+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LEFT, NOTRAN
      INTEGER            I, I1, I2, I3, IC, JC, MI, NI, NQ
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      EXTERNAL           lsame
!     ..
!     .. External Subroutines ..
!      EXTERNAL            xerbla, dlarf1f
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          max
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      info = 0
      left = lsame( side, 'L' )
      notran = lsame( trans, 'N' )
!
!     NQ is the order of Q
!
      IF( left ) THEN
         nq = m
      ELSE
         nq = n
      END IF
      IF( .NOT.left .AND. .NOT.lsame( side, 'R' ) ) THEN
         info = -1
      ELSE IF( .NOT.notran .AND. .NOT.lsame( trans, 'T' ) ) THEN
         info = -2
      ELSE IF( m.LT.0 ) THEN
         info = -3
      ELSE IF( n.LT.0 ) THEN
         info = -4
      ELSE IF( k.LT.0 .OR. k.GT.nq ) THEN
         info = -5
      ELSE IF( lda.LT.max( 1, nq ) ) THEN
         info = -7
      ELSE IF( ldc.LT.max( 1, m ) ) THEN
         info = -10
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DORM2R', -info )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( m.EQ.0 .OR. n.EQ.0 .OR. k.EQ.0 ) &
     &   RETURN
!
      IF( ( left .AND. .NOT.notran ) .OR. ( .NOT.left .AND. notran ) ) &
     &     THEN
         i1 = 1
         i2 = k
         i3 = 1
      ELSE
         i1 = k
         i2 = 1
         i3 = -1
      END IF
!
      IF( left ) THEN
         ni = n
         jc = 1
      ELSE
         mi = m
         ic = 1
      END IF
!
      DO 10 i = i1, i2, i3
         IF( left ) THEN
!
!           H(i) is applied to C(i:m,1:n)
!
            mi = m - i + 1
            ic = i
         ELSE
!
!           H(i) is applied to C(1:m,i:n)
!
            ni = n - i + 1
            jc = i
         END IF
!
!        Apply H(i)
!
         CALL dlarf1f( side, mi, ni, a( i, i ), 1, tau( i ), c( ic, &
     &               jc ), &
     &               ldc, work )
10    CONTINUE
      RETURN
!
!     End of DORM2R
!

      END


!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION DNRM2(N,X,INCX)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION X(*)
!       ..
!
!
!  =============
!
!  Arguments:
!  ==========
!
!
!  Authors:
!  ========
!
!
!
!
!  ==================
!
!  =====================
!  =====================================================================

function dnrm2( n, x, incx ) 
   integer, parameter :: wp = kind(1.d0)
   real(wp) :: dnrm2
!
!  -- Reference BLAS level1 routine (version 3.9.1) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     March 2021
!
!  .. Constants ..
   real(wp), parameter :: zero = 0.0_wp
   real(wp), parameter :: one  = 1.0_wp
   real(wp), parameter :: maxn = huge(0.0_wp)
!  ..
!  .. Blue's scaling constants ..
   real(wp), parameter :: tsml = real(radix(0._wp), wp)**ceiling( &
       (minexponent(0._wp) - 1) * 0.5_wp)
   real(wp), parameter :: tbig = real(radix(0._wp), wp)**floor( &
       (maxexponent(0._wp) - digits(0._wp) + 1) * 0.5_wp)
   real(wp), parameter :: ssml = real(radix(0._wp), wp)**( - floor( &
       (minexponent(0._wp) - digits(0._wp)) * 0.5_wp))
   real(wp), parameter :: sbig = real(radix(0._wp), wp)**( - ceiling( &
       (maxexponent(0._wp) + digits(0._wp) - 1) * 0.5_wp))
!  ..
!  .. Scalar Arguments ..
   integer :: incx, n
!  ..
!  .. Array Arguments ..
   real(wp) :: x(*)
!  ..
!  .. Local Scalars ..
   integer :: i, ix
   logical :: notbig
   real(wp) :: abig, amed, asml, ax, scl, sumsq, ymax, ymin
!
!  Quick return if possible
!
   dnrm2 = zero
   if( n <= 0 ) return
!
   scl = one
   sumsq = zero
!
!  Compute the sum of squares in 3 accumulators:
!     abig -- sums of squares scaled down to avoid overflow
!     asml -- sums of squares scaled up to avoid underflow
!     amed -- sums of squares that do not require scaling
!  The thresholds and multipliers are
!     tbig -- values bigger than this are scaled down by sbig
!     tsml -- values smaller than this are scaled up by ssml
!
   notbig = .true.
   asml = zero
   amed = zero
   abig = zero
   ix = 1
   if( incx < 0 ) ix = 1 - (n-1)*incx
   do i = 1, n
      ax = abs(x(ix))
      if (ax > tbig) then
         abig = abig + (ax*sbig)**2
         notbig = .false.
      else if (ax < tsml) then
         if (notbig) asml = asml + (ax*ssml)**2
      else
         amed = amed + ax**2
      end if
      ix = ix + incx
   end do
!
!  Combine abig and amed or amed and asml if more than one
!  accumulator was used.
!
   if (abig > zero) then
!
!     Combine abig and amed if abig > 0.
!
      if ( (amed > zero) .or. (amed > maxn) .or. (amed /= amed) ) then
         abig = abig + (amed*sbig)*sbig
      end if
      scl = one / sbig
      sumsq = abig
   else if (asml > zero) then
!
!     Combine amed and asml if asml > 0.
!
      if ( (amed > zero) .or. (amed > maxn) .or. (amed /= amed) ) then
         amed = sqrt(amed)
         asml = sqrt(asml) / ssml
         if (asml > amed) then
            ymin = amed
            ymax = asml
         else
            ymin = asml
            ymax = amed
         end if
         scl = one
         sumsq = ymax**2*( one + (ymin/ymax)**2 )
      else
         scl = one / ssml
         sumsq = asml
      end if
   else
!
!     Otherwise all values are mid-range
!
      scl = one
      sumsq = amed
   end if
   dnrm2 = scl*sqrt( sumsq )
   return

end function

end module lapack_in_fpp













