#include "CESR_platform.inc"

module bmad_basic_mod

  use dcslib
  use precision_def

  integer, parameter :: n_pole_maxx = 20  ! maximum multipole order

contains

!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!+
! Function field_interpolate_3d (position, field_mesh, deltas, position0)
!
! Function to interpolate a 3d field.
! The interpolation is such that the derivative is continuous.
!
! Note: For "interpolation" outside of the region covered by the field_mesh
! it is assumed that the field is constant, Equal to the field at the
! boundary.
!
! Modules needed:
!
! Input:
!   position(3)       -- Real(rp): (x, y, z) position.
!   field_mesh(:,:,:) -- Real(rp): Grid of field points.
!   deltas(3)         -- Real(rp): (dx, dy, dz) distances between mesh points.
!   position0(3)      -- Real(rp), optional:  position at (ix0, iy0, iz0) where
!                            (ix0, iy0, iz0) is the lower bound of the
!                            filed_mesh(i, j, k) array. If not present then
!                            position0 is taken to be (0.0, 0.0, 0.0)
! Output:
!   field_interpolate_3d -- Real(rp): interpolated field.
!-

real(rp) function field_interpolate_3d &
                             (position, field_mesh, deltas, position0)

  implicit none

  real(rp), optional, intent(in) :: position0(3)
  real(rp), intent(in) :: position(3), field_mesh(0:,0:,0:), deltas(3)

  real(rp) r(3), f(-1:2), g(-1:2), h(-1:2), r_frac(3)

  integer i0(3), ix, iy, iz, iix, iiy, iiz

!

  if (present(position0)) then
    r = (position - position0) / deltas
  else
    r = position / deltas
  endif

  i0 = int(r)
  r_frac = r - i0

  do ix = -1, 2
   iix = min(max(ix + i0(1), 0), ubound(field_mesh, 1))
   do iy = -1, 2
      iiy = min(max(iy + i0(2), 0), ubound(field_mesh, 2))
      do iz = -1, 2
        iiz = min(max(iz + i0(3), 0), ubound(field_mesh, 3))
        f(iz) = field_mesh(iix, iiy, iiz)
      enddo
      g(iy) = interpolate_1d (r_frac(3), f)
    enddo
    h(ix) = interpolate_1d (r_frac(2), g)
  enddo
  field_interpolate_3d = interpolate_1d (r_frac(1), h)

!---------------------------------------------------------------

contains

! interpolation in 1 dimension using 4 equally spaced points: P1, P2, P3, P4.
!   x = interpolation point.
!           x = 0 -> point is at P2.
!           x = 1 -> point is at P3.
! Interpolation is done so that the derivative is continuous.
! The interpolation uses a cubic polynomial

real function interpolate_1d (x, field_1d)

  implicit none

  real(rp) x, field_1d(4), df_2, df_3
  real(rp) c0, c1, c2, c3

!

  df_2 = (field_1d(3) - field_1d(1)) / 2   ! derivative at P2
  df_3 = (field_1d(4) - field_1d(2)) / 2   ! derivative at P3

  c0 = field_1d(2)
  c1 = df_2
  c2 = 3 * field_1d(3) - df_3 - 3 * field_1d(2) - 2 * df_2
  c3 = df_3 - 2 * field_1d(3) + 2 * field_1d(2) + df_2

  interpolate_1d = c0 + c1 * x + c2 * x**2 + c3 * x**3

end function

end function

!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!+
! Subroutine FitPoly(coe, x, y, order, samples)
!
! Subroutine to fit a polynomial, y = coe(0) + coe(1)*x + coe(2)*x^2 + ...,
! to the input data of x and y via least squares.
!
! Input:
!     x(:) -- Real(rp): vector of sample 'x' data
!     y(:) -- Real(rp): vector of sample 'y' data
!     order -- Integer: order of fitted polynomial
!     samples -- Integer: how many 'x, y' data samples
!
! Output:
!     coe(0:) -- Real(rp): array of polynomial coefficients
!-

subroutine fitpoly(coe, x, y, order, samples)

  implicit none

  integer maxcoe, maxsamp
  parameter(maxcoe=10, maxsamp=100)
  integer order, samples, numcoe
  real(rp) coe(0:), x(:), y(:)
  real(rp) Xmat(maxsamp,maxcoe), XtX(maxcoe,maxcoe), Xty(maxcoe)
  integer coe_index, sam_index, i, j, k

  numcoe = order + 1

  do sam_index = 1, samples
   Xmat(sam_index, 1) = 1.0
   do coe_index = 2, numcoe
    Xmat(sam_index, coe_index) = x(sam_index) *  &
    Xmat(sam_index, coe_index-1)
   enddo
  enddo

  do i = 1, numcoe
   Xty(i) = 0.0
   do j = 1, samples
    Xty(i) = Xty(i) + Xmat(j,i) * y(j)
   enddo
  enddo

  do i = 1, numcoe
   do j = 1, numcoe
    XtX(i,j) = 0.0
    do k = 1, samples
     XtX(i,j) = XtX(i,j) + Xmat(k,i) * Xmat(k,j)
    enddo
   enddo
  enddo

  call solvlin(XtX, Xty, coe, numcoe, maxcoe)

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Function c_multi (n, m)
!
! Subroutine to compute multipole factors:
!          c_multi(n, m) =  +/- ("n choose m")/n!
!
! Input:
!   n,m -- Integer: For n choose m
!
! Output:
!   c_multi  -- Real(rp): Multipole factor
!-

function c_multi (n, m) result (c_out)

  implicit none

  integer, intent(in) :: n, m
  integer in, im

  real(rp) c_out, factorial_n
  real(rp), save :: c(0:n_pole_maxx, 0:n_pole_maxx)

  logical, save :: init_needed = .true.

! The magnitude of c(n, m) is number of combinations normalized by n!

  if (init_needed) then

    c(0, 0) = 1

    do in = 1, n_pole_maxx
      c(in, 0) = 1
      c(in, in) = 1
      do im = 1, in-1
        c(in, im) = c(in-1, im-1) + c(in-1, im)
      enddo
    enddo

    factorial_n = 1

    do in = 0, n_pole_maxx
      if (in > 0) factorial_n = in * factorial_n
      do im = 0, in
        c(in, im) = c(in, im) / factorial_n
        if (mod(im, 4) == 0) c(in, im) = -c(in, im)
        if (mod(im, 4) == 3) c(in, im) = -c(in, im)
      enddo
    enddo

    init_needed = .false.

  endif

!

  c_out = c(n, m)

end function

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

function mexp (x, m)

  implicit none

  real(rp) x, mexp
  integer m

!

  if (m < 0) then
    mexp = 0
  elseif (m == 0) then
    mexp = 1
  else
    mexp = x**m
  endif

end function

end module
