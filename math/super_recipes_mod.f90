module super_recipes_mod

use precision_def
use output_mod
use re_allocate_mod
use swap_mod

! super_mrqmin_storage_struct is used to:
!   1) Make super_mrqmin thread safe.
!   2) Make super_mrqmin recursive. Needed, for example, in the Tao program:
!       Tao -> super_mrqmin -> tao_funcs -> tao_merit -> closed_orbit_calc -> super_mrqmin

type super_mrqmin_storage_struct
  ! Used by super_mrqmin
  real(rp), allocatable :: covar(:, :)   ! Covariance matrix. See mrqmin in NR for more details.
  real(rp), allocatable :: alpha(:, :)   ! Curvature matrix. See mrqmin in NR for more details.
  real(rp), allocatable :: atry(:), beta(:)
  real(rp), allocatable :: da(:,:)
  real(rp) :: ochisq
  logical, allocatable :: mask(:)
  ! Used by super_mrqcof
  real(rp), allocatable :: dyda(:, :)
  real(rp), allocatable :: old_dy(:), dy(:), wt(:), ymod(:)
end type

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine super_bicubic_interpolation(y, y1, y2, y12, x1l, x1u, x2l, x2u, x1, x2, ansy, ansy1, ansy2)
!
! Routine to do bicubic interpolation.
! This is from NR bcuint.
!
! Note! The four grid points are arrayed in counter-clockwise order beginning from the lower left.
! So, for example, y = [y_ll, y_lu, y_uu, y_ul] where "l" = lower, "u" = upper index.
!
! Input:
!   y(4)    -- real(rp): Function values at grid points.
!   y1(4)   -- real(rp): dy/dx1 derivatives.
!   y2(4)   -- real(rp): dy/dx2 derivatives.
!   y12(4)  -- real(rp): d2y/dx1*dx2 second derivatives.
!   x1l     -- real(rp): 1-direction coordinate at lower points.
!   x1u     -- real(rp): 1-direction coordinate at upper points
!   x2l     -- real(rp): 2-direction coordinate at lower points.
!   x2u     -- real(rp): 2-direction coordinate at upper points
!   x1      -- real(rp): 1-direction coordinate at point to evaluate.
!   x2      -- real(rp): 2-direction coordinate at point to evaluate.
!
! Output:
!   ansy    -- real(rp): Interpolation value.
!   ansy1   -- real(rp): 1-direction derivative at interpolation point.
!   ansy2   -- real(rp): 2-direction derivative at interpolation point.
!-

subroutine super_bicubic_interpolation(y, y1, y2, y12, x1l, x1u, x2l, x2u, x1, x2, ansy, ansy1, ansy2)

implicit none

real(rp), dimension(4), intent(in) :: y, y1, y2, y12
real(rp), intent(in) :: x1l, x1u, x2l, x2u, x1, x2
real(rp), intent(out) :: ansy, ansy1, ansy2
integer :: i
real(rp) :: t, u
real(rp), dimension(4,4) :: c

!

call super_bicubic_coef(y, y1, y2, y12, x1u-x1l, x2u-x2l, c)

if (x1u == x1l .or. x2u == x2l) call err_exit ('super_bicubic_interpolation: problem with input values - boundary pair equal?')
t = (x1-x1l)/(x1u-x1l)
u = (x2-x2l)/(x2u-x2l)
ansy = 0.0
ansy2 = 0.0
ansy1 = 0.0

do i = 4, 1, -1
  ansy = t*ansy + ((c(i, 4)*u + c(i, 3))*u + c(i, 2))*u + c(i, 1)
  ansy2 = t*ansy2 + (3.0_rp*c(i,4)*u + 2.0_rp*c(i, 3))*u + c(i, 2)
  ansy1 = u*ansy1 + (3.0_rp*c(4,i)*t + 2.0_rp*c(3, i))*t + c(2, i)
end do

ansy1 = ansy1/(x1u-x1l)
ansy2 = ansy2/(x2u-x2l)

end subroutine super_bicubic_interpolation

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine super_bicubic_coef(y, y1, y2, y12, d1, d2, c)
!
! Routine to compute coefficients for bicubic interpolation.
! This is from NR bcucof.
!
! Input:
!   y(4)    -- real(rp): Function values at grid points.
!   y1(4)   -- real(rp): dy/dx1 derivatives.
!   y2(4)   -- real(rp): dy/dx2 derivatives.
!   y12(4)  -- real(rp): d2y/dx1*dx2 second derivatives.
!   d1      -- real(rp): Grid width in 1-direction.
!   d2      -- real(rp): Grid width in 2-direction.
!
! Output:
!   c(4,4)  -- real(rp): Coefficients.
!-

subroutine super_bicubic_coef(y, y1, y2, y12, d1, d2, c)

implicit none

real(dp), intent(in) :: d1, d2
real(dp), dimension(4), intent(in) :: y, y1, y2, y12
real(dp), dimension(4,4), intent(out) :: c
real(dp), dimension(16) :: x
real(dp), dimension(16,16) :: wt

!

data wt /1, 0, -3, 2, 4*0, -3, 0, 9, -6, 2, 0, -6, 4, &
  8*0, 3, 0, -9, 6, -2, 0, 6, -4, 10*0, 9, -6, 2*0, -6, 4, 2*0, 3, -2, 6*0, -9, 6, &
  2*0, 6, -4, 4*0, 1, 0, -3, 2, -2, 0, 6, -4, 1, 0, -3, 2, 8*0, -1, 0, 3, -2, 1, 0, -3, &
  2, 10*0, -3, 2, 2*0, 3, -2, 6*0, 3, -2, 2*0, -6, 4, 2*0, 3, -2, 0, 1, -2, 1, 5*0, &
  -3, 6, -3, 0, 2, -4, 2, 9*0, 3, -6, 3, 0, -2, 4, -2, 10*0, -3, 3, 2*0, 2, -2, 2*0, &
  -1, 1, 6*0, 3, -3, 2*0, -2, 2, 5*0, 1, -2, 1, 0, -2, 4, -2, 0, 1, -2, 1, 9*0, -1, 2, &
  -1, 0, 1, -2, 1, 10*0, 1, -1, 2*0, -1, 1, 6*0, -1, 1, 2*0, 2, -2, 2*0, -1, 1/

!

x(1:4) = y
x(5:8) = y1*d1
x(9:12) = y2*d2
x(13:16) = y12*d1*d2
x = matmul(wt, x)
c = reshape(x, [4,4], order = [2,1])

end subroutine super_bicubic_coef

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine super_sort(arr)
!
! Routine to sort an integer array in place.
! This is the NR routine sort modified to sort integers.
!
! Input:
!   arr(:)      -- integer: Array of integers.
!
! Output:
!   arr(:)      -- integer: Sorted array.
!-

subroutine super_sort(arr)

implicit none

integer, dimension(:), intent(inout) :: arr
integer, parameter :: nn = 15, nstack = 50
integer :: a
integer :: n, k, i, j, jstack, l, r
integer, dimension(nstack) :: istack

!

n = size(arr)
jstack = 0
l = 1
r = n
do
  if (r-l < NN) then
    do j = l+1, r
      a = arr(j)
      do i = j-1, l, -1
        if (arr(i) <= a) exit
        arr(i+1) = arr(i)
      end do
      arr(i+1) = a
    end do
    if (jstack == 0) RETURN
    r = istack(jstack)
    l = istack(jstack-1)
    jstack = jstack-2
  else
    k = (l+r)/2
    call swap(arr(k), arr(l+1))
    if (arr(l)>arr(r))   call swap(arr(l), arr(r))
    if (arr(l+1)>arr(r)) call swap(arr(l+1), arr(r))
    if (arr(l)>arr(l+1)) call swap(arr(l), arr(l+1))
    i = l+1
    j = r
    a = arr(l+1)
    do
      do
        i = i+1
        if (arr(i) >= a) exit
      end do
      do
        j = j-1
        if (arr(j) <= a) exit
      end do
      if (j < i) exit
      call swap(arr(i), arr(j))
    end do
    arr(l+1) = arr(j)
    arr(j) = a
    jstack = jstack+2
    if (jstack > NSTACK) call err_exit('sort: NSTACK too small')
    if (r-i+1 >= j-l) then
      istack(jstack) = r
      istack(jstack-1) = i
      r = j-1
    else
      istack(jstack) = j-1
      istack(jstack-1) = l
      l = i
    end if
  end if
end do

end subroutine super_sort

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Function super_rtsafe (funcs, x1, x2, rel_tol, abs_tol, status) result (x_zero)
!
! Routine find the root of a function. Use this method when computing derivatives is easy.
! Otherwise, use super_zbrent.
!
! The root is considered found if the root is known to a tolerance x_tol given by:
!   x_tol = |x_zero| * rel_tol + abs_tol
!
! This routine is essentially rtsafe from Numerical Recipes with the feature that it 
! returns a status flag.
!
! Input:
!   funcs     -- Function whose root is to be found. The interface is:
!                  subroutine funcs(x, fval, fderiv, status)
!                    real(rp), intent(in) :: x
!                    real(rp), intent(out) :: fval, fderiv
!                    integer status
!                  end subroutine funcs
!
!   x1, x2    -- Real(rp): Bracket values.
!   rel_tol   -- real(rp): Relative tolerance for the error of the root.
!   abs_tol   -- real(rp): Absolute tolerance for the error of the root. This is xacc arg in rtsafe.
!
! Output:
!   x_zero    -- Real(rp): Root found.
!   status    -- Integer: Calculation status:
!                      -2    => Max iterations exceeded.
!                      -1    => Root not bracketed.
!                       0    => Normal.
!                       Other => Set by funcs. 
!-

function super_rtsafe (funcs, x1, x2, rel_tol, abs_tol, status) result (x_zero)

implicit none

interface
  subroutine funcs(x, fval, fderiv, status)
  import
  implicit none
  real(rp), intent(in) :: x
  real(rp), intent(out) :: fval, fderiv
  integer status
  end subroutine funcs
end interface

real(rp), intent(in) :: x1, x2, rel_tol, abs_tol
real(rp) :: x_zero
real(rp) :: df, dx, dxold, f, fh, fl, temp, xh, xl
integer, parameter :: maxit = 100
integer :: status, j

!

status = 0
x_zero = real_garbage$

call funcs(x1, fl, df, status); if (status /= 0) return
call funcs(x2, fh, df, status); if (status /= 0) return
if ((fl > 0.0 .and. fh > 0.0) .or. (fl < 0.0 .and. fh < 0.0)) then
  status = -1
  return
endif

if (fl == 0.0) then
  x_zero=x1
  return
else if (fh == 0.0) then
  x_zero = x2
  return
else if (fl < 0.0) then
  xl = x1
  xh = x2
else
  xh = x1
  xl = x2
end if

x_zero = 0.5_rp*(x1+x2)
dxold = abs(x2-x1)
dx = dxold

call funcs(x_zero, f, df, status); if (status /= 0) return

do j = 1, maxit
  if (((x_zero-xh)*df-f)*((x_zero-xl)*df-f) > 0.0 .or. abs(2.0_rp*f) > abs(dxold*df) ) then
    dxold = dx
    dx = 0.5_rp*(xh-xl)
    x_zero = xl+dx
    if (xl == x_zero) return
  else
    dxold = dx
    dx = f/df
    temp = x_zero
    x_zero = x_zero-dx
    if (temp == x_zero) return
  end if
  if (abs(dx) < rel_tol * abs(x_zero) + abs_tol) return
  call funcs(x_zero, f, df, status); if (status /= 0) return
  if (f < 0.0) then
    xl = x_zero
  else
    xh = x_zero
  end if
end do

status = -2

end function super_rtsafe

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Function super_mnbrak (ax, bx, cx, fa, fb, fc, func, status)
!
! Given initial range [ax, bx], this routine finds new points ax, bx, cx that bracket a minimum of 
! the function func.
!
! This routine is essentially mnbrak from Numerical Recipes with the feature that it returns
! a status integer.
!
! Input:
!   ax, bx     -- real(rp): Range of x to search. It is permitted that bx < ax.
!   func       -- function whose minimum is to be bracketed. The interface is:
!                    function func(x, status) result (value)
!                      real(rp), intent(in) :: x
!                      integer, optional :: status  ! If non-zero return value, super_mnbrak will terminate.
!                      real(rp) :: value
!                    end function func
!
! Output:
!   ax, bx      -- real(rp): Range over which there is a minimum.
!   cx          -- real(rp): Value in the range [ax, bx] such that func(cx) < min(func(ax), func(cx)).
!   fa, fb, fc  -- real(rp): values of func(ax), func(bx), and func(cx) respectively.
!   status      -- integer: Calculation status:
!                       0    => Normal.
!                       Other => If set by func. 
!-

subroutine super_mnbrak(ax, bx, cx, fa, fb, fc, func, status)

implicit none

real(rp), intent(inout) :: ax, bx
real(rp), intent(out) :: cx, fa, fb, fc

interface
  function func(x, status) result (value)
  import
  implicit none
  real(rp), intent(in) :: x
  real(rp) :: value
  integer, optional :: status
  end function func
end interface

real(rp), parameter :: gold = 1.618034_rp, glimit = 100.0_rp, tiny = 1.0e-20_rp
real(rp) :: fu, q, r, u, ulim
integer status

!

status = 0

fa = func(ax, status)
fb = func(bx, status)
if (fb > fa) then
  call swap(ax, bx)
  call swap(fa, fb)
end if
cx = bx+gold*(bx-ax)
fc = func(cx, status); if (status /= 0) return
do
  if (fb < fc) return
  r = (bx-ax)*(fb-fc)
  q = (bx-cx)*(fb-fa)
  u = bx-((bx-cx)*q-(bx-ax)*r)/(2.0_rp*sign(max(abs(q-r), tiny), q-r))
  ulim = bx+glimit*(cx-bx)
  if ((bx-u)*(u-cx) > 0.0) then
    fu = func(u, status); if (status /= 0) return
    if (fu < fc) then
      ax = bx
      fa = fb
      bx = u
      fb = fu
      return
    else if (fu > fb) then
      cx = u
      fc = fu
      return
    end if
    u = cx+gold*(cx-bx)
    fu = func(u, status); if (status /= 0) return
  else if ((cx-u)*(u-ulim) > 0.0) then
    fu = func(u, status); if (status /= 0) return
    if (fu < fc) then
      bx = cx
      cx = u
      u = cx+gold*(cx-bx)
      call shft(fb, fc, fu, func(u, status))
    end if
  else if ((u-ulim)*(ulim-cx) >= 0.0) then
    u = ulim
    fu = func(u, status); if (status /= 0) return
  else
    u = cx+gold*(cx-bx)
    fu = func(u, status); if (status /= 0) return
  end if
  call shft(ax, bx, cx, u)
  call shft(fa, fb, fc, fu)
end do

!------------------------
contains

subroutine shft(a, b, c, d)
real(rp), intent(out) :: a
real(rp), intent(inout) :: b, c
real(rp), intent(in) :: d
a = b
b = c
c = d
end subroutine shft

end subroutine super_mnbrak

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Function super_brent (ax, bx, cx, func, rel_tol, abs_tol, xmin, status) result (ymin)
!
! Routine to find the minimum of a function.
! The x-tolerance is:
!   x-tolerance = |xmin| * rel_tol + abs_tol
!
! This routine is essentially brent from Numerical Recipes with the feature that it returns
! a status integer if something goes wrong instead of bombing.
!
! Input:
!   ax, cx   -- real(rp): Range of x to search. It is permitted that cx < ax.
!   bx       -- real(rp): x-value between ax and cx such that func(bx) < min(func(ax), func(cx)).
!   func     -- function whose minimum is to be found. The interface is:
!                  function func(x, status) result (value)
!                    real(rp), intent(in) :: x
!                    integer, optional :: status  ! If non-zero return value, super_brent will terminate.
!                    real(rp) :: value
!                  end function func
!   rel_tol  -- real(rp): Relative tolerance for the error xmin. This is tol in NR brent routine.
!   abs_tol  -- real(rp): Absolute tolerance for the error xmin. Equivalent to 1e-18*abs(ax) in NR brent.
!
! Output:
!   xmin     -- real(rp): x-coordinate at minimum.
!   status   -- integer: Calculation status:
!                      -2    => Max iterations exceeded.
!                       0    => Normal.
!                       Other => Set by funcs. 
!   ymin     -- real(rp) value of func(xmin).
!-

function super_brent(ax, bx, cx, func, rel_tol, abs_tol, xmin, status) result (ymin)

implicit none

real(rp), intent(in) :: ax, bx, cx, rel_tol, abs_tol
real(rp), intent(out) :: xmin
real(rp) ymin
real(rp) :: brent

interface
  function func(x, status)
  import
  implicit none
  real(rp), intent(in) :: x
  real(rp) :: func
  integer, optional :: status
  end function func
end interface

real(rp), parameter :: cgold = 0.3819660_dp, zeps = 1.0e-3_dp*epsilon(ax)
integer, parameter :: itmax = 100

real(rp) tol1, tol2, a, b, d, e, etemp, fu, fv, fw, fx, p, q, r, u, v, w, x, xm
integer :: iter, status

!

status = 0
ymin = real_garbage$

a = min(ax, cx)
b = max(ax, cx)
v = bx
w = v
x = v
e = 0.0
fx = func(x, status); if (status /= 0) return
fv = fx
fw = fx
do iter = 1, itmax
  xm = 0.5_dp*(a+b)
  tol1 = rel_tol*abs(x)+abs_tol
  tol2 = 2.0_dp*tol1
  if (abs(x-xm) <= (tol2-0.5_dp*(b-a))) then
    xmin = x
    ymin = fx
    return
  end if
  if (abs(e) > tol1) then
    r = (x-w)*(fx-fv)
    q = (x-v)*(fx-fw)
    p = (x-v)*q-(x-w)*r
    q = 2.0_dp*(q-r)
    if (q > 0.0) p = -p
    q = abs(q)
    etemp = e
    e = d
    if (abs(p) >= abs(0.5_dp*q*etemp) .or. &
      p <= q*(a-x) .or. p >= q*(b-x)) then
      e = merge(a-x, b-x, x >= xm )
      d = cgold*e
    else
      d = p/q
      u = x+d
      if (u-a < tol2 .or. b-u < tol2) d = sign(tol1, xm-x)
    end if
  else
    e = merge(a-x, b-x, x >= xm )
    d = cgold*e
  end if
  u = merge(x+d, x+sign(tol1, d), abs(d) >= tol1 )
  fu = func(u, status); if (status /= 0) return
  if (fu <= fx) then
    if (u >= x) then
      a = x
    else
      b = x
    end if
    call shft(v, w, x, u)
    call shft(fv, fw, fx, fu)
  else
    if (u < x) then
      a = u
    else
      b = u
    end if
    if (fu <= fw .or. w == x) then
      v = w
      fv = fw
      w = u
      fw = fu
    else if (fu <= fv .or. v == x .or. v == w) then
      v = u
      fv = fu
    end if
  end if
end do

status = -2
return

!--------------------------------------
contains
subroutine shft(a, b, c, d)
real(rp), intent(out) :: a
real(rp), intent(inout) :: b, c
real(rp), intent(in) :: d
a = b
b = c
c = d
end subroutine shft

end function super_brent

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Function super_dbrent (ax, bx, cx, func, dfunc, rel_tol, abs_tol, xmin, status) result (func_min)
!
! Rotine to find the minimum of a function in the range [ax, cx].
! This if based on the NR routine dbrent.
!
! Input:
!   ax        -- real(rp): Minimum of the range to search.
!   bx        -- real(rp): Value in the range [ax, cx] such that func(bx) < min(func(ax), func(bx)).
!   cx        -- real(rp): Maximum of the range to search.
!   func      -- function whose minimum is to be found. The interface is:
!                  function func(x, status) result (value)
!                    real(rp), intent(in) :: x
!                    integer, optional :: status  ! If non-zero return value, super_dbrent will terminate.
!                    real(rp) :: value
!                  end function func
!   dfunc     -- function derivative. The interface is:
!                  function dfunc(x, status) result (value)
!                    real(rp), intent(in) :: x
!                    integer, optional :: status  ! If non-zero return value, super_dbrent will terminate.
!                    real(rp) :: dvalue
!                  end function dfunc
!   rel_tol   -- real(rp): Relative tolerance for the error xmin. This is tol in NR brent routine.
!   abs_tol   -- real(rp): Absolute tolerance for the error xmin. Equivalent to 1e-18*abs(ax) in NR dbrent.
!
! Output:
!   xmin      -- real(rp): x-coordinate at minimum.
!   status    -- integer: Calculation status:
!                      -2    => Max iterations exceeded.
!                       0    => Normal.
!                       Other => Set by funcs or dfuncs. 
!   func_min  -- real(rp) value of func(xmin).
!-

function super_dbrent(ax, bx, cx, func, dfunc, rel_tol, abs_tol, xmin, status) result (func_min)

implicit none
real(rp), intent(in) :: ax, bx, cx
real(rp), intent(out) :: xmin
real(rp) :: rel_tol, abs_tol, func_min

interface
  function func(x, status) result (dvalue)
  import
  implicit none
  real(rp) :: x
  real(rp) :: dvalue
  integer, optional :: status
  end function func

  function dfunc(x, status) result (value)
  import
  implicit none
  real(rp) :: x
  real(rp) :: value
  integer, optional :: status
  end function dfunc
end interface

integer, parameter :: itmax = 100
integer :: status, iter
real(rp) a, b, d, d1, d2, du, dv, dw, dx, e, fu, fv, fw, fx, olde, tol1, tol2
real(rp) u, u1, u2, v, w, x, xm
logical :: ok1, ok2

!

a = min(ax, cx)
b = max(ax, cx)
v = bx
w = v
x = v
e = 0.0
fx = func(x, status); if (status /= 0) return
fv = fx
fw = fx
dx = dfunc(x, status); if (status /= 0) return
dv = dx
dw = dx

do iter = 1, itmax
  xm = 0.5_rp*(a+b)
  tol1 = rel_tol*abs(x)+abs_tol
  tol2 = 2.0_rp*tol1
  if (abs(x-xm) <= (tol2-0.5_rp*(b-a))) exit
  if (abs(e) > tol1) then
    d1 = 2.0_rp*(b-a)
    d2 = d1
    if (dw /= dx) d1 = (w-x)*dx/(dx-dw)
    if (dv /= dx) d2 = (v-x)*dx/(dx-dv)
    u1 = x+d1
    u2 = x+d2
    ok1 = ((a-u1)*(u1-b) > 0.0) .and. (dx*d1 <= 0.0)
    ok2 = ((a-u2)*(u2-b) > 0.0) .and. (dx*d2 <= 0.0)
    olde = e
    e = d
    if (ok1 .or. ok2) then
      if (ok1 .and. ok2) then
        d = merge(d1, d2, abs(d1) < abs(d2))
      else
        d = merge(d1, d2, ok1)
      end if
      if (abs(d) <= abs(0.5_rp*olde)) then
        u = x+d
        if (u-a < tol2 .or. b-u < tol2) &
          d = sign(tol1, xm-x)
      else
        e = merge(a, b, dx >= 0.0)-x
        d = 0.5_rp*e
      end if
    else
      e = merge(a, b, dx >= 0.0)-x
      d = 0.5_rp*e
    end if
  else
    e = merge(a, b, dx >= 0.0)-x
    d = 0.5_rp*e
  end if
  if (abs(d) >= tol1) then
    u = x+d
    fu = func(u, status); if (status /= 0) return
  else
    u = x+sign(tol1, d)
    fu = func(u, status); if (status /= 0) return
    if (fu > fx) exit
  end if
  du = dfunc(u, status); if (status /= 0) return
  if (fu <= fx) then
    if (u >= x) then
      a = x
    else
      b = x
    end if
    call mov3(v, fv, dv, w, fw, dw)
    call mov3(w, fw, dw, x, fx, dx)
    call mov3(x, fx, dx, u, fu, du)
  else
    if (u < x) then
      a = u
    else
      b = u
    end if
    if (fu <= fw .or. w  == x) then
      call mov3(v, fv, dv, w, fw, dw)
      call mov3(w, fw, dw, u, fu, du)
    else if (fu <= fv .or. v  == x .or. v  == w) then
      call mov3(v, fv, dv, u, fu, du)
    end if
  end if
end do

status = -2
xmin = x
func_min = fx

!---------------------------------------------
contains

subroutine mov3(a, b, c, d, e, f)
real(rp), intent(in) :: d, e, f
real(rp), intent(out) :: a, b, c
a = d
b = e
c = f
end subroutine mov3

end function super_dbrent

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Function super_zbrent (func, x1, x2, rel_tol, abs_tol, status) result (x_zero)
!
! Routine to find the root of a function.
!
! The x-tolerance is:
!   x-tolerance = |x_root| * rel_tol + abs_tol
!
! Consider using super_rtsafe if computing derivatives is easy.
!
! This routine is essentially zbrent from Numerical Recipes with the feature that it returns
! a status integer if something goes wrong instead of bombing.
!
! Input:
!   func     -- function whose root is to be found. The interface is:
!                  function func(x, status) result (value)
!                    real(rp), intent(in) :: x
!                    integer :: status    ! If non-zero return value, super_zbrent will terminate.
!                    real(rp) :: value
!                  end function func
!   x1, x2   -- real(rp): Bracket values. Note: x2 does not have to be larger than x1
!   rel_tol  -- real(rp): Relative tolerance for x_zero. This is 4d-16 in NR zbrent.
!   abs_tol  -- real(rp): Absolute tolerance for x_zero. This is tol in NR zbrent.
!
! Output:
!   x_zero   -- real(rp): Root found.
!   status   -- integer: Calculation status:
!                      -2    => Max iterations exceeded.
!                      -1    => Root not bracketed.
!                       0    => Normal.
!                       Other => Set by func. 
!-

function super_zbrent (func, x1, x2, rel_tol, abs_tol, status) result (x_zero)

implicit none

real(rp), intent(in) :: x1, x2, rel_tol, abs_tol
real(rp) :: x_zero

interface
  function func(x, status) result (value)
    import
    implicit none
    real(rp), intent(in) :: x
    integer status
    real(rp) :: value
  end function func
end interface

integer, parameter :: itmax = 100
integer :: status, iter
real(rp) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm

character(*), parameter :: r_name = 'super_zbrent'

!

status = 0
x_zero = real_garbage$

a = x1
b = x2

fa = func(a, status); if (status /= 0) return
fb = func(b, status); if (status /= 0) return

if ((fa > 0.0 .and. fb > 0.0) .or. (fa < 0.0 .and. fb < 0.0)) then
  call out_io (s_fatal$, r_name, 'ROOT NOT BRACKETED!, \es12.4\ at \es12.4\ and \es12.4\ at \es12.4\ ', &
                  r_array = [fa, a, fb, b])
  x_zero = 1d100
  status = -1
  return
endif

c = b
fc = fb

do iter = 1,ITMAX
  if ((fb > 0.0 .and. fc > 0.0) .or. (fb < 0.0 .and. fc < 0.0)) then
    c = a
    fc = fa
    d = b-a
    e = d
  end if

  if (abs(fc) < abs(fb)) then
    a = b
    b = c
    c = a
    fa = fb
    fb = fc
    fc = fa
  end if

  tol1 = 0.5_rp * (rel_tol * abs(b) + abs_tol)
  xm = 0.5_rp*(c-b)

  if (abs(xm) <= tol1 .or. fb == 0.0) then
    !! x_zero = b
    if (fb == 0) then
      x_zero = b
    else
      x_zero = (b * fc - c * fb) / (fc - fb)   ! Linear interpolation.
    endif
    return
  end if

  if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
    s = fb/fa
    if (a == c) then
      p = 2.0_rp*xm*s
      q = 1.0_rp-s
    else
      q = fa/fc
      r = fb/fc
      p = s*(2.0_rp*xm*q*(q-r)-(b-a)*(r-1.0_rp))
      q = (q-1.0_rp)*(r-1.0_rp)*(s-1.0_rp)
    end if
    if (p > 0.0) q = -q
    p = abs(p)
    if (2.0_rp*p  <  min(3.0_rp*xm*q-abs(tol1*q),abs(e*q))) then
      e = d
      d = p/q
    else
      d = xm
      e = d
    end if
  else
    d = xm
    e = d
  end if

  a = b
  fa = fb
  b = b+merge(d,sign(tol1,xm), abs(d) > tol1 )
  fb = func(b, status); if (status /= 0) return
end do

call out_io (s_fatal$, r_name, 'EXCEEDED MAXIMUM ITERATIONS!')
status = -2
x_zero = b

end function super_zbrent 

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! subroutine super_mrqmin (y, weight, a, chisq, funcs, storage, alamda, status, maska)
!
! Routine to do non-linear optimizations. 
! This routine is essentially mrqmin from Numerical Recipes with some added features and
! some code tweaking to make the code run faster.
!
! Input:
!   y(:)        -- Real(rp): Data to fit to. See mrqmin in NR for more details.
!   weight(:)   -- Real(rp): This is equivalent to the 1/sig^2 of mrqmin in NR.
!   a(:)        -- Real(rp): Variables to vary. See mrqmin in NR for more details.
!   funcs       -- Function: User supplied subroutine. See mrqmin in NR for more details.
!                   The interface is:
!                        subroutine funcs(a, yfit, dyda, status)
!                          real(rp), intent(in) :: a(:)
!                          real(rp), intent(out) :: yfit(:)
!                          real(rp), intent(out) :: dyda(:, :)
!                          integer status
!                        end subroutine funcs
!                   Note: If funcs sets the status argument to anything non-zero, 
!                   super_mrqmin will halt the calculation and return back to the 
!                   calling routine. funcs should use positive values for status to
!                   avoid conflict with gaussj. 
!   storage     -- super_mrqmin_storage_struct: Storage for work space variables needed to be 
!                    retained by super_mrqmin between calls.
!   alamda      -- Real(rp): See mrqmin in NR for more details.
!   maska(:)    -- Logical, optional: Variable mask. See mrqmin in NR for more details.
!                    Default is True for all elements of the array.
!
! Output:
!   a(:)        -- Real(rp): Variables to vary. See mrqmin in NR for more details.
!   chisq       -- Real(rp): Chi^2 figure of merit. See mrqmin in NR for more details.
!   storage     -- super_mrqmin_storage_struct: 
!     %alpha        -- Curvature matrix. See mrqmin in NR for more details.
!     %covar        -- Covariance matrix. See mrqmin in NR for more details.
!   alamda      -- Real(rp): See mrqmin in NR for more details.
!   status      -- Integer: Calculation status. Any status /= 0 will cause super_mrqmin to stop.
!                      -2    => Singular matrix error in gaussj routine.
!                      -1    => Singular matrix error in gaussj routine.
!                       0    => Normal.
!                       Other => Set by funcs. 
!-

recursive subroutine super_mrqmin (y, weight, a, chisq, funcs, storage, alamda, status, maska)

implicit none

type (super_mrqmin_storage_struct) storage

real(rp) :: y(:), weight(:)
real(rp) :: a(:)
real(rp) :: chisq
real(rp) :: alamda
integer :: ma, ndata
integer status, mfit, j

logical, intent(in), optional :: maska(:)

character(*), parameter :: r_name = 'super_mrqmin'

interface
  subroutine funcs(a, yfit, dyda, status)
    import
    real(rp), intent(in) :: a(:)
    real(rp), intent(out) :: yfit(:)
    real(rp), intent(out) :: dyda(:, :)
    integer status
  end subroutine funcs
end interface

!

ndata = assert_equal([size(y), size(weight)], 'super_mrqmin: ndata')
ma = size(a)

call re_allocate(storage%mask, size(a))
if (present(maska)) then
  ma = assert_equal([ma, size(maska)], 'super_mrqmin: maska')
  storage%mask = maska
else
  storage%mask = .true.
endif

status = 0
mfit = count(storage%mask)

if (alamda < 0.0) then
  call re_allocate(storage%atry, ma)
  call re_allocate(storage%beta, ma)
  call re_allocate2d (storage%da, ma, 1)
  call re_allocate2d(storage%covar, ma, ma)
  call re_allocate2d(storage%alpha, ma, ma)
  call re_allocate2d(storage%dyda, ndata, ma)
  call re_allocate(storage%old_dy, ndata)
  call re_allocate(storage%dy, ndata, init_val = 0.0_rp)
  call re_allocate(storage%wt, ndata)
  call re_allocate(storage%ymod, ndata)
  alamda = 0.001_rp
  call super_mrqcof(a, y, storage%alpha, storage%beta, weight, chisq, funcs, storage, status)
  if (status /= 0) then
    deallocate(storage%atry, storage%beta, storage%da, storage%covar, storage%alpha)
    deallocate(storage%dyda, storage%dy, storage%wt, storage%ymod)
    return
  endif
  storage%ochisq = chisq
  storage%old_dy = storage%dy
  storage%atry = a
else
  ! If the calling program is changing weights then this will potentially confuse things.
  ! To prevent this, recalculate the original chisq.
  storage%ochisq = dot_product(storage%dy**2, weight)
end if

storage%covar(1:mfit, 1:mfit) = storage%alpha(1:mfit, 1:mfit)
forall (j = 1:mfit) storage%covar(j,j) =  storage%covar(j,j) * (1.0_rp+alamda)
storage%da(1:mfit, 1) = storage%beta(1:mfit)
call super_gaussj(storage%covar(1:mfit, 1:mfit), storage%da(1:mfit, 1:1), status)
if (status /= 0) then
  call out_io (s_error$, r_name, 'Note: Generally a singular matrix means that one or more datum values are', &
                                 'not affected by any variation of any variable.')
  return
endif

if (alamda == 0.0) then
  call covar_expand(storage%covar, storage%mask)
  call covar_expand(storage%alpha, storage%mask)
  deallocate(storage%atry, storage%beta, storage%da)
  deallocate(storage%dyda, storage%dy, storage%wt, storage%ymod)
  return
end if

storage%atry = a+unpack(storage%da(1:mfit, 1), storage%mask, 0.0_rp)
call super_mrqcof(storage%atry, y, storage%covar, storage%da(1:mfit, 1), weight, chisq, funcs, storage, status)

! Increase alamda by 2 (Instead of 10 as in NR version) gives better convergence. See:
!   "The Geometry of Nonlinear Least Squares, with applications to Sloppy Models and Optimization"
!   Mark K Transtrum, et. al.

if (chisq < storage%ochisq .and. status == 0) then
  alamda = 0.1_rp * alamda
  storage%ochisq = chisq
  storage%alpha(1:mfit, 1:mfit) = storage%covar(1:mfit, 1:mfit)
  storage%beta(1:mfit) = storage%da(1:mfit, 1)
  a = storage%atry
else
  alamda = 2.0_rp * alamda
  chisq = storage%ochisq
  storage%dy = storage%old_dy
end if

end subroutine super_mrqmin

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine super_mrqcof (a, y, co_alpha, da_beta, weight, chisq, funcs, storage, status)
! 
! Routine used by super_mrqmin. Not meant for general use.
!-

recursive subroutine super_mrqcof (a, y, co_alpha, da_beta, weight, chisq, funcs, storage, status)

implicit none

type (super_mrqmin_storage_struct) storage

real(rp) :: y(:), a(:), weight(:)
real(rp) :: da_beta(:)
real(rp) :: co_alpha(:, :)
real(rp) chisq

integer :: j, k, l, m, nv, nd
integer status

interface
  subroutine funcs(a, yfit, dyda, status)
    import
    real(rp), intent(in) :: a(:)
    real(rp), intent(out) :: yfit(:)
    real(rp), intent(out) :: dyda(:, :)
    integer status
  end subroutine funcs
end interface

!

nd = size(weight)
nv = size(a)

call funcs(a, storage%ymod, storage%dyda, status)

storage%old_dy = storage%dy
storage%dy = y - storage%ymod

if (status /= 0) return

j = 0
do l = 1, nv
  if (.not. storage%mask(l)) cycle
  j = j+1
  storage%wt = storage%dyda(:, l) * weight
  k = 0
  do m = 1, l
    k = k+1
    co_alpha(j, k) = dot_product(storage%wt, storage%dyda(:, m))
    co_alpha(k, j) = co_alpha(j, k)
  end do
  da_beta(j) = dot_product(storage%dy, storage%wt)
end do

chisq = dot_product(storage%dy**2, weight)

end subroutine super_mrqcof

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine super_gaussj (a, b, status)
! 
! Routine to solve a set of linear equations:
!   a * x = b
! Where a(N,N) is a square matrix and b(N,M) is rectangular with M being
! the number of linear equations to solve..
!
! This is the gaussj routine from Num Rec with an added status argument.
!
! Input:
!   a(:,:) -- Real(rp): matrix.
!   b(:,:) -- Real(rp): matrix.
!
! Output:
!   a(:,:) -- Real(rp): Inverse to input a(:,:) matrix.
!   b(:,:) -- Real(rp): Solutions to a * x = b euqations.
!   status -- Integer: Status. Set to -1 or -2 if there is an error.
!               Set to 0 otherwise.
!-

subroutine super_gaussj (a, b, status)

implicit none

real(rp), dimension(:, :), intent(inout) :: a, b
real(rp) :: pivinv

real(rp), dimension(size(a, 1)) :: dumc
integer, dimension(size(a, 1)) :: ipiv, indxr, indxc
logical, dimension(size(a, 1)) :: lpiv

integer status
integer, target :: irc(2)
integer :: i, l, n
integer, pointer :: irow, icol

character(*), parameter :: r_name = 'super_gaussj'

!

n = assert_equal([size(a, 1), size(a, 2), size(b, 1)], 'gaussj')
irow => irc(1)
icol => irc(2)
ipiv = 0

do i = 1, n
   lpiv = (ipiv == 0)
   irc = maxloc(abs(a), (spread(lpiv, dim = 2, ncopies = n) .and. spread(lpiv, dim = 1, ncopies = n)))
   ipiv(icol) = ipiv(icol)+1
   if (ipiv(icol) > 1) then
      status = -1
      call out_io (s_error$, r_name, 'SINGULAR MATRIX! (1)')
      return
   end if
   if (irow /= icol) then
      call swap(a(irow, :), a(icol, :))
      call swap(b(irow, :), b(icol, :))
   end if
   indxr(i) = irow
   indxc(i) = icol
   if (a(icol, icol) == 0.0) then
      status = -2
      call out_io (s_error$, r_name, 'SINGULAR MATRIX! (2)')
      return
   end if
   pivinv = 1.0_rp/a(icol, icol)
   a(icol, icol) = 1.0
   a(icol, :) = a(icol, :)*pivinv
   b(icol, :) = b(icol, :)*pivinv
   dumc = a(:, icol)
   a(:, icol) = 0.0
   a(icol, icol) = pivinv
   a(1:icol-1, :) = a(1:icol-1, :)-outer_product(dumc(1:icol-1), a(icol, :))
   b(1:icol-1, :) = b(1:icol-1, :)-outer_product(dumc(1:icol-1), b(icol, :))
   a(icol+1:, :) = a(icol+1:, :)-outer_product(dumc(icol+1:), a(icol, :))
   b(icol+1:, :) = b(icol+1:, :)-outer_product(dumc(icol+1:), b(icol, :))
end do

do l = n, 1, -1
   call swap(a(:, indxr(l)), a(:, indxc(l)))
end do

status = 0

end subroutine super_gaussj

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine super_ludcmp (a, indx, d, err)
!
! This routine is essentially ludcmp from Numerical Recipes with the added feature
! that an error flag is set instead of bombing the program when there is a problem.
!
! Input:
!   a(:,:) -- Real(rp): Input matrix.
!
! Output
!   indx(:) -- Integer: See NR.
!   d       -- Real(rp): See NR.
!   err     -- Logical: Error flag set True if there is a problem. False otherwise.
!-

subroutine super_ludcmp(a,indx,d, err)

implicit none
real(rp), dimension(:,:), intent(inout) :: a
integer, dimension(:), intent(out) :: indx
real(rp), intent(out) :: d
real(rp), dimension(size(a,1)) :: vv
real(rp), parameter :: tiny = 1.0e-20_rp
integer :: j,n,imax
character(*), parameter :: r_name = 'super_ludcmp'
logical err

!

err = .true.
n = assert_equal([size(a,1),size(a,2),size(indx)],'ludcmp')
d = 1.0
vv = maxval(abs(a),dim = 2)
if (any(vv == 0.0)) then
  call out_io (s_error$, r_name, 'singular matrix')
  return
endif
vv = 1.0_rp/vv
do j = 1,n
  imax = (j-1)+maxloc(vv(j:n)*abs(a(j:n,j)), 1)
  if (j /= imax) then
    call swap(a(imax,:),a(j,:))
    d = -d
    vv(imax) = vv(j)
  end if
  indx(j) = imax
  if (a(j,j) == 0.0) a(j,j) = tiny
  a(j+1:n,j) = a(j+1:n,j)/a(j,j)
  a(j+1:n,j+1:n) = a(j+1:n,j+1:n)-outer_product(a(j+1:n,j),a(j,j+1:n))
end do
err = .false.

end subroutine super_ludcmp

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Function super_qromb (func, a, b, rel_tol, abs_tol, k_order, err_flag) result (integral)
!
! Routine to use Romberg integration to integrate a function.
!
! This is essentially qromb from Numerical Recipes.
!
! The function func should satisfy the following interface:
!   function func(x) result (value)
!     real(rp), intent(in) :: x(:) ! Array of points to evaluate at.
!     real(rp) :: value(size(x))   ! Array of evaluation values.
!   end function func
!
! Input:
!   func    -- function: Function to integrate.
!   a       -- real(rp): Lower bound of integration region.
!   b       -- real(rp): Upper bound of integration region.
!   rel_tol -- real(rp): Relative tolerance.
!   abs_tol -- real(rp): Absolute tolerance.
!   k_order -- integer: Integration order. For smooth functions 5 is a good value. Use 2 if not smooth.
!
! Output:
!   integral  -- real(rp): Integral.
!   err_flag  -- logical: Set True if there is an error.
!-

function super_qromb (func, a, b, rel_tol, abs_tol, k_order, err_flag) result (integral)

implicit none

integer, parameter :: jmax = 50, jmaxp = jmax+1

real(rp)  :: a, b, rel_tol, abs_tol, integral
real(rp) :: d_int, h(jmaxp), s(jmaxp)

integer :: k_order
integer :: j, k, km

logical err_flag

interface
  function func(x)
  import
  real(rp), intent(in) :: x(:)
  real(rp) :: func(size(x))
  end function func
end interface

!

err_flag = .false.
k = k_order
km = k - 1
h(1) = 1.0

!

if (a == b) then
  integral = 0
  return
endif

!

do j = 1, jmax
  call super_trapzd(func, a, b, s(j), j)
  if (j >= k) then
    call super_polint(h(j-km:j), s(j-km:j), 0.0_rp, integral, d_int)
    if (abs(d_int) <= rel_tol*abs(integral) + abs_tol) return
  end if
  s(j+1) = s(j)
  h(j+1) = 0.25_rp*h(j)
end do

err_flag = .true.

end function super_qromb

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Function super_polint (xa, ya, x, y, dy)
!
! This is essentially polint from Numerical Recipes.
!
! Input:
!   xa(:), ya(:)  -- real(rp):
!   x             -- real(rp):
!
! Output:
!   y             -- real(rp):
!   dy            -- real(rp):
!-

subroutine super_polint(xa, ya, x, y, dy)

implicit none

real(rp), intent(in) :: xa(:), ya(:)
real(rp), intent(in) :: x
real(rp), intent(out) :: y, dy
integer :: m, n, ns
real(rp), dimension(size(xa)) :: c, d, den, ho

!

n = assert_equal([size(xa), size(ya)], 'polint')
c = ya
d = ya
ho = xa-x
ns = minloc(abs(x-xa), 1)
y = ya(ns)
ns = ns-1

do m = 1, n-1
  den(1:n-m) = ho(1:n-m)-ho(1+m:n)
  if (any(den(1:n-m) == 0.0)) call err_exit('polint: calculation failure')
  den(1:n-m) = (c(2:n-m+1)-d(1:n-m))/den(1:n-m)
  d(1:n-m) = ho(1+m:n)*den(1:n-m)
  c(1:n-m) = ho(1:n-m)*den(1:n-m)

  if (2*ns < n-m) then
    dy = c(ns+1)
  else
    dy = d(ns)
    ns = ns-1
  end if

  y = y+dy
end do

end subroutine super_polint

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine super_trapzd (func, a, b, s, n)
!
! This is essentially trapzd from Numerical Recipes.
!
! Input:
!   func
!   a, b      -- real(rp):
!   s         -- real(rp):
!   n         -- integer:
!
! Output:
!   s         -- rel(rp):
!-

subroutine super_trapzd(func, a, b, s, n)

implicit none

real(rp), intent(in) :: a, b
real(rp), intent(inout) :: s
integer, intent(in) :: n
real(rp) :: del, fsum, start, vec(2**max(0,(n-2)))
integer :: it

interface
  function func(x) result (value)
  import
  real(rp), intent(in) :: x(:)
  real(rp) :: value(size(x))
  end function func
end interface

!

if (n == 1) then
  s = 0.5_rp*(b-a)*sum(func( [a, b] ))
else
  it = 2**(n-2)
  del = (b-a)/it
  start = a - 0.5_rp*del
  do it = 1, size(vec)
    vec(it) = start + it * del
  enddo
  fsum = sum(func(vec))
  s = 0.5_rp*(s+del*fsum)
end if

end subroutine super_trapzd

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Function super_qromb_2D (func, ax, bx, ay, by, rel_tol, abs_tol, k_order, err_flag) result (integral)
!
! Routine to use Romberg integration to do a double integral of a function.
!
! This is essentially qromb from Numerical Recipes.
!
! The function func should satisfy the following interface:
!   function func(x,y) result (value)
!     real(rp) :: x, y      ! point to evaluate at.
!     real(rp) :: value     ! evaluation value.
!   end function func
!
! Input:
!   func    -- function: Function to integrate.
!   a       -- real(rp): Lower bound of integration region.
!   b       -- real(rp): Upper bound of integration region.
!   rel_tol -- real(rp): Relative tolerance.
!   abs_tol -- real(rp): Absolute tolerance.
!   k_order -- integer: Integration order. For smooth functions 5 is a good value. Use 2 if not smooth.
!
! Output:
!   integral  -- real(rp): Integral.
!   err_flag  -- logical: Set True if there is an error.
!-

function super_qromb_2D (func, ax, bx, ay, by, rel_tol, abs_tol, k_order, err_flag) result (integral)

implicit none

integer, parameter :: jmax = 50, jmaxp = jmax+1

real(rp)  :: ax, bx, ay, by, rel_tol, abs_tol, integral
real(rp) :: d_int, h(jmaxp), s(jmaxp)

integer :: k_order
integer :: j, k, km

logical err_flag

interface
  function func(x,y) result (r)
  import
  real(rp) :: x, y
  real(rp) :: r
  end function func
end interface

! This is the same a qromb except func is two dimensional.
! It is in trapzd_2D where the code had to be changed.

err_flag = .false.
k = k_order
km = k - 1
h(1) = 1.0

!

if (ax == bx .or. ay == by) then
  integral = 0
  return
endif

do j = 1, jmax
  call trapzd_2D(ax, bx, ay, by, s(j), j)
  if (j >= k) then
    call super_polint(h(j-km:j), s(j-km:j), 0.0_rp, integral, d_int)
    if (abs(d_int) <= rel_tol*abs(integral) + abs_tol) return
  end if
  s(j+1) = s(j)
  h(j+1) = 0.25_rp*h(j)
end do

err_flag = .true.

!----------------------------------------------------
contains
! This second integral is over x

subroutine trapzd_2D (ax, bx, ay, by, s, n_step)

implicit none
real(rp) :: ax, bx, ay, by, rel_tol, abs_tol
real(rp) :: s
real(rp) :: delx, dely, fsum
integer :: n_step, it, n, nx, ny

! Corner points get weight 1/4

if (n_step == 1) then
  s = 0.25_rp * (bx-ax) * (by-ay) * (func(ax,ay) + func(ax,by) + func(bx,ay) + func(bx,by))

else
  it=2**(n_step-1)
  delx = (bx-ax) / it
  dely = (by-ay) / it
  s =  0.25_rp * s

  ! Side points get weight 1/2
  do n = 1, it, 2
    s = s + 0.5_rp * delx * dely * &
          (func(ax+n*delx,ay) + func(ax+n*delx,by) + func(ax,ay+n*dely) + func(bx,ay+n*dely))
  enddo

  ! Interior points get weight 1
  do nx = 1, it-1, 2
    do ny = 1, it-1
      s = s + delx * dely * func(ax+nx*delx,ay+ny*dely)
    enddo
  enddo

  do nx = 2, it-2, 2
    do ny = 1, it-1, 2
      s = s + delx * dely * func(ax+nx*delx,ay+ny*dely)
    enddo
  enddo
end if

end subroutine trapzd_2D

end function super_qromb_2D

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Function super_poly (x, coef) result (value)
!
! Routine to compute Sum: coef(i)*x^i
!
! Input:
!   x       -- real(rp): Variable.
!   coef(:) -- real(rp): Coefficients.
!
! Output:
!   value   -- real(rp): Polynomial value.
!-

function super_poly (x, coeffs) result(value)

real(rp), intent(in) :: x
real(rp), intent(in) :: coeffs(:)
real(rp) :: value
integer :: i, n

!

n = size(coeffs)

if (n <= 0) then
  value = 0.0_rp

else
  value = coeffs(n)
  do i = n-1, 1, -1
    value = x*value + coeffs(i)
  end do
end if

end function super_poly

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine covar_expand(covar,maska)
!
! Routine to expand the covariance matrix covar so as to take into account parameters that have been masked.
! The rows and columns corresponding to these masked parameters will be set to zero.
! This is essentially covsrt from NR.
!
! Input:
!   covar(:,:)    -- real(rp): Covariance matrix.
!   maska(:)      -- real(rp): Mask vector.
!
! Output:
!   covar(:,:)    -- real(rp): Expanded matrix
!-

subroutine covar_expand(covar, maska)

implicit none

real(rp), intent(inout) :: covar(:, :)
logical, intent(in) :: maska(:)
integer :: ma, mfit, j, k

!

ma =  size(maska)
mfit = count(maska)
covar(mfit+1:ma, 1:ma) = 0.0
covar(1:ma, mfit+1:ma) = 0.0
k = mfit
do j = ma, 1, -1
  if (maska(j)) then
    call swap(covar(1:ma, k), covar(1:ma, j))
    call swap(covar(k, 1:ma), covar(j, 1:ma))
    k = k-1
  end if
end do

end subroutine covar_expand

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------



end module
