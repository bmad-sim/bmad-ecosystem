!+
! Subroutine FitPoly(coe, x, y, order, samples)
!
! Subroutine to fit a polynomial, y = coe(0) + coe(1)*x + coe(2)*x^2 + ...,
! to the input data of x and y via least squares.
!
! Input:
!     x(*) -- Real: vector of sample 'x' data
!     y(*) -- Real: vector of sample 'y' data
!     order -- Integer: order of fitted polynomial
!     samples -- Integer: how many 'x, y' data samples
!
! Output:
!     coe(0:*) -- Real: array of polynomial coefficients
!-

subroutine fitpoly(coe, x, y, order, samples)
  implicit none
  integer maxcoe, maxsamp
  parameter(maxcoe=10, maxsamp=100)
  integer order, samples, numcoe
  real coe(0:*), x(*), y(*)
  real Xmat(maxsamp,maxcoe), XtX(maxcoe,maxcoe), Xty(maxcoe)
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

  return
  end
