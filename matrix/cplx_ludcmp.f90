SUBROUTINE cplx_ludcmp(a,indx,d)

  USE sim_utils, dummy => cplx_ludcmp

  IMPLICIT NONE

  COMPLEX(rp), DIMENSION(:,:), INTENT(INOUT) :: a
  INTEGER, DIMENSION(:), INTENT(OUT) :: indx
  COMPLEX(rp), DIMENSION(size(a,1),size(a,1)) :: outerproduct
  REAL(rp), INTENT(OUT) :: d
  INTEGER :: k,n
  INTEGER i

  n=assert_equal([size(a,1),size(a,2),size(indx)],'cplx_ludcmp')

  d=1.0

  !- The pivoting approach used by nr does not work when the columns are complex
  !- conjugate pairs.  Set indx to 1,2,3,....
  do i=1,n
    indx(i) = i
  enddo

  do k=1,n-1
    a(k+1:n,k)=a(k+1:n,k)/a(k,k)
    outerproduct=spread(a(k+1:n,k),dim=2,ncopies=n-k)*spread(a(k,k+1:n),dim=1,ncopies=n-k)
    a(k+1:n,k+1:n)=a(k+1:n,k+1:n)-outerproduct
  end do

END SUBROUTINE cplx_ludcmp
