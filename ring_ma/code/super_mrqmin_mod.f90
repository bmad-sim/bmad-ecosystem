! This is an adaptation of the MRQMIN routine from NR.
! It has been modified to trap singular matrices and...
module super_mrqmin_mod
  implicit none
  logical super_mrqmin_error

contains
  SUBROUTINE super_mrqmin(x,y,sig,a,maska,covar,alpha,chisq,funcs,alamda)
  use precision_def
  USE nrtype; USE nrutil, ONLY : assert_eq,diagmult
  USE nr, ONLY : covsrt
  IMPLICIT NONE
  REAL(rp), DIMENSION(:), INTENT(IN) :: x,y,sig
  REAL(rp), DIMENSION(:), INTENT(INOUT) :: a
  REAL(rp), DIMENSION(:,:), INTENT(OUT) :: covar,alpha
  REAL(rp), INTENT(OUT) :: chisq
  REAL(rp), INTENT(INOUT) :: alamda
  LOGICAL, DIMENSION(:), INTENT(IN) :: maska
  INTERFACE
     SUBROUTINE funcs(x,a,yfit,dyda)
       USE precision_def
       REAL(rp), DIMENSION(:), INTENT(IN) :: x,a
       REAL(rp), DIMENSION(:), INTENT(OUT) :: yfit
       REAL(rp), DIMENSION(:,:), INTENT(OUT) :: dyda
     END SUBROUTINE funcs
  END INTERFACE
  INTEGER(I4B) :: ma,ndata
  INTEGER(I4B), SAVE :: mfit


  super_mrqmin_error = .false.
  call mrqmin_private

CONTAINS
!BL
  SUBROUTINE mrqmin_private
    REAL(rp), SAVE :: ochisq
    REAL(rp), DIMENSION(:), ALLOCATABLE, SAVE :: atry,beta
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, SAVE :: da
    ndata=assert_eq(size(x),size(y),size(sig),'mrqmin: ndata')
    ma=assert_eq((/size(a),size(maska),size(covar,1),size(covar,2),&
         size(alpha,1),size(alpha,2)/),'mrqmin: ma')
    mfit=count(maska)
    if (alamda < 0.0) then
       allocate(atry(ma),beta(ma),da(ma,1))
       alamda=0.001_rp
       call mrqcof(a,alpha,beta)
       ochisq=chisq
       atry=a
    end if
    covar(1:mfit,1:mfit)=alpha(1:mfit,1:mfit)
    call diagmult(covar(1:mfit,1:mfit),1.0_rp+alamda)
    da(1:mfit,1)=beta(1:mfit)
    call gaussj(covar(1:mfit,1:mfit),da(1:mfit,1:1))

    if (alamda == 0.0) then
       call covsrt(covar,maska)
       call covsrt(alpha,maska)
       deallocate(atry,beta,da)
       RETURN
    end if

! Here check to see if GAUSSJ worked.  If not, move on.
    if (super_mrqmin_error) return

    atry=a+unpack(da(1:mfit,1),maska,0.0_rp)
    call mrqcof(atry,covar,da(1:mfit,1))
    if (chisq < ochisq) then
       alamda=0.1_rp*alamda
       ochisq=chisq
       alpha(1:mfit,1:mfit)=covar(1:mfit,1:mfit)
       beta(1:mfit)=da(1:mfit,1)
       a=atry
    else
       alamda=10.0_rp*alamda
       chisq=ochisq
    end if
  END SUBROUTINE mrqmin_private
!BL
  SUBROUTINE mrqcof(a,alpha,beta)
    REAL(rp), DIMENSION(:), INTENT(IN) :: a
    REAL(rp), DIMENSION(:), INTENT(OUT) :: beta
    REAL(rp), DIMENSION(:,:), INTENT(OUT) :: alpha
    INTEGER(I4B) :: j,k,l,m
    REAL(rp), DIMENSION(size(x),size(a)) :: dyda
    REAL(rp), DIMENSION(size(x)) :: dy,sig2i,wt,ymod
    call funcs(x,a,ymod,dyda)
    sig2i=1.0_rp/(sig**2)
    dy=y-ymod
    j=0
    do l=1,ma
       if (maska(l)) then
          j=j+1
          wt=dyda(:,l)*sig2i
          k=0
          do m=1,l
             if (maska(m)) then
                k=k+1
                alpha(j,k)=dot_product(wt,dyda(:,m))
                alpha(k,j)=alpha(j,k)
             end if
          end do
          beta(j)=dot_product(dy,wt)
       end if
    end do
    chisq=dot_product(dy**2,sig2i)
  END SUBROUTINE mrqcof

  SUBROUTINE gaussj(a,b)
    USE nrtype; USE nrutil, ONLY : assert_eq,nrerror,outerand,outerprod,swap
    IMPLICIT NONE
    REAL(rp), DIMENSION(:,:), INTENT(INOUT) :: a,b
    INTEGER(I4B), DIMENSION(size(a,1)) :: ipiv,indxr,indxc
    LOGICAL(LGT), DIMENSION(size(a,1)) :: lpiv
    REAL(rp) :: pivinv
    REAL(rp), DIMENSION(size(a,1)) :: dumc
    INTEGER(I4B), TARGET :: irc(2)
    INTEGER(I4B) :: i,l,n
    INTEGER(I4B), POINTER :: irow,icol
    n=assert_eq(size(a,1),size(a,2),size(b,1),'gaussj')
    irow => irc(1)
    icol => irc(2)
    ipiv=0
    do i=1,n
       lpiv = (ipiv == 0)
       irc=maxloc(abs(a),outerand(lpiv,lpiv))
       ipiv(icol)=ipiv(icol)+1
       if (ipiv(icol) > 1) then
          super_mrqmin_error = .true.
          return
       end if
       if (irow /= icol) then
          call swap(a(irow,:),a(icol,:))
          call swap(b(irow,:),b(icol,:))
       end if
       indxr(i)=irow
       indxc(i)=icol
       if (a(icol,icol) == 0.0) then
          super_mrqmin_error = .true.
          return
       end if
       pivinv=1.0_rp/a(icol,icol)
       a(icol,icol)=1.0
       a(icol,:)=a(icol,:)*pivinv
       b(icol,:)=b(icol,:)*pivinv
       dumc=a(:,icol)
       a(:,icol)=0.0
       a(icol,icol)=pivinv
       a(1:icol-1,:)=a(1:icol-1,:)-outerprod(dumc(1:icol-1),a(icol,:))
       b(1:icol-1,:)=b(1:icol-1,:)-outerprod(dumc(1:icol-1),b(icol,:))
       a(icol+1:,:)=a(icol+1:,:)-outerprod(dumc(icol+1:),a(icol,:))
       b(icol+1:,:)=b(icol+1:,:)-outerprod(dumc(icol+1:),b(icol,:))
    end do
    do l=n,1,-1
       call swap(a(:,indxr(l)),a(:,indxc(l)))
    end do
  END SUBROUTINE gaussj
END SUBROUTINE super_mrqmin
end module super_mrqmin_mod
