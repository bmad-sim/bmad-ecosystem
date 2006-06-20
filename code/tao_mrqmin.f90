subroutine tao_mrqmin(x,y,sig,a,maska,covar,alpha,chisq,alamda, limited)
use nrtype; use nrutil, only : assert_eq,diagmult
use nr, only : covsrt,gaussj
use tao_interface, only: tao_mrq_func
use precision_def

implicit none

real(rp), dimension(:), intent(in) :: x,y,sig
real(rp), dimension(:), intent(inout) :: a
real(rp), dimension(:,:), intent(out) :: covar,alpha
real(rp), intent(out) :: chisq
real(rp), intent(inout) :: alamda
integer(i4b) :: ma,ndata
integer(i4b), save :: mfit
logical(lgt), dimension(:), intent(in) :: maska
logical limited

!

call tao_mrqmin_private

!----------------------------------------------
contains

subroutine tao_mrqmin_private

real(rp), save :: ochisq
real(rp), dimension(:), allocatable, save :: atry,beta
real(rp), dimension(:,:), allocatable, save :: da

!

ndata=assert_eq(size(x),size(y),size(sig),'mrqmin: ndata')
ma=assert_eq((/size(a),size(maska),size(covar,1),size(covar,2), &
              size(alpha,1),size(alpha,2)/),'mrqmin: ma')
mfit=count(maska)

if (alamda < 0.0) then
  allocate(atry(ma),beta(ma),da(ma,1))
  alamda=0.001_rp
  call tao_mrqcof(a,alpha,beta)
  if (limited) then
    deallocate(atry,beta,da)
    return
  endif
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
  return
end if

atry=a+unpack(da(1:mfit,1),maska,0.0_rp)
call tao_mrqcof(atry,covar,da(1:mfit,1))
if (limited) return

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

end subroutine tao_mrqmin_private

!-----------------------------------------------------------

subroutine tao_mrqcof(a,alpha,beta)

real(rp), dimension(:), intent(in) :: a
real(rp), dimension(:), intent(out) :: beta
real(rp), dimension(:,:), intent(out) :: alpha
integer(i4b) :: j,k,l,m
real(rp), dimension(size(x),size(a)) :: dyda
real(rp), dimension(size(x)) :: dy,sig2i,wt,ymod

!

call tao_mrq_func(x,a,ymod,dyda, limited)
if (limited) return

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

end subroutine tao_mrqcof

end subroutine tao_mrqmin