
subroutine fft_plane(n,co, Q)

use bmad
 use nr
 use nrutil, only: dpc

implicit none

type(coord_struct) co(:)

integer n
integer i, imx(1), j, ixy, turn(8192)
integer isign

 real(rp) Q(3)
 real(rp) pie, range
 real(rp) ftx(8192), fty(8192), ftz(8192)
 real(rp), allocatable :: ra(:,:),a1(:)
 real(rp) d,b,c,A,tune(10,3)

 complex(dpc) cdata(1,8192)

 pie=atan(1.)*4

 allocate(ra(1:n,1:3))
 allocate(a1(1:n-10))


    range = 390.1
    do j = 1,n

     cdata(1,j) = cmplx(co(j)%vec(1), co(j)%vec(2))*(2*(sin((pie*j)/n))*(sin((pie*j)/n)))
!     ftx(j) = co(j)%vec(1)*(2*(sin((pie*j)/n))*(sin((pie*j)/n)))
    end do
    call fourrow(cdata(:,1:n),1)

    forall(i=1:n) ftx(i)=sqrt(cdata(1,i)* conjg(cdata(1,i)))
!    print '(11e12.4)',co(n-10:n)%vec(1)
!    call realft(ftx(1:n),1)

    isign=1
    do j = 1,n
     cdata(1,j) = cmplx(co(j)%vec(3),co(j)%vec(4)) *(2*(sin((pie*j)/n))*(sin((pie*j)/n)))
!     fty(j) = co(j)%vec(3) *(2*(sin((pie*j)/n))*(sin((pie*j)/n)))
    end do
    call fourrow(cdata(:,1:n), isign)

     forall(i=1:n)fty(i)=sqrt(cdata(1,i)*conjg(cdata(1,i)))
!    call realft(fty(1:n),1)
    isign=1

    do j = 1,n
     cdata(1,j) = cmplx(co(j)%vec(5),-co(j)%vec(6))*(2*(sin((pie*j)/n))*(sin((pie*j)/n)))
    end do
    call fourrow(cdata(:,1:n), isign)
    forall(i=1:n)ftz(i)=sqrt(cdata(1,i)*conjg(cdata(1,i)))

  do i = 1,n
    turn(i) = i
    ra(i,1) = abs(ftx(i))
    ra(i,2) = abs(fty(i))
    ra(i,3) = abs(ftz(i))
  end do


  do ixy = 1,3
   do i = 1,10

   imx = maxloc(ra(5:n-5,ixy)) + 4

    d= ra(imx(1),ixy)
    b=ra(imx(1)+1,ixy)
    if(d /= 0. .or. b /= 0.)then
    c=cos(twopi/n)
    A= (-(d+b*c)*(d-b)+b*sqrt(c*c*(d+b)*(d+b)-2*d*b*(2*c*c-c-1)))/(d*d+b*b+2*d*b*c)
    endif
    tune(i,ixy)= imx(1)/float(n)+(1/twopi)*asin(A*sin(twopi/float(n)))
!    print '(1x,3i5,2e12.4)',i,ixy,imx(1),ra(imx(1),ixy), tune(i,ixy)
    do j = max(1,imx(1)-20),min(n,imx(1)+20)
     ra(j,ixy)=0.
    end do
   end do
  end do

!  print '(1x,2a17,a19)','  vertical(kHz)  ','  horizontal(kHz) ',' longitudinal(kHz) ' 
!  do i = 1,10
!    print '(1x,1(3x,3f10.3,4x))', (tune(i,2))*range,(tune(i,1))*range,(tune(i,3))*range
!  end do

!  Q = 1.-0.5*tune(1,1:3)
  Q(1:3) = tune(1,1:3)
!  print *,' Q=',Q

  deallocate(ra)
  deallocate(a1)

 return
end
