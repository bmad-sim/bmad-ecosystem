!........................................................................
!+
! Subroutine histogram_new (ele, coord, in_file, sig)
!
! Description:
!
! Arguments  :
!
! Mod/Commons:
!
! Calls      :
!
! Author     :
!
! Modified   :
!-
!........................................................................
!
! $Id$
!
! $Log$
! Revision 1.2  2007/01/30 16:14:31  dcs
! merged with branch_bmad_1.
!
! Revision 1.1.1.1.2.1  2006/12/22 20:30:42  dcs
! conversion compiles.
!
! Revision 1.1.1.1  2005/06/14 14:59:02  cesrulib
! Beam Simulation Code
!
!
!........................................................................
!
#include "CESR_platform.h"

 subroutine histogram_new (ele, coord, in_file, sig)

  use bmad_struct
  use bmad_interface

 implicit none

  type (coord_struct) coord(:)
  type (ele_struct) ele
         real(rp) sum(6), avg(6), rms(6)
         real(rp) x
         real(rp) a(6),y,f,fsum,ysum
         real(rp), dimension(size(coord),3) :: amp
         real(rp) sig(3)

         integer i,j,k, n
         integer number(0:60,3)
         integer nbin, nsig
         integer psize
         integer unit

         character*60 in_file

         nbin=10
         nsig=3

         do j=1,3
           do i=0,nbin*nsig
            number(i,j)=0
           end do
         end do

!         i=0
!         do while(.true.)
!           i=i+1
!           read(12,*,end=99)(vec(i,j),j=1,6)
!         end do
! 99      continue

         psize =size(coord)

! amplitude
         do j=1,psize
           amp(j,1) = ele%a%gamma* coord(j)%vec(1)**2 &
                    - 2*ele%a%alpha*coord(j)%vec(1)*coord(j)%vec(2) &
                    + ele%a%beta * coord(j)%vec(2)**2
           amp(j,2) = ele%b%gamma* coord(j)%vec(3)**2 &
                    - 2*ele%b%alpha*coord(j)%vec(3)*coord(j)%vec(4) &
                    + ele%b%beta * coord(j)%vec(4)**2
           amp(j,3) = ele%z%gamma* coord(j)%vec(5)**2 &
                    - 2*ele%z%alpha*coord(j)%vec(5)*coord(j)%vec(6) &
                    + ele%z%beta * coord(j)%vec(6)**2
         end do

         do k = 1,3
!          sum(k) = 0.
!          do j=1,psize
!            sum(k) = sum(k)+amp(j,k)
!           end do
!          avg(k) = sum(k)/float(psize)
          avg(k) = 0.
         end do

         do k=1,3
          sum(k)=0.
          do j=1,psize
             sum(k) = sum(k)+(amp(j,k)-avg(k))**2
           end do
          rms(k)=sqrt(sum(k)/float(psize))
        end do

        do k=1,3
         do j=1,psize
          x=amp(j,k)/rms(k)*nbin
!          if(x .gt. 0)n = (amp(j,k)/rms(k)*nbin) +.5
!          if(x .lt. 0)n = (amp(j,k)/rms(k)*nbin) -.5
           n = (amp(j,k)/rms(k)*nbin) +.5
 3        format(1x,a24,e12.4,a7,i4)
         if(n .le. nsig*nbin)then
            number(n,k) = number(n,k) + 1
          endif
          end do 
        end do       

         do k=1,3
           fsum = 0.
           ysum = 0.
           do j = 0, nsig*nbin
              f= exp(-(float(j)*rms(k)/nbin-avg(k))**2/2./(rms(k))**2)
              y = float(number(j,k))
              fsum = fsum + f*f
              ysum = ysum + y*f
!             write(6,*)'   j,f,y,rms(k)', j,f,y,rms(k)
           end do
           A(k) = ysum/fsum
         end do

!           write(6,4)a(k),avg(k),rms(k)

           open (unit=2, file = in_file, carriagecontrol='list')
           write(2,'(a8,a1,e12.4,a1)') '   Ax = ',"`",a(1),"'"
           write(2,'(a11,a1,e12.4,a1)') '   sig_x = ',"`",rms(1),"'"
           write(2,'(a8,a1,e12.4,a1)') '   Ay = ',"`",a(2),"'"
           write(2,'(a11,a1,e12.4,a1)') '   sig_y = ',"`",rms(2),"'"
           write(2,'(a8,a1,e12.4,a1)') '   Az = ',"`",a(3),"'"
           write(2,'(a11,a1,e12.4,a1)') '   sig_z = ',"`",rms(3),"'"
           write(2,*) 'return'
           write(2,*)

         do j = 0, nsig*nbin
           write(2,2)(j*rms(k)/nbin, number(j,k),k=1,3)
 2         format(1x,3(e12.4,i4))
         end do
       
         sig(1:3) = rms(1:3)

         end
