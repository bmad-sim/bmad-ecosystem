!........................................................................
!+
! Subroutine subroutine histogram (ele, coord, in_file, sig, a_out)
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
! Revision 1.3  2007/01/30 16:14:31  dcs
! merged with branch_bmad_1.
!
! Revision 1.2.2.1  2006/12/22 20:30:42  dcs
! conversion compiles.
!
! Revision 1.2  2005/09/21 20:59:07  dcs
! more changes to get around compiler bug.
!
! Revision 1.1.1.1  2005/06/14 14:59:02  cesrulib
! Beam Simulation Code
!
!
!........................................................................
!
subroutine histogram (ele, coord, in_file, sig, a_out)
  
  use bmad_interface
  
  implicit none
  
  type (coord_struct) coord(:)
  type (ele_struct) ele
  real(rp) sum(6), avg(6), rms(6)
  real(rp) disp_sum(6), disp_avg(6), disp_rms(6)
  real(rp) x, disp_x
  real(rp) a(6),a_out(1:3),y,f,fsum,ysum
  real(rp) disp_a(6),disp_fsum,disp_ysum
  real(rp), dimension(size(coord),3) :: amp, disp_amp
  real(rp) V_inv(4,4), temp_vec(4)
  real(rp) sig(3)

  integer i,j,k, n, disp_n
  integer number(0:60,3),disp_number(-60:60,3)
  integer nbin, nsig
  integer psize
  integer unit
  
  character(*) in_file
  
  nbin=10
  nsig=3
  
  do j=1,3
     do i=-nbin*nsig,nbin*nsig
        if(i>=0)number(i,j)=0
        disp_number(i,j)=0
     end do
  end do
  
!         i=0
!         do while(.true.)
!           i=i+1
!           read(12,*,end=99)(vec(i,j),j=1,6)
!         end do
! 99      continue

  psize =size(coord)-1
  
  V_inv(1:4,1:4)=0
  forall (j=1:4) V_inv(j,j)=ele%gamma_c
  V_inv(1:2,3:4) = -ele%c_mat
  V_inv(3,1) = ele%c_mat(2,2)
  V_inv(4,2) = ele%c_mat(1,1)
  V_inv(3,2) = -ele%c_mat(1,2)
  V_inv(4,1) = -ele%c_mat(2,1)

! amplitude
  do j=1,psize
     
     temp_vec(1:4) = matmul(V_inv, coord(j)%vec(1:4)) !switch back to normal modes
     
     amp(j,1) = ele%a%gamma* temp_vec(1)**2 &
          + 2*ele%a%alpha*temp_vec(1)*temp_vec(2) &
          + ele%a%beta * temp_vec(2)**2
     amp(j,2) = ele%b%gamma* temp_vec(3)**2 &
          + 2*ele%b%alpha*temp_vec(3)*temp_vec(4) &
          + ele%b%beta * temp_vec(4)**2
     amp(j,3) = ele%z%gamma* coord(j)%vec(5)**2 &
          + 2*ele%z%alpha*coord(j)%vec(5)*coord(j)%vec(6) &
          + ele%z%beta * coord(j)%vec(6)**2
     
     
     disp_amp(j,1)=temp_vec(1)
     disp_amp(j,2)=temp_vec(3)
     disp_amp(j,3)=coord(j)%vec(5)
     
  end do
  
  do k = 1,3
     sum(k) = 0.
     disp_sum(k) =0.
     do j=1,psize
        sum(k) = sum(k)+amp(j,k)
        disp_sum(k) = disp_sum(k)+disp_amp(j,k)
     end do
     avg(k) = sum(k)/float(psize)
     disp_avg(k) = disp_sum(k)/float(psize)
     avg(k) = 0.
  end do
  
  do k=1,3
     sum(k)=0.
     disp_sum(k)=0.
     do j=1,psize
        sum(k) = sum(k)+(amp(j,k)-avg(k))**2
        disp_sum(k) = disp_sum(k)+(disp_amp(j,k)-disp_avg(k))**2
     end do
     rms(k)=sqrt(sum(k)/float(psize))
     disp_rms(k)=sqrt(disp_sum(k)/float(psize))
  end do

  do k=1,3
     do j=1,psize
        x=amp(j,k)/rms(k)*nbin
        disp_x=(disp_amp(j,k)-disp_avg(k))/disp_rms(k)*nbin
        if(disp_x .gt. 0)disp_n = ((disp_amp(j,k)-disp_avg(k))/disp_rms(k)*nbin) +.5
        if(disp_x .lt. 0)disp_n = ((disp_amp(j,k)-disp_avg(k))/disp_rms(k)*nbin) -.5
        n = (amp(j,k)/rms(k)*nbin) +.5
3       format(1x,a24,e12.4,a7,i4)
        if(n .le. nsig*nbin)then
           number(n,k) = number(n,k) + 1
        endif
        if(abs(disp_n) .le. nsig*nbin)then
           disp_number(disp_n,k) = disp_number(disp_n,k) + 1
        endif
     end do
  end do

  do k=1,3
     fsum = 0.
     ysum = 0.
     disp_fsum=0.
     disp_ysum=0.
     do j = 0, nsig*nbin
        if(j >= 0)then             
           f= exp(-(float(j)*rms(k)/nbin-avg(k))**2/2./(rms(k))**2)
           y = float(number(j,k))
           fsum = fsum + f*f
           ysum = ysum + y*f
        endif
        
        f= exp(-(float(j)*disp_rms(k)/nbin-disp_avg(k))**2/2./(disp_rms(k))**2)
        y = float(disp_number(j,k))
        disp_fsum = disp_fsum + f*f
        disp_ysum = disp_ysum + y*f
!             write(6,*)'   j,f,y,rms(k)', j,f,y,rms(k)
     end do
     A(k) = ysum/fsum
     disp_A(k) = disp_ysum/disp_fsum
  end do

!           write(6,4)a(k),avg(k),rms(k)

  if(index(in_file,'junk') == 0)then
     open (unit=2, file = in_file)
     write(2,'(a8,a1,e12.4,a1)') '   Ax = ',"`",a(1),"'"
     write(2,'(a13,a1,e12.4,a1)')'   disp_Ax = ',"`",disp_a(1),"'"
     write(2,'(a15,a1,e12.4,a1)')'  disp_avg_x = ',"`",disp_avg(1),"'"
     write(2,'(a10,a1,e12.4,a1)')'  sig_x = ',"`",rms(1),"'"
     write(2,'(a15,a1,e12.4,a1)')'  disp_sig_x = ',"`",disp_rms(1),"'"
     
     write(2,'(a8,a1,e12.4,a1)') '   Ay = ',"`",a(2),"'"
     write(2,'(a13,a1,e12.4,a1)')'   disp_Ay = ',"`",disp_a(2),"'"
     write(2,'(a15,a1,e12.4,a1)')'  disp_avg_y = ',"`",disp_avg(2),"'"
     write(2,'(a10,a1,e12.4,a1)')'  sig_y = ',"`",rms(2),"'"
     write(2,'(a15,a1,e12.4,a1)')'  disp_sig_y = ',"`",disp_rms(2),"'"
           
     write(2,'(a8,a1,e12.4,a1)') '   Az = ',"`",a(3),"'"
     write(2,'(a13,a1,e12.4,a1)')'   disp_Az = ',"`",disp_a(3),"'"
     write(2,'(a15,a1,e12.4,a1)')'  disp_avg_z = ',"`",disp_avg(3),"'"
     write(2,'(a10,a1,e12.4,a1)')'  sig_z = ',"`",rms(3),"'"
     write(2,'(a15,a1,e12.4,a1)')'  disp_sig_z = ',"`",disp_rms(3),"'"


     write(2,*) 'return'
     write(2,*)


     do j = -nsig*nbin, nsig*nbin
!           if(j>=0)write(2,2)(j*rms(k)/nbin, number(j,k),k=1,3)
        write(2,2)(j*disp_rms(k)/nbin, disp_number(j,k),k=1,3)
2       format(1x,3(e12.4,i4))
     end do
  endif

  sig(1:3) = disp_rms(1:3)
  a_out(1:3) = A(1:3)
  
end subroutine histogram


