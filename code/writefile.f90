!........................................................................
!+
! Subroutine : writefile(in_file, parameters)
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
! Revision 1.2  2005/09/21 20:59:07  dcs
! more changes to get around compiler bug.
!
! Revision 1.1.1.1  2005/06/14 14:59:02  cesrulib
! Beam Simulation Code
!
!
!........................................................................
!
subroutine writefile(in_file, parameters)
  
  use bmad_interface
  
  implicit none
  real(rp), dimension(1:,1:), intent(in) ::  parameters
  character(*), intent(in) ::  in_file
  


  if(index(in_file,'junk') == 0)then
     open (unit=2, file = in_file)
     write(2,'(a8,a1,e12.4,a1)') '   Ax = ',"`",parameters(1,1),"'"
!??     write(2,'(a13,a1,e12.4,a1)')'   disp_Ax = ',"`",disp_a(1),"'"
     write(2,'(a15,a1,e12.4,a1)')'  disp_avg_x = ',"`",parameters(2,1),"'"
     write(2,'(a10,a1,e12.4,a1)')'  sig_x = ',"`",parameters(3,1),"'"
!??     write(2,'(a15,a1,e12.4,a1)')'  disp_sig_x = ',"`",disp_rms(1),"'"
     
     write(2,'(a8,a1,e12.4,a1)') '   Ay = ',"`",parameters(1,2),"'"
!??     write(2,'(a13,a1,e12.4,a1)')'   disp_Ay = ',"`",disp_a(2),"'"
     write(2,'(a15,a1,e12.4,a1)')'  disp_avg_y = ',"`",parameters(2,2),"'"
     write(2,'(a10,a1,e12.4,a1)')'  sig_y = ',"`",parameters(3,2),"'"
!??     write(2,'(a15,a1,e12.4,a1)')'  disp_sig_y = ',"`",disp_rms(2),"'"
           
     write(2,'(a8,a1,e12.4,a1)') '   Az = ',"`",parameters(3,2),"'"
!??     write(2,'(a13,a1,e12.4,a1)')'   disp_Az = ',"`",disp_a(3),"'"
     write(2,'(a15,a1,e12.4,a1)')'  disp_avg_z = ',"`",parameters(2,3),"'"
     write(2,'(a10,a1,e12.4,a1)')'  sig_z = ',"`",parameters(3,3),"'"
!??     write(2,'(a15,a1,e12.4,a1)')'  disp_sig_z = ',"`",disp_rms(3),"'"


     write(2,*) 'return'
     write(2,*)


!     do j = -nsig*nbin, nsig*nbin
!           if(j>=0)write(2,2)(j*rms(k)/nbin, number(j,k),k=1,3)
!        write(2,2)(j*disp_rms(k)/nbin, disp_number(j,k),k=1,3)
!2       format(1x,3(e12.4,i4))
!     end do
  endif

end subroutine writefile
