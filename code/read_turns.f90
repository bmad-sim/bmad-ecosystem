!........................................................................
!+
! Subroutine : read_turns()
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
! Revision 1.1  2005/06/14 14:59:02  cesrulib
! Initial revision
!
!
!........................................................................
!
#include "CESR_platform.inc"

subroutine read_turns()
  use precision_def
  implicit none

  real,allocatable :: vals(:,:,:)
  integer,allocatable :: n_out(:,:)
  integer :: n_tot,n_past,n_curr,n_cycle(1:50),i,j,inputstatus1,n_full,n_inc,i_cycle
  character(128) :: file_name

!  write(*,'(A)',advance="no") 'Please enter a .turns file name to be read: '
!  read '(A)', file_name
 
  file_name = 'beambeam.turns'
  open(unit=15,file=file_name,action="read",position="rewind")

! read off step size
  rewind(15)
  read(15,'(118x,i6)',iostat=inputstatus1) n_inc
  IF(inputstatus1 > 0) STOP "*****Input Error--n_inc Read*****"

! count total number of turns
  i = 0
  rewind(15)
  DO 
     i = i+1
     READ(15,'()', IOSTAT = inputstatus1)
     IF(inputstatus1 > 0) STOP "*****Input Error--n_tot Count*****"
     IF(inputstatus1 < 0) THEN
        n_tot = (i-1)*n_inc
        EXIT
     END IF
  END DO

! count number of turns in a cycle(n_cycle) & number of cycles (i_cycle)
  n_past = 0
  n_curr = 0
  n_cycle(:) = 0
  i_cycle = 1
  rewind(15)
     do
        n_past = n_curr
        read(15,'(118x,i6)',iostat=inputstatus1) n_curr
        if(inputstatus1 > 0) STOP "*****Input Error--n_cycle Count*****"
        if(inputstatus1 < 0) then
           n_cycle(i_cycle) = n_past
           exit
        end if
        if(n_curr.lt.n_past) then
           n_cycle(i_cycle) = n_past
           i_cycle = i_cycle + 1
        end if
     end do

     n_full = 0
     do i=1,i_cycle
        n_curr = n_cycle(i)
        if(n_curr.gt.n_full) n_full = n_curr
     end do

! read out n,luminocities,and n_out from the .turns file
     allocate(vals(1:n_full,1:i_cycle,1:3),n_out(1:n_full,1:i_cycle))

     rewind(15)
     do i=1,i_cycle
        do j=1,n_cycle(i)/n_inc
           read(15,'(126x,e12.4,2x,i5,2x,e12.4)',IOSTAT = inputstatus1) vals(j,i,1),n_out(j,i),vals(j,i,3)
           if(inputstatus1 > 0) STOP "*****Input Error--Read Cols*****"
        end do
     end do

     close(15)

  open(unit=16,file='beambeam.lums',action='write')     
  do j=1,n_full/n_inc
     write(16,'(1x,i6)',advance="no") j*n_inc
     do i=1,i_cycle
        if(vals(j,i,1)/=0) then
           write(16,'(1x,e12.4,1x,i5,1x,e12.4)',advance="no") vals(j,i,1),n_out(j,i),vals(j,i,3)
        else
           write(16,'(32x)',advance="no")
        end if
     end do
     write(16,*) ''
  end do
  close(16)

  
  deallocate(vals,n_out)
end subroutine read_turns
