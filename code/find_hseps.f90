!........................................................................
!+
!
! Subroutine  subroutine find_hseps(ring, ix)
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
#include "CESR_platform.inc"

  subroutine find_hseps(ring, ix)

  use bmadz_mod

  implicit none

  type (lat_struct) ring

  integer ix(4),i

  ix(1:4) = 0

  i=0
   do while( any(ix(1:4) == 0))
   if(ring%ele(i)%name == 'H_SEP_08W')ix(1)=i
   if(ring%ele(i)%name == 'H_SEP_08E')ix(4)=i
   if(ring%ele(i)%name == 'H_SEP_45W')ix(2)=i
   if(ring%ele(i)%name == 'H_SEP_45E')ix(3)=i
   i=i+1
   if( i > ring%n_ele_max) then
    print *,' FIND_HSEPS: cannot find horizontal separators'
    stop
   endif
  enddo
  return
  end
