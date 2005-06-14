!........................................................................
!+
!
! Subroutine  subroutine find_hseps(ring, ix_)
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
#include "CESR_platform.h"

  subroutine find_hseps(ring, ix_)

  use bmadz_mod

  implicit none

  type (ring_struct) ring

  integer ix_(4),i

  ix_(1:4) = 0

  i=0
   do while( any(ix_(1:4) == 0))
   if(ring%ele_(i)%name == 'H_SEP_08W')ix_(1)=i
   if(ring%ele_(i)%name == 'H_SEP_08E')ix_(4)=i
   if(ring%ele_(i)%name == 'H_SEP_45W')ix_(2)=i
   if(ring%ele_(i)%name == 'H_SEP_45E')ix_(3)=i
   i=i+1
   if( i > ring%n_ele_max) then
    type *,' FIND_HSEPS: cannot find horizontal separators'
    stop
   endif
  enddo
  return
  end
