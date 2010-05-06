!........................................................................
!+
! Subroutine : query_real( parameter, default, fmt)
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
! $Id: query.f90,v 1.1.1.1 2005/06/14 14:59:02 cesrulib Exp $
!
! $Log: query.f90,v $
! Revision 1.1.1.1  2005/06/14 14:59:02  cesrulib
! Beam Simulation Code
!
!
!........................................................................
!
#include "CESR_platform.h"

subroutine query_real( parameter, default, fmt)

  use bmadz_mod

  implicit none
  character*(*) parameter, fmt
  character*10 string
  character*30 fm
  character prompt*10/' (default='/, close*3/') ?'/
  real(rdef) default

   fm = '(1x,a,a10,'//fmt//',a3,$)'
!  fm = a,' (default=',fmt,') ?',$)
  print fm,parameter,prompt,default,close
  read(5, '(a)')string
  if(string(1:10) == '          ' )return
  read(string,*)default

  return
end

subroutine query_int( parameter, default, fmt)

  implicit none
  character*(*) parameter, fmt
  character*10 string
  character*30 fm
  character prompt*10/' (default='/, close*3/') ?'/
  integer default

   fm = '(1x,a,a10,'//fmt//',a3,$)'
!  fm = a,' (default=',fmt,') ?',$)
  print fm,parameter,prompt,default,close
  read(5, '(a)')string
  if(string(1:10) == '          ' )return
  read(string,*)default

  return
end


subroutine query_character( parameter, default, fmt)

  implicit none
  character*(*) parameter, fmt
  character*72 string
  character*30 fm
  character prompt*10/' (default='/, close*3/') ?'/
  character*(*) default

   fm = '(1x,a,a10,'//fmt//',a3,$)'
!  fm = a,' (default=',fmt,') ?',$)
  print fm,parameter,prompt,default,close
  read(5, '(a)')string
  if(string(1:10) == '          ' )return
  default = string

  return
end
