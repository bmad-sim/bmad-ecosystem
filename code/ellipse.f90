!........................................................................
!+
! Subroutine ellipse(ele, major,minor,theta)
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
 subroutine ellipse(ele, major,minor,theta)

  use bmad_struct
  use bmad_interface

      implicit none

      type(ele_struct) ele

      real(rdef) bet, rootb, alp, gamma, A,B, q, r, s
      real(rdef) theta, s2, c2, sc,  major, minor

      integer ix,n
 
      bet = ele%x%beta
      rootb = sqrt(ele%x%beta)
      alp  = ele%x%alpha
      gamma = ele%gamma_c
      A= (-ele%c_mat(2,2)+ele%c_mat(1,2)*alp/bet)/gamma
      B = ele%c_mat(1,2)/bet/gamma
      q=1./B**2
      r = -2.*A/B**2
      s = (A**2+B**2)/B**2
      theta = 0.5*atan2(r, q-s)
      s2 = sin(theta)**2
      c2 = cos(theta)**2
      sc = sin(theta)*cos(theta)
      major = 1./sqrt(q*s2 - r*sc + s*c2)
      minor = 1./sqrt(q*c2 + r*sc + s*s2)

      return
      
  end
