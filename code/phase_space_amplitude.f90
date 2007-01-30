!........................................................................
!+
! Subroutine : phase_space_amplitude (ele, orb, amp_x, amp_y)
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

subroutine phase_space_amplitude (ele, orb, amp_x, amp_y)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ele_struct) ele
  type (coord_struct) orb

  real(rdef) amp_x, amp_y
  real(rdef) V_inv(4,4), temp_vec(4)

  integer j

  V_inv(1:4,1:4)=0
  forall (j=1:4) V_inv(j,j)=ele%gamma_c
  V_inv(1:2,3:4) = -ele%c_mat
  V_inv(3,1) = ele%c_mat(2,2)
  V_inv(4,2) = ele%c_mat(1,1)
  V_inv(3,2) = -ele%c_mat(1,2)
  V_inv(4,1) = -ele%c_mat(2,1)
  temp_vec(1:4) = matmul(V_inv,orb%vec(1:4))  !switch lab coordinates
  orb%vec(1:4) = temp_vec(1:4)

  amp_x = sqrt(ele%a%gamma * orb%vec(1)**2 &
                + 2*ele%a%alpha * orb%vec(1) * orb%vec(2) &
                + ele%a%beta * orb%vec(2)**2)

  amp_y = sqrt(ele%b%gamma * orb%vec(3)**2 &
                + 2*ele%b%alpha * orb%vec(3) * orb%vec(4) &
                + ele%b%beta * orb%vec(4)**2)

  return

end subroutine
