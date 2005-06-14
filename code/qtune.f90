!........................................................................
!+
! Subroutine : q_tune (ring, Q_x, Q_y, ok)
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

  subroutine q_tune (ring, Q_x, Q_y, ok)

!  set fractional tunes to Q_x, Q_y

  use bmad_struct
  use bmadz_interface
  use bsim_interface

  implicit none

  type (ring_struct) ring
  type (coord_struct), allocatable :: orb_(:)

  real(rdef) Q_x, Q_y
  real(rdef), allocatable :: dk1(:) 
  real(rdef) int_Q_x, int_Q_y, phi_x, phi_y

  logical ok

       allocate(orb_(0:ring%n_ele_max))       
       allocate(dk1(ring%n_ele_max))

       call choose_quads(ring, dk1)
       int_Q_x = int(ring%ele_(ring%n_ele_ring)%x%phi / twopi)
       int_Q_y = int(ring%ele_(ring%n_ele_ring)%y%phi / twopi)
       phi_x = (int_Q_x + Q_x) * twopi
       phi_y = (int_Q_y + Q_y) * twopi
       call custom_set_tune (phi_x, phi_y, dk1, ring, orb_, ok)

       deallocate(dk1)
       deallocate(orb_)


     return
     end






