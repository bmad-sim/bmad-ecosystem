!........................................................................
!+
! Subroutine : track_a_particle (closed_orbit, orb0, n_turn, ring, Q_x, Q_y, Q_z)
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

subroutine track_a_particle (closed_orbit, orb0, n_turn, ring, Q_x, Q_y, Q_z)

  use bmad_struct
  use bmad_interface
  use bmad_utils_mod

  implicit none

  type (lat_struct) ring
  type (coord_struct), allocatable,save :: orb(:)
  type (coord_struct)  orb0, closed_orbit, coord

  integer n_turn
  integer n_ele, i

  real(rp) x_max, y_max, amp_x_max, amp_y_max, amp_x, amp_y, amp_x0, amp_y0
  real(rp) Q_x, Q_y, Q_z

  call reallocate_coord(orb,ring%n_ele_max)

  x_max = 0.
  y_max = 0.
  amp_x_max = 0.
  amp_y_max = 0.

  orb(0) = orb0
  coord%vec = orb0%vec - closed_orbit%vec
  
  call phase_space_amplitude(ring%ele(0), coord, amp_x0, amp_y0)
  n_ele = ring%n_ele_track
  type *,' TRACK_A_PARTICLE:'
  type ('(a15,6e10.2)'),'  closed_orbit ', closed_orbit%vec
  type ('(a15,6e10.2)'),'  orb(0)      ', orb(0)%vec
  do i = 1, n_turn
    call track_all (ring, orb)
    if (ring%param%lost)goto 10
      orb(0) = orb(n_ele)
      coord%vec = orb(0)%vec - closed_orbit%vec
      call phase_space_amplitude (ring%ele(0), coord, amp_x, amp_y)
      x_max = max(x_max, abs(orb(n_ele)%vec(1)))
      y_max = max(y_max, abs(orb(n_ele)%vec(3)))
      amp_x_max = max(amp_x_max, amp_x)
      amp_y_max = max(amp_y_max, amp_y)
  end do

10 continue
!

  type *, 'Track: Q_x, Q_y, Q_z ', Q_x, Q_y, Q_z
  type *, 'Track: start', orb0%vec(1), orb0%vec(3)
  type *, 'Track: end  ', orb(0)%vec(1), orb(0)%vec(3),'   ',i,'turns'
  type *, 'Track: max displacement x,y ',x_max, y_max
  type *, 'Track: max amplitude ratio x,y ',amp_x_max/amp_x0, amp_y_max/amp_y0

  type('(3f10.5,2e12.4,2x,2e12.4,i10)'),Q_x,Q_y,Q_z, orb(0)%vec(1), orb(0)%vec(3), &
                 amp_x_max/amp_x0, amp_y_max/amp_y0, i
  write(21,'(3f10.5,2e12.4,2x,2e12.4,i10)')Q_x,Q_y,Q_z, orb(0)%vec(1), orb(0)%vec(3), &
                 amp_x_max/amp_x0, amp_y_max/amp_y0, i


end subroutine







































