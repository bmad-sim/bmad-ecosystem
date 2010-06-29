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

   interface
    subroutine fft_plane(n, co, Q)
    use bmad
    implicit none
    type(coord_struct) co(:)
    integer n
    real(rp) Q(3)
    end subroutine
   end interface

  type (lat_struct) ring
  type (coord_struct), allocatable,save :: orb(:), traj(:), co(:)
  type (coord_struct)  orb0, closed_orbit, coord

  integer n_turn
  integer n_ele, i
  integer n_fft

  real(rp) x_max, y_max, amp_x_max, amp_y_max, amp_x, amp_y, amp_x0, amp_y0
  real(rp) Q_x, Q_y, Q_z
  real(rp) Q_begin(3), Q_end(3)

  call reallocate_coord(orb,ring%n_ele_max)
  call reallocate_coord(traj,n_turn)
  call reallocate_coord(co,512)

  x_max = 0.
  y_max = 0.
  amp_x_max = 0.
  amp_y_max = 0.

  orb(0) = orb0
  coord%vec = orb0%vec - closed_orbit%vec
  
  call phase_space_amplitude(ring%ele(0), coord, amp_x0, amp_y0)
  n_ele = ring%n_ele_track
  print *,' TRACK_A_PARTICLE:'
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

      traj(i)%vec = coord%vec
  end do
  if(ring%param%ixx == 0)goto 10  
  if(n_turn >=1024)then
    n_fft = 512
     forall(i=1:n_fft)co(i)%vec = traj(i)%vec
     call fft_plane(n_fft,  co, Q_begin)
     forall(i=1:n_fft)co(i)%vec = traj(n_turn-n_fft+i)%vec
     call fft_plane(n_fft,  co, Q_end)
   else
    print '(a)',' N_turn < 1024. No fft '
  endif
       
10 continue
!

  print *, 'Track: Q_x, Q_y, Q_z ', Q_x, Q_y, Q_z
  print *, 'Track: start', orb0%vec(1), orb0%vec(3)
  print *, 'Track: end  ', orb(0)%vec(1), orb(0)%vec(3),'   ',i,'turns'
  print *, 'Track: max displacement x,y ',x_max, y_max
  print *, 'Track: max amplitude ratio x,y ',amp_x_max/amp_x0, amp_y_max/amp_y0


if(ring%param%ixx == 0)then
  type('(3f10.5,2e12.4,2x,2e12.4,i10)'),Q_x,Q_y,Q_z, orb(0)%vec(1), orb(0)%vec(3), &
                 amp_x_max/amp_x0, amp_y_max/amp_y0, i
  write(21,'(3f10.5,2e12.4,2x,2e12.4,i10)')Q_x,Q_y,Q_z, orb(0)%vec(1), orb(0)%vec(3), &
                 amp_x_max/amp_x0, amp_y_max/amp_y0, i
 else
  type('(3f10.5,2e12.4,2x,2e12.4,i10,6e12.4)'),Q_x,Q_y,Q_z, orb(0)%vec(1), orb(0)%vec(3), &
                 amp_x_max/amp_x0, amp_y_max/amp_y0, i,q_begin(1:3), q_end(1:3)
  write(21,'(3f10.5,2e12.4,2x,2e12.4,i10, 6e12.4)')Q_x,Q_y,Q_z, orb(0)%vec(1), orb(0)%vec(3), &
                 amp_x_max/amp_x0, amp_y_max/amp_y0, i, q_begin(1:3), q_end(1:3)
endif


end subroutine







































