!+
! Subroutine SET_TUNE (phi_x_set, phi_y_set, dk1, ring, orb_, ok)
!
! Subroutine to Q_tune a ring. Program will set the tunes to
! within 0.001 radian (0.06 deg).
!                                    
! Modules Needed:
!   use bmad
!
! Input:
!   phi_x_set       -- Real(rdef): Horizontal set tune (radians)
!   phi_y_set       -- Real(rdef): Vertical set tune (radians)
!   dk1(n_ele_maxx) -- Real(rdef): Relative amount to vary a quad in tuning.
!                       dk1(i) relates to ring%ele_(i). Those quads with a
!                       positive dk1(i) will be varied as one group and the
!                       quads with negative dk1(i) will be varied as another
!                       group.
!
! Output:
!   ring      -- Ring_struct: Q_tuned ring
!   orb_(0:)  -- Coord_struct: New closed orbit.
!   ok        -- Logical: Set True if everything is ok. False otherwise.
!-

!$Id$
!$Log$
!Revision 1.5  2003/01/27 14:40:43  dcs
!bmad_version = 56
!
!Revision 1.4  2002/07/16 20:44:02  dcs
!*** empty log message ***
!
!Revision 1.3  2002/02/23 20:32:25  dcs
!Double/Single Real toggle added
!
!Revision 1.2  2001/09/27 18:31:57  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine set_tune (phi_x_set, phi_y_set, dk1, ring, orb_, ok)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ring_struct) ring
  type (ele_struct) ave
  type (coord_struct) orb_(0:)

  real(rdef) phi_x_set, phi_y_set, dphi_x, dphi_y
  real(rdef) phi_x, phi_y, d_xx, d_xy, d_yx, d_yy, det
  real(rdef) l_beta_x, l_beta_y, dk_x, dk_y, dk1(:)

  integer i, j

  logical ok

! q_tune

  do i = 1, 10

    call ring_make_mat6 (ring, -1, orb_)

    call twiss_at_start(ring)
    ok = bmad_status%ok
    if (.not. ok) return

    call twiss_propagate_all(ring)
    call closed_orbit_at_start (ring, orb_(0), 4, .true.)
    call track_all (ring, orb_)

    phi_x = ring%ele_(ring%n_ele_ring)%x%phi
    phi_y = ring%ele_(ring%n_ele_ring)%y%phi
    dphi_x = phi_x_set - phi_x 
    dphi_y = phi_y_set - phi_y 
    if (abs(dphi_x) < 0.001 .and. abs(dphi_y) < 0.001) return

    d_xx = 0
    d_xy = 0
    d_yx = 0
    d_yy = 0

    do j = 1, ring%n_ele_max
      if (dk1(j) == 0) cycle
      call twiss_at_element (ring, j, average = ave)
      l_beta_x =  abs(dk1(j)) * ave%x%beta * ave%value(l$) / 2
      l_beta_y = -abs(dk1(j)) * ave%y%beta * ave%value(l$) / 2

      if (dk1(j) > 0) then
        d_xx = d_xx + l_beta_x
        d_yx = d_yx + l_beta_y
      else
        d_xy = d_xy + l_beta_x
        d_yy = d_yy + l_beta_y
      endif
    enddo

    det = d_xx * d_yy - d_xy * d_yx
    dk_x = (d_yy * dphi_x - d_xy * dphi_y) / det
    dk_y = (d_xx * dphi_y - d_yx * dphi_x) / det

! put in the changes

    do j = 1, ring%n_ele_max
      if (dk1(j) == 0) cycle
      if (dk1(j) > 0) then
        ring%ele_(j)%value(k1$) = ring%ele_(j)%value(k1$) + abs(dk1(j)) * dk_x
      else
        ring%ele_(j)%value(k1$) = ring%ele_(j)%value(k1$) + abs(dk1(j)) * dk_y
      endif
    enddo

  enddo

  type *, 'ERROR IN SET_TUNE: CANNOT GET TUNE RIGHT.'
  type *, '      CURRENT TUNE:', phi_x/twopi, phi_y/twopi
  type *, '      SET TUNE:    ', phi_x_set/twopi, phi_y_set/twopi
  ok = .false.
  call err_exit

end subroutine
