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
!   phi_x_set -- Real(rp): Horizontal set tune (radians)
!   phi_y_set -- Real(rp): Vertical set tune (radians)
!   dk1(:)    -- Real(rp): Relative amount to vary a quad in tuning.
!                  dk1(i) relates to ring%ele_(i). Those quads with a
!                  positive dk1(i) will be varied as one group and the
!                  quads with negative dk1(i) will be varied as another group.
!   orb_(0)%vec(6) -- Coord_struct: Energy dE/E at which the tune is computed.
!
! Output:
!   ring      -- Ring_struct: Q_tuned ring
!   orb_(0:)  -- Coord_struct: New closed orbit.
!   ok        -- Logical: Set True if everything is ok. False otherwise.
!-

#include "CESR_platform.inc"

subroutine set_tune (phi_x_set, phi_y_set, dk1, ring, orb_, ok)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ring_struct) ring
  type (ele_struct), save :: ave
  type (coord_struct), allocatable :: orb_(:)

  real(rp) phi_x_set, phi_y_set, dphi_x, dphi_y
  real(rp) phi_x, phi_y, d_xx, d_xy, d_yx, d_yy, det
  real(rp) l_beta_x, l_beta_y, dk_x, dk_y, dk1(:)

  integer i, j

  logical ok

  character(20) :: r_name = 'set_tune'
  real(rp), dimension(2) :: phi_array

! q_tune

  call closed_orbit_at_start (ring, orb_(0), 4, .true.)
  ok = bmad_status%ok
  if (.not. ok) return
  call track_all (ring, orb_)

  do i = 1, 10

    call ring_make_mat6 (ring, -1, orb_)

    call twiss_at_start(ring)
    ok = bmad_status%ok
    if (.not. ok) return

    call twiss_propagate_all(ring)
    call closed_orbit_at_start (ring, orb_(0), 4, .true.)
    ok = bmad_status%ok
    if (.not. ok) return
    call track_all (ring, orb_)

    phi_x = ring%ele_(ring%n_ele_use)%x%phi
    phi_y = ring%ele_(ring%n_ele_use)%y%phi
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
      if (associated(ring%ele_(j)%taylor(1)%term)) &
                                        call kill_taylor (ring%ele_(j)%taylor)
    enddo

  enddo

  phi_array(1) = phi_x/twopi
  phi_array(2) = phi_y/twopi
  phi_array(1) = phi_x_set/twopi
  phi_array(2) = phi_y_set/twopi
  call out_io (s_error$, r_name, (/ 'CANNOT GET TUNE RIGHT.', &
        'CURRENT TUNE: \2f\ ', &
        'SET TUNE:     \2f\ ' /), &
        (/ phi_x/twopi, phi_y/twopi, phi_x_set/twopi, phi_y_set/twopi /) )
  ok = .false.

end subroutine
