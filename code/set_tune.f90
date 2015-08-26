!+
! Subroutine set_tune (phi_a_set, phi_b_set, dk1, lat, orb, ok)
!
! Subroutine to Q_tune a lat. Program will set the tunes to
! within 0.001 radian (0.06 deg).
! Note: The tune is computed with reference to the closed orbit.
!                                    
! Modules Needed:
!   use bmad
!
! Input:
!   phi_a_set -- Real(rp): Horizontal set tune (radians)
!   phi_b_set -- Real(rp): Vertical set tune (radians)
!   dk1(:)    -- Real(rp): Relative amount to vary a quad in tuning.
!                  That is, the variation will be proportional to dk1.
!                  dk1(i) relates to lat%ele(i). Those quads with a
!                  positive dk1(i) will be varied as one group and the
!                  quads with negative dk1(i) will be varied as another group.
!   orb(0)%vec(6) -- Coord_struct: If RF is off: Energy dE/E at which the tune is computed.
!
! Output:
!   lat      -- lat_struct: Q_tuned lat
!   orb(0:)  -- Coord_struct: New closed orbit.
!   ok       -- Logical: Set True if everything is ok. False otherwise.
!-

subroutine set_tune (phi_a_set, phi_b_set, dk1, lat, orb, ok)

use bmad_interface, except_dummy => set_tune
use bookkeeper_mod, only: lattice_bookkeeper, set_flags_for_changed_attribute, rf_is_on

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (ele_struct) :: ave
type (coord_struct), allocatable :: orb(:)

real(rp) phi_a_set, phi_b_set, dphi_a, dphi_b, dQ_max
real(rp) phi_a, phi_b, d_xx, d_xy, d_yx, d_yy, det
real(rp) l_beta_a, l_beta_b, dk_x, dk_y, dk1(:)

integer i, j, status

logical ok, err, rf_on

character(20) :: r_name = 'set_tune'
real(rp), dimension(2) :: phi_array

! Init

dQ_max = 0.001
ok = .false.
rf_on = rf_is_on(lat%branch(0))

do j = 1, lat%n_ele_max
  if (dk1(j) == 0) cycle
  ele => lat%ele(j)
  if (attribute_free(ele, 'K1', .false.)) cycle
  call out_io (s_warn$, r_name, 'K1 ATTRIBUTE NOT FREE TO VARY OF ELEMENT: ' // ele%name, 'WILL NOT USE THIS!')
enddo

! Q tune

do i = 1, 10

  if (.not. bmad_com%auto_bookkeeper) call lattice_bookkeeper(lat)

  if (rf_on) then
    call closed_orbit_calc (lat, orb, 6, err_flag = err)
  else
    call closed_orbit_calc (lat, orb, 4, err_flag = err)
  endif

  if (err) return

  call lat_make_mat6 (lat, -1, orb)

  call twiss_at_start(lat, status = status)
  if (status /= ok$) return

  call twiss_propagate_all (lat, err_flag = err)
  if (err) return

  phi_a = lat%ele(lat%n_ele_track)%a%phi
  phi_b = lat%ele(lat%n_ele_track)%b%phi
  dphi_a = phi_a_set - phi_a 
  dphi_b = phi_b_set - phi_b 
  if (abs(dphi_a) < dQ_max .and. abs(dphi_b) < dQ_max) then
    ok = .true.
    return
  endif

  d_xx = 0
  d_xy = 0
  d_yx = 0
  d_yy = 0

  do j = 1, lat%n_ele_max
    if (dk1(j) == 0) cycle
    call twiss_at_element (lat, j, average = ave)
    l_beta_a =  abs(dk1(j)) * ave%a%beta * ave%value(l$) / 2
    l_beta_b = -abs(dk1(j)) * ave%b%beta * ave%value(l$) / 2

    if (dk1(j) > 0) then
      d_xx = d_xx + l_beta_a
      d_yx = d_yx + l_beta_b
    else
      d_xy = d_xy + l_beta_a
      d_yy = d_yy + l_beta_b
    endif
  enddo

  det = d_xx * d_yy - d_xy * d_yx
  dk_x = (d_yy * dphi_a - d_xy * dphi_b) / det
  dk_y = (d_xx * dphi_b - d_yx * dphi_a) / det

  ! put in the changes

  do j = 1, lat%n_ele_max
    if (dk1(j) == 0) cycle
    ele => lat%ele(j)
    if (.not. attribute_free(ele, 'K1', .false.)) cycle
    if (dk1(j) > 0) then
      ele%value(k1$) = ele%value(k1$) + abs(dk1(j)) * dk_x
    else
      ele%value(k1$) = ele%value(k1$) + abs(dk1(j)) * dk_y
    endif
    call set_flags_for_changed_attribute (ele, ele%value(k1$))
  enddo

enddo

phi_array(1) = phi_a/twopi
phi_array(2) = phi_b/twopi
phi_array(1) = phi_a_set/twopi
phi_array(2) = phi_b_set/twopi
call out_io (s_error$, r_name, 'CANNOT GET TUNE RIGHT.', &
      'CURRENT TUNE: \2f\ ', &
      'SET TUNE:     \2f\ ', &
      r_array = [phi_a/twopi, phi_b/twopi, phi_a_set/twopi, phi_b_set/twopi ])


end subroutine
