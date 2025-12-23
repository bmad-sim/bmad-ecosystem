!+
! Function set_tune (phi_a_set, phi_b_set, dk1, eles, branch, orb, print_err) result (ok)
!
! Function to Q_tune a lattice branch. The tunes will be set to within 0.001 radian (0.06 deg).
! Note: The tune is computed with reference to the closed orbit.
!                                    
! Input:
!   phi_a_set     -- real(rp): Horizontal set tune (radians)
!   phi_b_set     -- real(rp): Vertical set tune (radians)
!   dk1(:)        -- real(rp): Relative amount to vary a quad in tuning. The variation will be proportional to dk1. 
!                      Those quads with a positive dk1(i) will be varied as one group and the quads with negative 
!                      dk1(i) will be varied as another group. 
!                      The routine choose_quads_for_set_tune can be used to calculate values for dk1.
!   eles(:)       -- ele_pointer_struct: eles(i)%ele points to quadrupole corresponding to dk1(i).
!   branch        -- branch_struct: Lattice branch to tune.
!   orb(0)%vec(6) -- Coord_struct: If RF is off: Energy dE/E at which the tune is computed.
!   print_err     -- logical, optional: Print error message if there is a problem? Default is True.
!
! Output:
!   branch        -- branch_struct: Q_tuned lattice branch
!   orb(0:)       -- coord_struct: New closed orbit.
!   ok            -- logical: Set True if everything is ok. False otherwise.
!-

function set_tune (phi_a_set, phi_b_set, dk1, eles, branch, orb, print_err) result (ok)

use bmad_interface, except_dummy => set_tune

implicit none

type (ele_pointer_struct) :: eles(:)
type (branch_struct), target :: branch
type (ele_struct), pointer :: ele
type (ele_struct) :: ave
type (coord_struct), allocatable :: orb(:)

real(rp) phi_a_set, phi_b_set, dphi_a, dphi_b, dQ_max
real(rp) phi_a, phi_b, d_xx, d_xy, d_yx, d_yy, det
real(rp) l_beta_a, l_beta_b, dk_x, dk_y, dk1(:)

integer i, j, status

logical, optional :: print_err
logical ok, err, rf_on, master_saved

character(*), parameter :: r_name = 'set_tune'
real(rp), dimension(2) :: phi_array

! Init

dQ_max = 0.001
ok = .false.
rf_on = rf_is_on(branch)

do j = 1, size(dk1)
  if (dk1(j) == 0) cycle
  ele => pointer_to_ele(branch%lat, eles(j)%loc)

  master_saved = ele%field_master
  ele%field_master = .false.
  if (attribute_free(ele, 'K1', .false.)) cycle
  call out_io (s_warn$, r_name, 'K1 ATTRIBUTE NOT FREE TO VARY OF ELEMENT: ' // ele%name, 'WILL NOT USE THIS!')
  ele%field_master = master_saved
enddo

! Q tune

do i = 1, 10
  if (rf_on) then
    call closed_orbit_calc (branch%lat, orb, 6, 1, branch%ix_branch, err, print_err)
  else
    call closed_orbit_calc (branch%lat, orb, 4, 1, branch%ix_branch, err, print_err)
  endif

  if (err) return

  call lat_make_mat6 (branch%lat, -1, orb, branch%ix_branch)

  call twiss_at_start(branch%lat, status, branch%ix_branch, print_err)
  if (status /= ok$) return

  call twiss_propagate_all (branch%lat, branch%ix_branch, err)
  if (err) return

  phi_a = branch%ele(branch%n_ele_track)%a%phi
  phi_b = branch%ele(branch%n_ele_track)%b%phi
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

  do j = 1, size(dk1)
    if (dk1(j) == 0) cycle
    ele => pointer_to_ele(branch%lat, eles(j)%loc)

    call twiss_at_element (ele, average = ave)
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

  do j = 1, size(dk1)
    if (dk1(j) == 0) cycle
    ele => pointer_to_ele(branch%lat, eles(j)%loc)

    master_saved = ele%field_master
    ele%field_master = .false.
    if (.not. attribute_free(ele, 'K1', .false.)) cycle
    if (dk1(j) > 0) then
      ele%value(k1$) = ele%value(k1$) + abs(dk1(j)) * dk_x
    else
      ele%value(k1$) = ele%value(k1$) + abs(dk1(j)) * dk_y
    endif
    call set_flags_for_changed_attribute (ele, ele%value(k1$))
    call attribute_bookkeeper (ele, .true.)
    ele%field_master = master_saved
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

end function
