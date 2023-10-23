!+
! Function set_tune_via_group_knobs (phi_set, branch, group_knobs, orb, print_err) result (ok)
!
! Function to Q_tune a lattice branch. The tunes will be set to within 0.001 radian (0.06 deg).
! Note: The tune is computed with reference to the closed orbit.
!                                    
! Input:
!   phi_set(2)      -- real(rp): Set tunes (radians).
!   branch          -- branch_struct: Lattice branch to tune.
!   group_knobs(2)  -- character(*): Names of group knobs to vary.
!   orb(0)%vec(6)   -- Coord_struct: If RF is off: Energy dE/E at which the tune is computed.
!   print_err       -- logical, optional: Print error message if there is a problem? Default is True.
!
! Output:
!   branch          -- branch_struct: Q_tuned lattice branch
!   orb(0:)         -- coord_struct: New closed orbit.
!   ok              -- logical: Set True if everything is ok. False otherwise.
!-

function set_tune_via_group_knobs (phi_set, branch, group_knobs, orb, print_err) result (ok)

use bmad_interface, except_dummy => set_tune_via_group_knobs

implicit none

type (branch_struct), target :: branch
type (ele_struct), pointer :: group1, group2
type (ele_struct) :: ave
type (ele_struct), pointer :: ele
type (coord_struct), allocatable :: orb(:)

real(rp) phi_set(2), dphi_a, dphi_b, dQ_max
real(rp) phi_a, phi_b, d_a1, d_a2, d_b1, d_b2, det
real(rp) d1, d2, del0
real(rp) :: phi_array(2)
real(rp), allocatable :: deriv1(:), deriv2(:), kinit(:)

integer i, j, status, n_loc

logical, optional :: print_err
logical ok, err, rf_on, master_saved

character(*) group_knobs(2)
character(*), parameter :: r_name = 'set_tune_via_group_knobs'

! Init

dQ_max = 0.001
del0 = 0.001
ok = .false.
rf_on = rf_is_on(branch)

allocate(kinit(branch%n_ele_track), deriv1(branch%n_ele_track), deriv2(branch%n_ele_track))
deriv1 = 0;  deriv2 = 0

group1 => find_group(group_knobs(1), branch, del0, kinit, deriv1, err); if (err) return
group2 => find_group(group_knobs(2), branch, del0, kinit, deriv2, err); if (err) return

! Q tune

do i = 1, 10
  call lattice_bookkeeper(branch%lat)

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
  dphi_a = phi_set(1) - phi_a 
  dphi_b = phi_set(2) - phi_b 
  if (abs(dphi_a) < dQ_max .and. abs(dphi_b) < dQ_max) then
    ok = .true.
    return
  endif

  d_a1 = 0
  d_a2 = 0
  d_b1 = 0
  d_b2 = 0

  do j = 1, branch%n_ele_track
    if (deriv1(j) == 0 .and. deriv2(j) == 0) cycle
    ele => branch%ele(j)
    call twiss_at_element (ele, average = ave)

    if (deriv1(j) /= 0) then
      d_a1 = d_a1 + deriv1(j) * ave%a%beta * ave%value(l$) / 2
      d_b1 = d_b1 - deriv1(j) * ave%b%beta * ave%value(l$) / 2
    endif

    if (deriv2(j) /= 0) then
      d_a2 = d_a2 + deriv2(j) * ave%a%beta * ave%value(l$) / 2
      d_b2 = d_b2 - deriv2(j) * ave%b%beta * ave%value(l$) / 2
    endif
  enddo

  det = d_a1 * d_b2 - d_a2 * d_b1
  d1 = (d_b2 * dphi_a - d_b1 * dphi_b) / det
  d2 = (d_a1 * dphi_b - d_a2 * dphi_a) / det

  ! Put in the changes

  group1%control%var(1)%value = group1%control%var(1)%value + d1
  group2%control%var(1)%value = group2%control%var(1)%value + d2

  call set_flags_for_changed_attribute (group1)
  call set_flags_for_changed_attribute (group2)
enddo

phi_array(1) = phi_a/twopi
phi_array(2) = phi_b/twopi
phi_array(1) = phi_set(1)/twopi
phi_array(2) = phi_set(2)/twopi
call out_io (s_error$, r_name, 'CANNOT GET TUNE RIGHT.', &
      'CURRENT TUNE: \2f\ ', &
      'SET TUNE:     \2f\ ', &
      r_array = [phi_a/twopi, phi_b/twopi, phi_set(1)/twopi, phi_set(2)/twopi ])

!--------------------------------------------------------------------------------
contains

function find_group(name, branch, del0, kinit, deriv, err) result (group_ele)

type (branch_struct) branch
type (ele_struct), pointer :: group_ele, ele
type (ele_pointer_struct), allocatable ::eles(:)
character(*) name
real(rp) del0, kinit(:), deriv(:)
integer i, j, n_loc
logical err

!

call lat_ele_locator(name, branch%lat, eles, n_loc, err, ix_dflt_branch = branch%ix_branch)
if (err) return
err = .true.

if (n_loc == 0) then
  call out_io(s_error$, r_name, 'Group element not found: ' // ele%name)
  return
endif

if (n_loc > 1) then
  call out_io(s_error$, r_name, 'Multiple lattice elements match group name: ' // ele%name)
  return
endif

group_ele => eles(1)%ele
if (group_ele%key /= group$) then
  call out_io(s_error$, r_name, 'Element is not a group type element which is needed for tune variation: ' // ele%name)
  return
endif

!

kinit = branch%ele(1:branch%n_ele_track)%value(k1$)

group_ele%control%var(1)%value = group_ele%control%var(1)%value + del0
call control_bookkeeper(branch%lat, group_ele)

do i = 1, branch%n_ele_track
  ele => branch%ele(i)
  if (ele%key /= quadrupole$) cycle
  if (ele%value(tilt$) /= 0) cycle
  if (ele%value(k1$) == kinit(i)) cycle
  deriv(i) = (ele%value(k1$) - kinit(i)) / del0
enddo

group_ele%control%var(1)%value = group_ele%control%var(1)%value - del0
call control_bookkeeper(branch%lat, group_ele)

err = .false.

end function find_group

end function
