!+
! function set_tune_3d (branch, target_tunes, use_phase_trombone, quad_mask, z_tune_set, group_knobs, print_err) result (everything_ok)
!
! Wrapper for set_tune and set_z_tune together.
!
! Input:
!   branch              -- branch_struct:
!   target_tunes(1:3)   -- real(rp): tunes for a, b, z modes (rad/2pi). Must include integer part.
!   quad_mask           -- character(*), optional: List of quads to not use in qtuneing.
!   use_phase_trombone  -- logical, optional: Default False. If true, use a match element in phase trombone mode to adjust the tunes.
!                            The match element must be the first element in the lattice. Use insert_phase_trombone to insert one.
!   z_tune_set          -- logical, optional: Default True. If false, do not try to set the synch tune.
!   group_knobs(2)      -- character(*), optional: If set non-blank, use these group elements for tuning.
!   print_err           -- logical, optional: Print error message if there is a problem? Default is True.
!
! Output:
!   branch              -- branch_struct: with adjusted quads and RF to match desired tunes.
!   everything_ok       -- logical: Returns true or false if set was successful.  
!-

function set_tune_3d (branch, target_tunes, quad_mask, use_phase_trombone, z_tune_set, group_knobs, print_err) result (everything_ok)

use bmad

implicit none

type(branch_struct), target :: branch
type (ele_struct), pointer :: ele
type(coord_struct), allocatable :: co(:)
type (ele_pointer_struct), allocatable :: eles(:)

real(rp) target_tunes(3), dQ_a, dQ_b, dQ_max
real(rp), allocatable :: dk1(:)
integer n, status, i
logical, optional :: use_phase_trombone, z_tune_set, print_err
logical everything_ok, err, use_groups

character(*), optional :: quad_mask, group_knobs(2)
character(*), parameter :: r_name = 'set_tune_3d'

!

everything_ok = .true.
dQ_max = 0.0001

if (all(target_tunes < 1)) then
  call out_io (s_fatal$, r_name, 'Only fractional tunes given for target_tunes!', &
                                 'Must supply integer + fractional tunes.', &
                                 'Stopping here...')
  stop
endif

! If user has not specified one or more tunes, set target

if (target_tunes(1) < 1.e-12) target_tunes(1) = branch%ele(branch%n_ele_track)%a%phi / twopi
if (target_tunes(2) < 1.e-12) target_tunes(2) = branch%ele(branch%n_ele_track)%b%phi / twopi
if (abs(target_tunes(3)) < 1.e-12) target_tunes(3) = branch%z%tune / twopi

! Phase trombone

if (logic_option(.false., use_phase_trombone)) then
  ele => branch%ele(1)
  n = branch%n_ele_track
  do i = 1, 20
    call twiss_and_track(branch%lat, co, status, branch%ix_branch, print_err = print_err)
    if (status == ok$) then
      dQ_a = target_tunes(1) - branch%ele(n)%a%phi/twopi
      dQ_b = target_tunes(2) - branch%ele(n)%b%phi/twopi
    else
      dQ_a = 0.01   ! Try to get off this resonance
      dQ_b = 0.01
    endif
    if (abs(dQ_a) <= dQ_max .and. abs(dQ_b) <= dQ_max) return
    ele%value(dphi_a$) = ele%value(dphi_a$) + twopi * dQ_a
    ele%value(dphi_b$) = ele%value(dphi_b$) + twopi * dQ_b
    call make_mat6(ele, branch%param, co(0))
  enddo
  everything_ok = .false.
  return
endif

!

use_groups = .false.
if (present(group_knobs)) use_groups = (group_knobs(1) /= '')

if (use_groups) then
  everything_ok = set_tune_via_group_knobs(twopi*target_tunes, branch, group_knobs, co, print_err)

else
  allocate(dk1(branch%n_ele_max))
  call choose_quads_for_set_tune(branch, dk1, eles, quad_mask, err)
  if (err) then
    call out_io (s_error$, r_name, &
      'CANNOT FIND A QUAD WITH BETA_A < BETA_B AND A QUAD WITH BETA_A > BETA_B (BOTH WITH NO TILT).')
    return
  endif

  everything_ok = set_tune(twopi*target_tunes(1), twopi*target_tunes(2), dk1, eles, branch, co, print_err)
endif

!

if (logic_option(.true., z_tune_set)) call set_z_tune(branch, twopi*target_tunes(3), print_err = print_err)

end function set_tune_3d
