!+
! Subroutine check_if_s_in_bounds (branch, s, err_flag, translated_s, print_err)
!
! Routine to check if a given longitudinal position s is within the bounds of a given branch of a lattice.
! For linear branches the bounds are normally [0, branch_length].
! For circular branches negative s values do make sense so the bounds 
!   are normally [-branch_length, branch_length].
!
! "Normally" means that starting s-position in the branch is zero. This routine does
! adjust for non-zero starting s-positions.
!
! This routine will bomb the program if global_com%exit_on_error is True.
!
! Optionally: translated_s is a translated longitudinal position which is normally
! in the range [0, branch_length].
!
! Input:
!   branch        -- branch_struct: Branch
!   s             -- Real(rp): longitudinal position in the given branch.
!   print_err     -- logical, optional: Print error message if there is an error? Default is True.
!   
! Output:
!   err_flag      -- Logical: Set True if s position is out-of-bounds. False otherwise.
!   translated_s  -- Real(rp), optional: position translated to the range [0, branch_length]
!-

subroutine check_if_s_in_bounds (branch, s, err_flag, translated_s, print_err)

use bmad_struct

implicit none

type (branch_struct) branch

real(rp) s, ss, s_min, s_max, ds_fudge, s_bound
real(rp), optional :: translated_s

logical err_flag
logical, optional :: print_err

character(*), parameter :: r_name = 'check_if_s_in_bounds'

! Setup

s_min = branch%ele(0)%s
s_max = branch%ele(branch%n_ele_track)%s 
ds_fudge = bmad_com%significant_length
err_flag = .false.
ss = s

! Check

if (s > s_max + ds_fudge) then
  err_flag = .true.
  s_bound = s_max
elseif (branch%param%geometry == closed$) then
  if (s < s_min - (s_max - s_min) - ds_fudge) then
    err_flag = .true.
    s_bound = s_min - (s_max - s_min)
  endif
  if (s < s_min) ss = s + (s_max - s_min)
elseif (s < s_min - ds_fudge) then
  err_flag = .true.
  s_bound = s_min
endif

! Finish

if (err_flag .and. logic_option(.true., print_err)) then
  call out_io (s_fatal$, r_name, &
        'S-POSITION \es20.12\ PAST EDGE OF LATTICE. ' , &
        'PAST LATTICE EDGE AT: \es20.12\ ', r_array = [s, s_bound])
endif

if (present(translated_s)) translated_s = ss

end subroutine check_if_s_in_bounds 
