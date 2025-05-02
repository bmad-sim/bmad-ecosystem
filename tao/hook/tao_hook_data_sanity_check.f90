!+
! Function tao_hook_data_sanity_check (found, datum, print_err, default_data_type, uni) result (is_valid)
!
! Hook routine to check if a custom datum is internally consistent.
! This routine is called by tao_data_sanity_check. See this routine for more details.
!
! Input:
!   datum               -- tao_hook_data_struct: Datum to check.
!   print_err           -- logical: Print error message if data is not valid?
!   default_data_type   -- character(*): Default data type associated with the datum's d2 structure.
!   uni                 -- tao_universe_struct, optional: Universe to use instead of datum%d1%d2%ix_universe
!
! Output:
!   found     -- logical: If set True, tao_data_sanity_check will not proform any checks.
!                 If set False, Not a custom datum and tao_data_sanity_check will proform the standard checks. 
!   is_valid  -- logical: True if internally consistent. Value ignored if found is set False.
!-

function tao_hook_data_sanity_check (found, datum, print_err, default_data_type, uni) result (is_valid)

use tao_interface

implicit none

type (tao_data_struct) datum
type (tao_universe_struct), optional, target :: uni
type (branch_struct), pointer :: branch
type (tao_universe_struct), pointer :: u

logical found, print_err, is_valid
character(*) default_data_type
character(*), parameter :: r_name = 'tao_hook_data_sanity_check'

!

found = .false.
return

!

if (present(uni)) then
  u => uni
else
  u => s%u(datum%d1%d2%ix_universe)
endif

branch => u%design%lat%branch(datum%ix_branch)

end function tao_hook_data_sanity_check

