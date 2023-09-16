!+
! Function tao_read_phase_space_index (name, ixc, print_err) result (ix_ps)
!
! Returns the integer value in the array <name> at position <ixc>. This is used
! for finding a phase space index reference in the range 1 to 6.
!
! Input:
!   name      -- Character(*): character array holding the index. Must be in the range 1-6.
!   ixc       -- Integer: location within <name> to evaluate index.
!   print_err -- logical, optional: If present and False then do not print an error message
!
! Output:
!   ix_ps    -- Integer: Index at <name>(<ixc>:<ixc>). Returns 0 if bad index.
!
! Example:
!      name = r:26
!      ixc  = 3
!
! Gives:
!      ix_ps = 2
!
! Example:
!      name = mat_94
!      ixc  = 7
! Gives an error.
!-

function tao_read_phase_space_index (name, ixc, print_err) result (ix_ps)

use sim_utils

implicit none

character(*) name
integer ix_ps, ixc
character(*), parameter :: r_name = 'tao_read_phase_space_index'
logical, optional :: print_err

!

if (ixc < 1 .or. ixc > len_trim(name)) then
  if (logic_option(.true., print_err)) call out_io (s_error$, r_name, 'BAD INDEX LOCATION: \i0\ FOR: ' // name, ixc)
  ix_ps = 0
  return
endif

ix_ps = index('123456', name(ixc:ixc))
if (ix_ps == 0) then
  if (logic_option(.true., print_err)) call out_io (s_error$, r_name, 'BAD PHASE SPACE INDEX: "' // name(ixc:ixc) // '" IN STRING: ' // name)
endif

end function tao_read_phase_space_index
