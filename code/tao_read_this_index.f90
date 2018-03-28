!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Function tao_read_this_index (name, ixc) result (ix)
!
! Returns the integer value in the array <name> at position <ixc>. This is used
! for finding a 6-dimensional index reference so any value at <ixc> greater than
! 6 returns an error.
!
! Input:
!  name     -- Character(*): character array holding the index
!  ixc      -- Integer: location within <name> to evaluate
!
! Output:
!  ix       -- Integer: Index at <name>(<ixc>:<ixc>)
!
! Example:
!      name = r:26
!      ixc  = 3
!
! Gives:
!      ix = 2
!
! Example:
!      name = mat_94
!      ixc  = 7
! Gives:
!      Error: "BAD INDEX CONSTRAINT: mat_94"
!-

function tao_read_this_index (name, ixc) result (ix)

use bmad_struct

implicit none

character(*) name
integer ix, ixc
character(20) :: r_name = 'tao_read_this_index'

ix = index('123456', name(ixc:ixc))
if (ix == 0) then
  call out_io (s_abort$, r_name, 'BAD INDEX CONSTRAINT: ' // name)
  call err_exit
endif

end function tao_read_this_index
