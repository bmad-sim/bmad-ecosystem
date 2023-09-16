!+
! Function coords_curvilinear_to_floor (xys, branch, err_flag) result (global)
!
! Routine to find the global position of a local lab (x, y, s) position.
! s = position from beginning of lattice branch.
!
! Input:
!   xys(3)      -- real(rp): (x, y, s) lab frame position vector.
!   branch      -- branch_struct: Lattice branch that defines the local reference coordinates.
!
! Output:
!   global      -- floor_position_struct: Global floor position corresponding to (x, y, s)
!               --    %w    -- W matrix to transform vectors: v_global = w_mat * v_local
!   err_flag    -- logical: Set True if global floor position cannot be computed.
!-

function coords_curvilinear_to_floor (xys, branch, err_flag) result (global)

use bmad_interface, dummy => coords_curvilinear_to_floor

implicit none

type (branch_struct), target :: branch
type (ele_struct), pointer :: ele
type (floor_position_struct) global, local

real(rp) xys(3), w_mat(3,3)
integer ix_ele
logical err_flag

!

ix_ele = element_at_s (branch%lat, xys(3), .true., branch%ix_branch, err_flag)
if (err_flag) return
ele => branch%ele(ix_ele)

if (ele%orientation == 1) then
  local%r = [xys(1), xys(2), xys(3) - ele%s_start]
else
  local%r = [xys(1), xys(2), ele%s - xys(3)]
endif

global = coords_local_curvilinear_to_floor (local, ele, w_mat = w_mat)
global%w = w_mat

end function coords_curvilinear_to_floor
