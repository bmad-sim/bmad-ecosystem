!+
! Subroutine synrad3d_wall_to_synrad_walls (wall3d, seg_len_max, branch, walls)
!
! Routine to convert from a synrad3d wall structure to a synrad wall structure.
!
! Input:
!   wall3d       -- sr3d_wall_struct: Synrad3d wall structure.
!   seg_len_max  -- Real(rp): Maximum length of wall segments.
!   branch       -- branch_struct: lattice branch to use
!
! Output:
!   walls -- Walls_struct: synrad wall structure.
!-

subroutine synrad3d_wall_to_synrad_walls (wall3d, seg_len_max, branch, walls)

use synrad_mod, except => synrad3d_wall_to_synrad_walls
use synrad3d_utils

implicit none

type (sr3d_wall_struct) wall3d
type (walls_struct), target :: walls
type (branch_struct) branch
type (wall_struct), pointer :: minus_side, plus_side

real(rp) seg_len_max, dr_dtheta
integer i, n_pt
logical in_ante

! init

plus_side => walls%positive_x_wall
minus_side  => walls%negative_x_wall

!

n_pt = ubound(wall3d%section, 1)
allocate (minus_side%pt(0:n_pt), plus_side%pt(0:n_pt))
minus_side%n_pt_max = n_pt
plus_side%n_pt_max = n_pt

do i = 0, n_pt

  minus_side%pt(i)%s = wall3d%section(i)%s
  minus_side%pt(i)%name = wall3d%section(i)%name
  minus_side%pt(i)%phantom = .false.

  plus_side%pt(i) = minus_side%pt(i)

  call sr3d_wall_section_params (wall3d%section(i), 1.0_rp, 0.0_rp, plus_side%pt(i)%x, dr_dtheta, in_ante)
  call sr3d_wall_section_params (wall3d%section(i), -1.0_rp, 0.0_rp, minus_side%pt(i)%x, dr_dtheta, in_ante)
  minus_side%pt(i)%x = -minus_side%pt(i)%x

enddo

!

call synrad_setup_walls (walls, branch, seg_len_max)

end subroutine
