!+
! Subroutine synrad3d_wall_to_synrad_walls (branch, seg_len_max, walls)
!
! Routine to convert from a synrad3d wall structure to a synrad wall structure.
!
! Input:
!   branch       -- branch_struct: lattice branch with wall3d.
!   seg_len_max  -- Real(rp): Maximum length of wall segments.
!
! Output:
!   walls -- Walls_struct: synrad wall structure.
!-

subroutine synrad3d_wall_to_synrad_walls (branch, seg_len_max, walls)

use synrad_mod, except => synrad3d_wall_to_synrad_walls
use wall3d_mod

implicit none

type (branch_struct), target :: branch
type (walls_struct), target :: walls
type (wall_struct), pointer :: minus_side, plus_side
type (wall3d_section_struct), pointer :: sec

real(rp) seg_len_max, dr_dtheta
integer i, n_pt

! init

plus_side => walls%positive_x_wall
minus_side  => walls%negative_x_wall

!

n_pt = ubound(branch%wall3d%section, 1)
allocate (minus_side%pt(0:n_pt), plus_side%pt(0:n_pt))
minus_side%n_pt_max = n_pt
plus_side%n_pt_max = n_pt

do i = 0, n_pt

  minus_side%pt(i)%s = branch%wall3d%section(i)%s
  minus_side%pt(i)%name = branch%wall3d%section(i)%name
  minus_side%pt(i)%phantom = .false.

  plus_side%pt(i) = minus_side%pt(i)
  sec => branch%wall3d%section(i)

  call calc_wall_radius (sec%v,  1.0_rp, 0.0_rp, plus_side%pt(i)%x, dr_dtheta)
  call calc_wall_radius (sec%v, -1.0_rp, 0.0_rp, minus_side%pt(i)%x, dr_dtheta)
  minus_side%pt(i)%x = -minus_side%pt(i)%x

enddo

!

call synrad_setup_walls (walls, branch, seg_len_max)

end subroutine
