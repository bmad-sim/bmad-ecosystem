!+
! Subroutine synrad3d_wall_to_synrad_walls (wall3d, seg_len_max, s_lat, geometry, walls)
!
! Routine to convert from a synrad3d wall structure to a synrad wall structure.
!
! Input:
!   wall3d       -- sr3d_wall_struct: Synrad3d wall structure.
!   seg_len_max  -- Real(rp): Maximum length of wall segments.
!   s_lat        -- Real(rp): Lattice length
!   geometry     -- Integer: Type of lattice. open$ or closed$
!
! Output:
!   walls -- Walls_struct: synrad wall structure.
!-

subroutine synrad3d_wall_to_synrad_walls (wall3d, seg_len_max, s_lat, geometry, walls)

use synrad_mod, except => synrad3d_wall_to_synrad_walls
use synrad3d_utils

implicit none

type (sr3d_wall_struct) wall3d
type (walls_struct), target :: walls
type (wall_struct), pointer :: inside, outside

real(rp) s_lat, seg_len_max, dr_dtheta
integer i, n_pt, geometry
logical in_ante

! init

outside => walls%positive_x_wall
inside  => walls%negative_x_wall

outside%side = positive_x$
inside%side = negative_x$

!

n_pt = ubound(wall3d%section, 1)
allocate (inside%pt(0:n_pt), outside%pt(0:n_pt))
inside%n_pt_tot = n_pt
outside%n_pt_tot = n_pt

do i = 0, n_pt

  inside%pt(i)%s = wall3d%section(i)%s
  inside%pt(i)%name = wall3d%section(i)%name
  inside%pt(i)%phantom = .false.
  inside%pt(i)%type = no_alley$
  inside%pt(i)%ix_pt = i

  outside%pt(i) = inside%pt(i)

  call sr3d_wall_section_params (wall3d%section(i), 1.0_rp, 0.0_rp, outside%pt(i)%x, dr_dtheta, in_ante)
  call sr3d_wall_section_params (wall3d%section(i), -1.0_rp, 0.0_rp, inside%pt(i)%x, dr_dtheta, in_ante)
  inside%pt(i)%x = -inside%pt(i)%x

enddo

!

call delete_overlapping_wall_points (outside)
call delete_overlapping_wall_points (inside)

! check that endpoints are correct

if (abs(outside%pt(outside%n_pt_tot)%s - s_lat) > 0.01) then
  print *, 'WARNING: OUTSIDE WALL ENDS AT:', outside%pt(outside%n_pt_tot)%s
  print *, '         AND NOT AT LATTICE END OF:', s_lat
endif

if (abs(inside%pt(inside%n_pt_tot)%s - s_lat) > 0.01) then
  print *, 'WARNING: INSIDE WALL ENDS AT:', inside%pt(inside%n_pt_tot)%s
  print *, '         AND NOT AT LATTICE END OF:', s_lat
endif

outside%pt(outside%n_pt_tot)%s = s_lat
inside%pt(inside%n_pt_tot)%s = s_lat

! do some checking

call check_wall (inside, s_lat, geometry)
call check_wall (outside, s_lat, geometry)

! segment wall

call break_wall_into_segments (inside, seg_len_max)
call break_wall_into_segments (outside, seg_len_max)

end subroutine
