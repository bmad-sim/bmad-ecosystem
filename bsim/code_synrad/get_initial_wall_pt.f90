!+
! subroutine get_initial_wall_pt (ray, wall, ix_pt)
!
! Routine to find a wall point on the wall that is near the ray but in back of it.
! If ray%direction = 1 then return a wall point that has
!   wall%pt(ix_pt-1)%s <= ray%now%s
! And for ray%direction = -1:
!   wall%pt(ix_pt)%s >= ray%now%s
!
! Wall points in an alley section are ignored. 
! That is, wall%pt(ix_pt)%next_to_alley is guaranteed to be False.
!
! Modules needed:
!   use synrad_mod
!
! Input:
!   ray    -- ray_struct:
!   wall   -- wall_struct: An outside or inside wall.
!
! Output:
!   ix_pt  -- Integer: wall%pt(:) index.
!-

subroutine get_initial_wall_pt (ray, wall, ix_pt)

use synrad_struct
use synrad_interface, except => get_initial_wall_pt

implicit none

type (ray_struct) ray
type (wall_struct) wall

integer ix_pt, ix0, ix1, ix2

!

if (wall%n_pt_max == 0) then
  print *, 'There are no points in the wall!'
  print *, 'You should check the wall first with check_wall!'
  if (global_com%exit_on_error) call err_exit
endif

! edge cases

if (ray%now%s == wall%pt(wall%n_pt_max)%s) then
  ix_pt = wall%n_pt_max
  return
endif

if (ray%now%s == 0) then
  ix_pt = 1
  return
endif

! normal case. divide and conquer.

ix0 = 0
ix2 = wall%n_pt_max

do
  ix1 = (ix2 + ix0) / 2

  if (wall%pt(ix1)%s < ray%now%s) then
    ix0 = ix1
  elseif (wall%pt(ix1)%s > ray%now%s) then
    ix2 = ix1
  elseif (ray%direction == 1) then   ! here wall%pt(ix1)%s == ray%now%s
    ix_pt = ix1
    exit
  else
    ix_pt = ix1 + 1
    exit
  endif

  if (ix2 - ix0 == 1) then
    ix_pt = ix2        ! pt ix2 contains ray%now%s.
    exit
  endif
enddo

! If in an alley then find a point just outside the alley

do
  if (.not. wall%pt(ix_pt)%next_to_alley) exit
  ix_pt = ix_pt - ray%direction
enddo

end subroutine
