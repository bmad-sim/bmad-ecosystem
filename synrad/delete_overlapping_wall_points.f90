subroutine delete_overlapping_wall_points (wall)

use synrad_struct
use synrad_interface, except => delete_overlapping_wall_points

implicit none

type (wall_struct), target :: wall
integer i

!

i = 1
do while (i .le. wall%n_pt_tot)
  if (wall%pt(i)%x == wall%pt(i-1)%x .and. wall%pt(i)%s == wall%pt(i-1)%s) then
    print *, 'Note: Deleting overlapping point on: ', wall_name(wall%side)
    print *,  i, wall%pt(i)%s, wall%pt(i)%x, '  ', wall%pt(i)%name
    wall%pt(i:wall%n_pt_tot-1) = wall%pt(i+1:wall%n_pt_tot)
    wall%n_pt_tot = wall%n_pt_tot - 1
  else
    i = i + 1
  endif
enddo

! Fix %ix_pt

do i = lbound(wall%pt, 1), ubound(wall%pt, 1)
  wall%pt(i)%ix_pt = i
enddo

end subroutine
