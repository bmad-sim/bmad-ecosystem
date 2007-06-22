subroutine delete_overlapping_wall_points (wall)

  use sr_struct
  use sr_interface

  implicit none

  type (wall_struct), target :: wall
  integer i

!

  i = 2
  do while (i .le. wall%n_pt_tot)
    if (wall%pt(i)%x == wall%pt(i-1)%x .and. &
                             wall%pt(i)%s == wall%pt(i-1)%s) then
      type *, 'Note: Deleting overlapping point on: ', wall_name(wall%side)
      type *,  i, wall%pt(i)%s, wall%pt(i)%x, '  ', wall%pt(i)%name
      wall%pt(i:wall%n_pt_tot-1) = wall%pt(i+1:wall%n_pt_tot)
      wall%n_pt_tot = wall%n_pt_tot - 1
    else
      i = i + 1
    endif
  enddo

end subroutine
