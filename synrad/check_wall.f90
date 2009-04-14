
subroutine check_wall (wall, s_lat)

use synrad_struct
use synrad_interface, except => check_wall

implicit none

type (wall_struct), target :: wall
type (lat_struct) lat
type (wall_pt_struct), pointer :: wptr(:)

real(rp) s_lat

integer n, n0, n1

!

wptr => wall%pt

do n = 1, wall%n_pt_tot

  if (wptr(n-1)%s > wptr(n)%s .and. wptr(n)%type == no_alley$) then
    print *, 'ERROR: ', wall_name(wall%side), &
                    ' WALL POINTS NOT ORDERED IN ASSENDING LONGITUDINAL POSITION'
    print '(6x, i4, f11.5, 2x, a)', n-1, wptr(n-1)%s, wptr(n-1)%name
    print '(6x, i4, f11.5, 2x, a)', n, wptr(n)%s, wptr(n)%name
  endif

  n0 = wptr(n-1)%ix_pt
  n1 = wptr(n)%ix_pt
  if (wptr(n0)%s > wptr(n1)%s) then
    print *, 'ERROR: ', trim(wall_name(wall%side)), ' WALL IX_PT ARRAY NOT ORDERED'
    print '(11x, a)', 'N Ix_pt          S    Name'
    print '(6x, 2i6, f11.5, 2x, a)', n-1, n0, wptr(n0)%s, wptr(n0)%name
    print '(6x, 2i6, f11.5, 2x, a)', n,   n1, wptr(n1)%s, wptr(n1)%name
    call err_exit
  endif

  if (n0 == n1) then
    print *, 'ERROR: ', wall_name(wall%side), ' WALL IX_PT DEGENERACY:', n-1, n
    print *, '        BOTH POINT TO: ', n0, n1
    call err_exit
  endif

enddo

if (wptr(0)%s /= 0) then
  print *, 'ERROR: ', wall_name(wall%side), ' WALL DOES NOT START AT S = 0'
  call err_exit
endif

if (wptr(wall%n_pt_tot)%s /= s_lat) then
  print *, 'ERROR: ', wall_name(wall%side), ' WALL DOES NOT GO A FULL TURN'
  call err_exit
endif

end subroutine
