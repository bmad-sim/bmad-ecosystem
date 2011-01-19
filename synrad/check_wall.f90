
subroutine check_wall (wall, s_lat, lat_type)

use synrad_struct
use synrad_interface, except => check_wall

implicit none

type (wall_struct), target :: wall
type (lat_struct) lat
type (wall_pt_struct), pointer :: pt(:)

real(rp) s_lat
integer lat_type, np, n, n0, n1

character(16) :: r_name = 'check_wall'
character(60) line1, line2

!

pt => wall%pt
np = wall%n_pt_tot

do n = 1, np

  if (pt(n-1)%s > pt(n)%s .and. pt(n)%type == no_alley$) then
    write (line1, '(6x, i4, f11.5, 2x, a)') n-1, pt(n-1)%s, pt(n-1)%name
    write (line2, '(6x, i4, f11.5, 2x, a)') n, pt(n)%s, pt(n)%name
    call out_io (s_fatal$, r_name, wall_name(wall%side), &
                    ' WALL POINTS NOT ORDERED IN ASSENDING LONGITUDINAL POSITION', line1, line2)
  endif

  n0 = pt(n-1)%ix_pt
  n1 = pt(n)%ix_pt
  if (pt(n0)%s > pt(n1)%s) then
    write (line1, '(6x, 2i6, f11.5, 2x, a)') n-1, n0, pt(n0)%s, pt(n0)%name
    write (line2, '(6x, 2i6, f11.5, 2x, a)') n,   n1, pt(n1)%s, pt(n1)%name
    call out_io (s_fatal$, r_name, &
                    trim(wall_name(wall%side)) // ' WALL IX_PT ARRAY NOT ORDERED', &
                    '           N Ix_pt          S    Name', line1, line2)
    call err_exit
  endif

  if (n0 == n1) then
    call out_io (s_fatal$, r_name, &
                    wall_name(wall%side) // ' WALL IX_PT DEGENERACY: \i0\  \i0\ ', &
                    'BOTH POINT TO: \i0\  \i0\ ', i_array = [n-1, n, n0, n1])
    call err_exit
  endif

enddo

if (pt(0)%s /= 0) then
  call out_io (s_fatal$, r_name, wall_name(wall%side) // ' WALL DOES NOT START AT S = 0')
  call err_exit
endif

if (pt(np)%s /= s_lat) then
  call out_io (s_fatal$, r_name, wall_name(wall%side) // ' WALL DOES NOT GO A FULL TURN')
  call err_exit
endif

! If the lattice is circular then the end x must be equal to the start x.

if (lat_type == circular_lattice$) then
  if (pt(np)%x /= pt(0)%x) then
    call out_io (s_fatal$, r_name, 'BEGINNING X DOES NOT EQUAL FINAL X FOR: ' // wall_name(wall%side))
    call err_exit
  endif
endif

end subroutine
