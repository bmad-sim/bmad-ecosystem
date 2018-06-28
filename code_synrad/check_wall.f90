
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
character(100) line1, line2

!

pt => wall%pt
np = wall%n_pt_max

do n = 1, np

  if (pt(n)%s < pt(n-1)%s - 1) then
    write (line1, '(6x, i4, f11.5, 2x, a)') n-1, pt(n-1)%s, pt(n-1)%name
    write (line2, '(6x, i4, f11.5, 2x, a)') n,   pt(n)%s,   pt(n)%name
    call out_io (s_warn$, r_name, wall_name(wall%side), &
                    ' LONGITUDINAL STEP BACK OF WALL POINT IS OVER 1 METER:', line1, line2)
  endif

enddo

!

if (pt(0)%s /= 0) then
  call out_io (s_fatal$, r_name, wall_name(wall%side) // ' WALL DOES NOT START AT S = 0')
  if (global_com%exit_on_error) call err_exit
endif

if (pt(np)%s /= s_lat) then
  call out_io (s_fatal$, r_name, wall_name(wall%side) // ' WALL DOES NOT GO A FULL TURN')
  if (global_com%exit_on_error) call err_exit
endif

! If the lattice is circular then the end x must be equal to the start x.

if (lat_type == closed$) then
  if (abs(pt(np)%x - pt(0)%x) > 1d-3 * bmad_com%significant_length) then
    call out_io (s_fatal$, r_name, 'BEGINNING X DOES NOT EQUAL FINAL X FOR: ' // wall_name(wall%side))
    if (global_com%exit_on_error) call err_exit
  endif
endif

end subroutine
