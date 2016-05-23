program co_test

use bmad

implicit none

bmad_com%radiation_damping_on = .true.

open (1, file = 'output.now')

call run1 ('ad05f_pos_RF.lat', 1)
call run1 ('ad05f_pos_RF.lat', 2)
call run1 ('ad05f_neg_RF.lat', 1)
call run1 ('ad05f_neg_RF.lat', 2)

!-------------------------------------------------------------------
contains

subroutine run1(lat_file, which)

type (lat_struct) lat
type (coord_struct), allocatable :: orb(:)

character(*) lat_file
integer which

!

call bmad_parser(lat_file, lat)

if (which == 2) then
  call closed_orbit_calc(lat, orb, 6)
  call lat_make_mat6(lat, -1, orb)
  call twiss_at_start(lat)
  write (1, '(3a, i0, a, t35, a, t50, f12.8)') '"', trim(lat_file), '-', which,  '-A Z"', 'ABS 1E-10', orb(0)%vec(5)
  deallocate(orb)
endif

call closed_orbit_calc(lat, orb, 6)
write (1, '(3a, i0, a, t35, a, t50, f12.8)') '"', trim(lat_file), '-', which,  '-B Z"', 'ABS 1E-10', orb(0)%vec(5)
call lat_make_mat6(lat, -1, orb)
call twiss_at_start(lat)
call twiss_propagate_all(lat)
call calc_z_tune(lat)
write (1, '(3a, i0, a, t35, a, t50, f14.10)') '"', trim(lat_file), '-', which,  '-B Z_tune"', 'ABS 1E-10', lat%z%tune / twopi

end subroutine

end program
