program rf_test

use bmad

implicit none

type (lat_struct), target :: lat
type (coord_struct), allocatable :: orbit(:)
type (ele_struct), pointer :: ele

integer n
character(100) file_name

!

open (1, file = 'output.now', recl = 200)


file_name = 'rf_test.bmad'
call bmad_parser(file_name, lat)

allocate(orbit(0:lat%n_ele_max))
call init_coord(orbit(0), lat%particle_start, lat%ele(0), downstream_end$) 
call track_all(lat, orbit)

n = lat%n_ele_track
write (1, '(a, 4es20.12)') '"Orbit%vec4" ABS 1e-6', orbit(n)%vec(1:4)
write (1, '(a, 2es20.12)') '"Orbit%z"    ABS 1e-6', orbit(n)%vec(5:6)
write (1, '(a, es20.12)')  '"Orbit%p0c"  REL 1e-8', orbit(n)%p0c
write (1, '(a, es20.12)')  '"Orbit%t"    REL 1e-8', orbit(n)%t
write (1, '(a, 3es20.12)') '"Orbit%spin" REL 1e-8', orbit(n)%spin

ele => lat%ele(n)
write (1, '(a, 3es20.12)') '"ele%energy" REL 1e-8', ele%value(E_tot$), ele%value(p0c$), ele%ref_time

end program
