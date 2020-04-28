program converter_test

use bmad

implicit none

type (lat_struct) lat
type (coord_struct) end_orb
real(rp) vec(6), x, y, r, dx_ds, dy_ds
integer i

call bmad_parser('lat.bmad', lat)

do i = 1, 10
  call track1 (lat%particle_start, lat%ele(1), lat%param, end_orb)
  vec = end_orb%vec
  x = vec(1);  y = vec(3)
  r = sqrt(x**2 + y**2)
  dx_ds = (x * vec(2) + y * vec(4)) / (r * (1 + vec(6)))
  dy_ds = (y * vec(2) - x * vec(4)) / (r * (1 + vec(6)))
  write (1, '(4es14.6)') end_orb%charge, end_orb%p0c * (1 + vec(6)), r, dx_ds, dy_ds
enddo

end program
