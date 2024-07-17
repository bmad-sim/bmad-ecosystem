program nonlin_test

use bmad

implicit none

type (lat_struct) lat
type (coord_struct), allocatable :: orb(:)

real(rp) ch_x, ch_y
integer status
logical err

!

open (1, file = 'output.now', recl = 200)

bmad_com%auto_bookkeeper = .false.
call bmad_parser('small_ring.bmad', lat)
call twiss_and_track (lat, orb, status)

call chrom_tune (lat, 1.0e-4_rp, 1.0_rp, 2.0_rp, 0.05_rp, err)
call chrom_calc (lat, 1.0e-4_rp, ch_x, ch_y, err)

write (1, '(a, 2f14.8)') '"Chrom" ABS 1e-7', ch_x, ch_y

close(1)

end program
