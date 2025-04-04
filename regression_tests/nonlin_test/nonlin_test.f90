program nonlin_test

use bmad

implicit none

type (lat_struct), target :: lat
type (coord_struct), allocatable :: orb(:)
type (ele_struct), pointer :: ele
type (ele_struct) ele0

real(rp) ch_x, ch_y
integer status
logical err

!

open (1, file = 'output.now', recl = 200)

bmad_com%auto_bookkeeper = .false.
call bmad_parser('small_ring.bmad', lat)
ele0 = lat%ele(0)

call twiss_and_track (lat, orb, status)
call chrom_tune (lat, 1.0e-4_rp, 1.0_rp, 2.0_rp, 0.05_rp, err)
call chrom_calc (lat, 1.0e-4_rp, ch_x, ch_y, err)

ele => lat%ele(2)
write (1, '(a, 2f14.8)') '"Chrom" ABS 1e-7', ch_x, ch_y
write (1, '(a, 4f14.8)') '"W-func-closed" REL 1e-4', ele%a%dbeta_dpz, ele%b%dbeta_dpz, ele%a%dalpha_dpz, ele%b%dalpha_dpz

lat%param%geometry = open$

call chrom_calc (lat, 1.0e-4_rp, ch_x, ch_y, err, orb0 = orb(0))
write (1, '(a, 4f14.8)') '"W-func-open" REL 1e-4', ele%a%dbeta_dpz, ele%b%dbeta_dpz, ele%a%dalpha_dpz, ele%b%dalpha_dpz

close(1)

end program
