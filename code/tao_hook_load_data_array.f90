!+
! Subroutine tao_hook_load_data_array (found, datum, lattice, orb, datum_value)
!-

subroutine tao_hook_load_data_array (found, datum, lattice, orb, datum_value)

use tao_mod

implicit none

type (tao_data_struct) datum
type (ring_struct) lattice
type (coord_struct) orb(0:)

real(rp) datum_value
logical found

!

found = .false.

end subroutine
