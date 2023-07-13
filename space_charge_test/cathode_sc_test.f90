program cathode_sc_test

use bmad
use space_charge_mod

implicit none

type (lat_struct) :: lat
type (ele_struct), pointer :: ele
type (beam_struct), target :: beam
type (bunch_struct), pointer :: bunch
type (coord_struct), pointer :: p
logical err
integer i


!Parse lattice
call bmad_parser('lat.bmad',lat)
!Read beam
call hdf5_read_beam('initial_particles.h5', beam, err, lat%ele(0))
bunch => beam%bunch(1)

!Track through elements
do i = 1, lat%n_ele_track
  ele => lat%ele(i)
  call track1_bunch_space_charge(bunch, ele, err)
enddo
  
open (1, file = 'output.now')
p => bunch%particle(1)
write (1, '(a, es20.10)') '"cathode_sc:vec(1)" ABS  1e-10', p%vec(1)
write (1, '(a, es20.10)') '"cathode_sc:vec(2)" ABS  1e-09', p%vec(2)
write (1, '(a, es20.10)') '"cathode_sc:vec(3)" ABS  1e-10', p%vec(3)
write (1, '(a, es20.10)') '"cathode_sc:vec(4)" ABS  1e-09', p%vec(4)
write (1, '(a, es20.10)') '"cathode_sc:vec(5)" ABS  2e-09', p%vec(5)
write (1, '(a, es20.10)') '"cathode_sc:vec(6)" ABS  2e-08', p%vec(6)
write (1, '(a, es20.10)') '"cathode_sc:s"      ABS  1e-10', p%s
write (1, '(a, es20.10)') '"cathode_sc:t"      ABS  1e-10', p%t

end program
