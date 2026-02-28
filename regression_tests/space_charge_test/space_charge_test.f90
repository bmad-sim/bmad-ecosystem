program cathode_sc_test

use bmad
use space_charge_mod
use beam_utils, only: init_beam_distribution
use beam_mod, only: track1_bunch

implicit none

type (lat_struct) :: lat
type (ele_struct), pointer :: ele
type (beam_struct), target :: beam
type (bunch_struct), pointer :: bunch
type (coord_struct), pointer :: p
type (coord_struct) :: p0, delp
type (beam_init_struct) beam_init
type (coord_struct), allocatable :: centroid(:)

real(rp) sigma(2), nk(2), dnk(2,2)
logical err
integer i

!

open (1, file = 'output.now')

! bbi_kick test

sigma = [10.0_rp, 10.0_rp]
call bbi_kick(0.01_rp, 0.02_rp, sigma, nk, dnk)
write (1, '(a, 2es20.10)') '"bbi_kick1-rnd-nk"   ABS 1e-10', nk
write (1, '(a, 4es20.10)') '"bbi_kick1-rnd-dnk"  ABS 1e-10', dnk

sigma = [10.0_rp, 10.1_rp]
call bbi_kick(0.01_rp, 0.02_rp, sigma, nk, dnk, .false.)
write (1, '(a, 2es20.10)') '"bbi_kick1-nk"   ABS 1e-10', nk
write (1, '(a, 4es20.10)') '"bbi_kick1-dnk"  ABS 1e-10', dnk


sigma = [10.0_rp, 10.1_rp]
call bbi_kick(20.0_rp, 40.0_rp, sigma, nk, dnk, .false.)
write (1, '(a, 2es20.10)') '"bbi_kick2-nk"   ABS 1e-10', nk
write (1, '(a, 4es20.10)') '"bbi_kick2-dnk"  ABS 1e-10', dnk

sigma = [10.0_rp, 10.0_rp]
call bbi_kick(20.0_rp, 40.0_rp, sigma, nk, dnk)
write (1, '(a, 2es20.10)') '"bbi_kick2-rnd-nk"   ABS 1e-10', nk
write (1, '(a, 4es20.10)') '"bbi_kick2-rnd-dnk"  ABS 1e-10', dnk


sigma = [10.0_rp, 10.1_rp]
call bbi_kick(200.0_rp, 400.0_rp, sigma, nk, dnk, .false.)
write (1, '(a, 2es20.10)') '"bbi_kick3-nk"   ABS 1e-10', nk
write (1, '(a, 4es20.10)') '"bbi_kick3-dnk"  ABS 1e-10', dnk

sigma = [10.0_rp, 10.0_rp]
call bbi_kick(200.0_rp, 400.0_rp, sigma, nk, dnk)
write (1, '(a, 2es20.10)') '"bbi_kick3-rnd-nk"   ABS 1e-10', nk
write (1, '(a, 4es20.10)') '"bbi_kick3-rnd-dnk"  ABS 1e-10', dnk


sigma = [20.0_rp, 10.0_rp]
call bbi_kick(0.01_rp, 0.04_rp, sigma, nk, dnk)
write (1, '(a, 2es20.10)') '"bbi_kick4-nk"   ABS 1e-10', nk
write (1, '(a, 4es20.10)') '"bbi_kick4-dnk"  ABS 1e-10', dnk

sigma = [20.0_rp, 10.1_rp]
call bbi_kick(0.01_rp, 0.04_rp, sigma, nk, dnk, .true.)
write (1, '(a, 2es20.10)') '"bbi_kick4-lin-nk"   ABS 1e-10', nk
write (1, '(a, 4es20.10)') '"bbi_kick4-lin-dnk"  ABS 1e-10', dnk

!sigma = [2.0_rp, 1.0_rp]
!call bbi_kick(1.0_rp, 0.5_rp, sigma, nk, dnk)
!write (1, '(a, 2es20.10)') '"bbi_kick-nk"   ABS 1e-10', nk
!write (1, '(a, 4es20.10)') '"bbi_kick-dnk"  ABS 1e-10', dnk

!sigma = [2.0_rp, 1.0_rp]
!call bbi_kick(1.0_rp, 0.5_rp, sigma, nk, dnk, .true.)
!write (1, '(a, 2es20.10)') '"bbi_kick-nk"   ABS 1e-10', nk
!write (1, '(a, 4es20.10)') '"bbi_kick-dnk"  ABS 1e-10', dnk




! Gun test

call bmad_parser('gun.bmad', lat)
call hdf5_read_beam('initial_particles.h5', beam, err, lat%ele(0))
bunch => beam%bunch(1)

do i = 1, lat%n_ele_track
  ele => lat%ele(i)
  call track1_bunch_space_charge(bunch, ele, err)
enddo
  
p => bunch%particle(1)
write (1, '(a, es20.10)') '"cathode_sc:vec(1)" ABS  2e-07', p%vec(1)
write (1, '(a, es20.10)') '"cathode_sc:vec(2)" ABS  2e-07', p%vec(2)
write (1, '(a, es20.10)') '"cathode_sc:vec(3)" ABS  2e-07', p%vec(3)
write (1, '(a, es20.10)') '"cathode_sc:vec(4)" ABS  2e-07', p%vec(4)
write (1, '(a, es20.10)') '"cathode_sc:vec(5)" ABS  5e-07', p%vec(5)
write (1, '(a, es20.10)') '"cathode_sc:vec(6)" ABS  3e-06', p%vec(6)
write (1, '(a, es20.10)') '"cathode_sc:s"      ABS  1e-10', p%s
write (1, '(a, es20.10)') '"cathode_sc:t"      ABS  1e-10', p%t

! Grid test

call bmad_parser('grid.bmad', lat)

call twiss_and_track(lat, centroid, orb_start = lat%particle_start)

space_charge_com%n_bin = 6
space_charge_com%ds_track_step = 0.1
space_charge_com%space_charge_mesh_size = [32, 32, 32]

beam_init%position_file = 'particles_grid.h5'
beam_init%n_particle = 100
call init_beam_distribution(lat%ele(0), lat%param, beam_init, beam)

p0 = bunch%particle(10)
p => bunch%particle(10)
call track1_bunch(bunch, lat%ele(1), err, centroid)
delp%vec = p%vec - p0%vec

write (1, '(a, es20.10)') '"grid:vec(1)" ABS  2e-12', delp%vec(1)
write (1, '(a, es20.10)') '"grid:vec(2)" ABS  2e-12', delp%vec(2)
write (1, '(a, es20.10)') '"grid:vec(3)" ABS  2e-12', delp%vec(3)
write (1, '(a, es20.10)') '"grid:vec(4)" ABS  2e-12', delp%vec(4)
write (1, '(a, es20.10)') '"grid:vec(5)" ABS  3e-12', delp%vec(5)
write (1, '(a, es20.10)') '"grid:vec(6)" ABS  3e-12', delp%vec(6)
write (1, '(a, es20.10)') '"grid:s"      ABS  1e-10', p%s
write (1, '(a, es20.10)') '"grid:t"      ABS  1e-17', p%t

end program
