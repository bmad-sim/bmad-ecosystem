program csr_and_space_charge_test

use beam_mod

!

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (bunch_struct), target :: bunch_init, bunch0, bunch, bunch2
type (coord_struct), allocatable :: centroid(:)
type (coord_struct), pointer :: p0, p, p2
type (beam_init_struct) beam_init
type (branch_struct), pointer :: branch

logical err

!

open (1, file = 'output.now')

call bmad_parser ('csr_and_space_charge_test.bmad', lat)

! Branch 1: e_gun with cathode

branch => lat%branch(1)
call ran_seed_put(1)
beam_init%n_particle = 100
beam_init%a_norm_emit = 1.0e-6
beam_init%b_norm_emit = 1.0e-6
beam_init%dPz_dz = 0.0
beam_init%bunch_charge = 77.0e-12
beam_init%sig_pz = 0e-9
beam_init%sig_z = 0.000899377 ! 3 ps * cLight

call init_bunch_distribution (branch%ele(0), branch%param, beam_init, 0, bunch_init)

bunch2 = bunch_init
!!call track1_bunch(bunch2, branch%ele(1), err)

! Branch 0: CSR

branch => lat%branch(0)
call ran_seed_put(1)
beam_init%n_particle = 1000
beam_init%a_norm_emit = 1.0e-6
beam_init%b_norm_emit = 1.0e-6
beam_init%dPz_dz = 0.0
beam_init%bunch_charge = 77.0e-12
beam_init%sig_pz = 0e-9
beam_init%sig_z = 0.000899377 ! 3 ps * cLight

call init_bunch_distribution (branch%ele(0), branch%param, beam_init, 0, bunch_init)

!

call reallocate_coord(centroid, lat, branch%ix_branch)
call init_coord (centroid(0), lat%particle_start, branch%ele(0), downstream_end$)
call track_all (lat, centroid, branch%ix_branch)

ele => branch%ele(1)
bunch0 = bunch_init
call track1_bunch(bunch0, ele, err, centroid)
p0 => bunch0%particle(10)
write (1, '(a, 6es16.8)') '"No-CSR-Space-Ch" ABS 1e-14', p0%vec

space_charge_com%beam_chamber_height = 0.0254 !1 inch full height
space_charge_com%n_shield_images = 0
space_charge_com%ds_track_step = 0.01
space_charge_com%n_bin = 40
space_charge_com%particle_bin_span = 2
space_charge_com%lsc_sigma_cutoff = 0.1

!

bmad_com%csr_and_space_charge_on = .true.
ele%csr_method = one_dim$
bunch = bunch_init
call track1_bunch(bunch, ele, err, centroid)
p => bunch%particle(10)
write (1, '(a, 6es16.8)') '"CSR"             ABS 1e-14', p%vec - p0%vec

!

space_charge_com%n_shield_images = 4
space_charge_com%beam_chamber_height = 0.01
bunch2 = bunch_init
call track1_bunch(bunch2, ele, err, centroid)
p2 => bunch2%particle(10)
write (1, '(a, 6es16.8)') '"CSR-with-Shield" ABS 1e-14', p2%vec(1:5) - p%vec(1:5)
write (1, '(a, 6es16.8)') '"CSR-with-Shield" ABS 1e-13', p2%vec(6) - p%vec(6)

!

ele%csr_method = off$
ele%space_charge_method = slice$
bunch = bunch_init
call track1_bunch(bunch, ele, err, centroid)
p => bunch%particle(10)
write (1, '(a, 6es16.8)') '"Space_Charge"    ABS 1e-14', p%vec - p0%vec

end program
