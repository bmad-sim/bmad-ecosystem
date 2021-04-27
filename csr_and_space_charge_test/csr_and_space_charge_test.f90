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

! Branch 1: e_gun

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

call track1_bunch(bunch_init, lat, branch%ele(1), bunch2, err)

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
call track1_bunch(bunch_init, lat, ele, bunch0, err, centroid)
p0 => bunch0%particle(10)
write (1, '(a, 6es16.8)') '"No-CSR-Space-Ch" ABS 1e-14', p0%vec

csr_param%beam_chamber_height = 0.0254 !1 inch full height
csr_param%n_shield_images = 0
csr_param%ds_track_step = 0.01
csr_param%n_bin = 40
csr_param%particle_bin_span = 2
csr_param%sigma_cutoff = 0.1

!

bmad_com%csr_and_space_charge_on = .true.
ele%csr_method = one_dim$
call track1_bunch(bunch_init, lat, ele, bunch, err, centroid)
p => bunch%particle(10)
write (1, '(a, 6es16.8)') '"CSR"             ABS 1e-14', p%vec - p0%vec

!

csr_param%n_shield_images = 4
csr_param%beam_chamber_height = 0.01
call track1_bunch(bunch_init, lat, ele, bunch2, err, centroid)
p2 => bunch2%particle(10)
write (1, '(a, 6es16.8)') '"CSR-with-Shield" ABS 1e-14', p2%vec - p%vec

!

ele%csr_method = off$
ele%space_charge_method = slice$
call track1_bunch(bunch_init, lat, ele, bunch, err, centroid)
p => bunch%particle(10)
write (1, '(a, 6es16.8)') '"Space_Charge"    ABS 1e-14', p%vec - p0%vec

end program
