


program hdf5_test

use beam_mod
use hdf5_openpmd_mod

implicit none

type (lat_struct) lat
type (bunch_struct), target :: bunch1, bunch2, bunches(3)
type (coord_struct), pointer :: p
type (beam_struct) beam, beam2
type (beam_init_struct) beam_init
type (ele_struct) ele
type (grid_field_struct), pointer :: gf1(:) => null(), gf0(:) => null()
type (pmd_unit_struct) unit

integer i, j, n_part
integer h5_err
integer(HID_T) f_id, r_id, z_id, complex_t
type (hdf5_info_struct) info
character(40) sname
real(rp) q_tot, q_live
real(rp) patch_lo(4,3), patch_hi(4,3), patch_offset(4), patch_extent(4)
integer(int64) id8(20)
character(1), parameter :: axes(3) = ['x', 'y', 'z']

complex(rp) c1_in(2), c1_out(2), c2_in(2,3), c2_out(2,3), c3_in(2,3,4), c3_out(2,3,4)

logical error, good1(10), good2(10), g1, g2, ignore_beta

character(8), allocatable :: strs(:)

!

open (1, file = 'output.now', recl = 200)

!---------------------
! General Write/Read

call hdf5_open_file ('general.h5', 'WRITE', f_id, error)
call pmd_init_compound_complex(complex_t)

call H5Gcreate_f(f_id, '/data', r_id, h5_err) 

call hdf5_write_attribute_string(r_id, 'attrib_name', [character(8):: 'Val1', 'Bval2', 'X'], error)

call hdf5_write_dataset_real(r_id, 'scalar_real', 1.25_rp, error)
call hdf5_write_dataset_int(r_id, 'scalar_int', 17, error)
patch_lo(:,1) = [-2.0_rp, 1.0e8_rp, 0.0_rp, 0.0_rp]
patch_lo(:,2) = [-1.0e8_rp, -1.0_rp, 0.0_rp, 0.0_rp]
patch_lo(:,3) = [-3.0_rp, 4.0_rp, 0.0_rp, 0.0_rp]
patch_hi = patch_lo
patch_hi(1,:) = [2.0_rp, 1.0e8_rp, 3.0_rp]
call pmd_write_particle_patches(r_id, [2, 1, 1, 0], [0, 2, 3, 4], patch_lo, patch_hi, error)

c1_in = [cmplx(1.0_rp, 2.0_rp), cmplx(3.0_rp, 4.0_rp)]
do i = 1, ubound(c2_in,2)
  c2_in(:,i) = i * c1_in
  do j = 1, ubound(c3_in, 3)
    c3_in(:,i,j) = (i + 10*j) * c1_in
  enddo
enddo

call pmd_write_complex_to_dataset(r_id, 'complex_rank1', complex_t, "", unit, c1_in, error)
call pmd_write_complex_to_dataset(r_id, 'complex_rank2', complex_t, "", unit, c2_in, error)
call pmd_write_complex_to_dataset(r_id, 'complex_rank3', complex_t, "", unit, c3_in, error)

call H5Gclose_f(r_id, h5_err)
call H5fclose_f(f_id, h5_err)

!

call hdf5_open_file ('general.h5', 'READ', f_id, error)
r_id = hdf5_open_group(f_id, 'data', error, .true.)

call hdf5_read_attribute_string(r_id, 'attrib_name', strs, error, .true.)

call pmd_read_complex_dataset(r_id, 'complex_rank3', complex_t, 1.0_rp, '', c3_out, error)
call pmd_read_complex_dataset(r_id, 'complex_rank2', complex_t, 1.0_rp, '', c2_out, error)
call pmd_read_complex_dataset(r_id, 'complex_rank1', complex_t, 1.0_rp, c1_out, error)

call hdf5_read_dataset_real(r_id, 'scalar_real', q_tot, error)
call hdf5_read_dataset_int(r_id, 'scalar_int', n_part, error)
write (1, '(a, 2l1, a)') '"Scalar-datasets" STR "', q_tot == 1.25_rp, n_part == 17, '"'
do i = 1, 3
  call pmd_read_real_dataset(r_id, 'particlePatches/offset/' // axes(i), 1.0_rp, patch_offset, error)
  call pmd_read_real_dataset(r_id, 'particlePatches/extent/' // axes(i), 1.0_rp, patch_extent, error)
  good1(i) = .not. error .and. all(patch_offset == patch_lo(:,i)) .and. &
                 all(patch_extent > 0) .and. all(patch_offset + patch_extent > patch_hi(:,i))
enddo
write (1, '(a, 3l1, a)') '"Patches-half-open-bounds" STR "', good1(1:3), '"'

call H5Gclose_f(r_id, h5_err)
call H5fclose_f(f_id, h5_err)

!

write (1, '(a, 10a8)') '"Strs" STR', (quote(strs(i)), i = 1, size(strs))

call pmd_kill_compound_complex(complex_t)

!---------------------
! Beam Read/Write

beam_init%position_file = 'beam_1.dat'
n_part = 10

call bmad_parser('lat.bmad', lat)
call init_bunch_distribution(lat%ele(0), lat%param, beam_init, 0, bunch1)
beam_init%use_t_coords = .true.
beam_init%use_z_as_t = .true.
call init_bunch_distribution(lat%ele(0), lat%param, beam_init, 0, bunch2)

! Some checks that the distribution has been correctly setup

p => bunch1%particle(7)
write (1, '(a, 6es18.10)') '"Bunch1-vec"    REL 1e-10', (p%vec(i), i = 1, 6)
write (1, '(a, 6es18.10)') '"Bunch1-other"  REL 1e-10', p%s, p%t, p%charge, p%p0c, p%beta
write (1, '(a, 6es18.10)') '"Bunch1-spin"   REL 1e-10', (p%spin(i), i = 1, 3)
write (1, '(a, 6i6)')      '"Bunch1-int"    ABS ', p%state, p%direction, p%species, p%location

p => bunch2%particle(4)
write (1, '(a, 6es18.10)') '"Bunch1-vec"    REL 1e-10', (p%vec(i), i = 1, 6)
write (1, '(a, 6es18.10)') '"Bunch1-other"  REL 1e-10', p%s, p%t, p%charge, p%p0c, p%beta
write (1, '(a, 6es18.10)') '"Bunch1-spin"   REL 1e-10', (p%spin(i), i = 1, 3)
write (1, '(a, 6i6)')      '"Bunch1-int"    ABS ', p%state, p%direction, p%species, p%location

!

do i = 1, n_part
  j = i + 2
  if (j > n_part) j = j - n_part
  bunch2%particle(i) = bunch1%particle(j)
  bunch2%particle(i)%species = photon$
  bunch2%particle(i)%spin = 0
  bunch2%particle(i)%field = [i, i+1]
  bunch2%particle(i)%dt_ref = 100_rp * (i+5.0_rp)
  bunch2%particle(i)%phase = 1d-3 * [j, j-1]
  bunch2%particle(i)%beta = 1
  bunch2%particle(i)%s = j
enddo

call hdf5_write_beam ('bunch.h5', [bunch1], .false., error)
call hdf5_write_beam ('bunch.h5', [bunch2], .true., error)
call hdf5_read_beam('bunch.h5', beam, error)

do i = 1, n_part
  good1(i) = coord_is_equal(bunch1%particle(i), beam%bunch(1)%particle(i))
  good2(i) = coord_is_equal(bunch2%particle(i), beam%bunch(2)%particle(i))  
enddo

write (1, '(a, 10l1, a)') '"Bunch1" STR "', good1, '"'
write (1, '(a, 10l1, a)') '"Bunch2" STR "', good2, '"'

!---------------------
! Bunch old format

call hdf5_read_beam('bunch_old.h5', beam, error)

do i = 1, n_part
  good1(i) = coord_is_equal(bunch1%particle(i), beam%bunch(1)%particle(i))
  good2(i) = coord_is_equal(bunch2%particle(i), beam%bunch(2)%particle(i))  
enddo

write (1, '(a, 10l1, a)') '"Old-Bunch1" STR "', good1, '"'
write (1, '(a, 10l1, a)') '"Old-Bunch2" STR "', good2, '"'

!---------------------
! Bunch ASCII

call write_beam_file ('bunch.ascii', beam)
call read_beam_file ('bunch.ascii', beam2, beam_init_struct(), error)

! Ignore beta for bunch2 since using time coords here.

do i = 1, n_part
  good1(i) = coord_is_equal(beam2%bunch(1)%particle(i), beam%bunch(1)%particle(i))
  good2(i) = coord_is_equal(beam2%bunch(2)%particle(i), beam%bunch(2)%particle(i), ignore_beta = .true.)  
enddo

write (1, '(a, 10l1, a)') '"ASCII-Bunch1" STR "', good1, '"'
write (1, '(a, 10l1, a)') '"ASCII-Bunch2" STR "', good2, '"'

!---------------------
! Patches: three bunches as one species record. The middle bunch holds no particles, the
! first has two lost particles that alive_only drops and a pre_born one it keeps, and every
! particle carries a label in ix_user that must come back.

bunches(1) = bunch1
bunches(2) = bunch1
call reallocate_bunch (bunches(2), 0)
bunches(3) = bunch1
do i = 1, n_part
  bunches(1)%particle(i)%ix_user = 100 + i
  bunches(3)%particle(i)%ix_user = 300 + i
enddo
bunches(1)%particle(2)%state = lost$
bunches(1)%particle(5)%state = lost$
bunches(1)%particle(7)%state = pre_born$

call hdf5_write_beam ('patches.h5', bunches, .false., error, alive_only = .true., as_patches = .true.)
call hdf5_read_beam ('patches.h5', beam, error)

write (1, '(a, 4i6)') '"Patches-counts" ABS 0 ', size(beam%bunch), size(beam%bunch(1)%particle), &
                                                 size(beam%bunch(2)%particle), size(beam%bunch(3)%particle)
j = 0
good1 = .true.
do i = 1, n_part
  if (bunches(1)%particle(i)%state == lost$) cycle
  j = j + 1
  good1(j) = coord_is_equal(bunches(1)%particle(i), beam%bunch(1)%particle(j)) .and. &
                                                   (beam%bunch(1)%particle(j)%ix_user == 100 + i)
enddo
good1(n_part-1:n_part) = (j == n_part - 2)
do i = 1, n_part
  good2(i) = coord_is_equal(bunches(3)%particle(i), beam%bunch(3)%particle(i)) .and. &
                                                   (beam%bunch(3)%particle(i)%ix_user == 300 + i)
enddo
write (1, '(a, 10l1, a)') '"Patches-Bunch1" STR "', good1, '"'
write (1, '(a, 10l1, a)') '"Patches-Bunch3" STR "', good2, '"'

! Without the writer's mark the patches are another code's partition and the species is
! one bunch.

call hdf5_open_file ('patches.h5', 'APPEND', f_id, error)
z_id = hdf5_open_group (f_id, 'data/00001/particles', error, .true.)
call hdf5_get_object_by_index (z_id, 0, sname, info, error)
r_id = hdf5_open_group (z_id, trim(sname), error, .true.)
call hdf5_read_attribute_real (r_id, 'totalCharge', q_tot, error, .true.)
call hdf5_read_attribute_real (r_id, 'chargeLive', q_live, error, .true.)
call H5Gclose_f (r_id, h5_err)
write (1, '(a, l1, a)') '"Patches-totalCharge-is-live" STR "', abs(q_tot - q_live) <= 1d-12 * abs(q_live), '"'
r_id = hdf5_open_group (z_id, trim(sname) // '/particlePatches', error, .true.)
call H5Adelete_f (r_id, 'patchesAreBunches', h5_err)
call H5Gclose_f (r_id, h5_err)
call H5Gclose_f (z_id, h5_err)
call H5Fclose_f (f_id, h5_err)
call hdf5_read_beam ('patches.h5', beam, error)
write (1, '(a, 2i6)') '"Patches-unmarked" ABS 0 ', size(beam%bunch), size(beam%bunch(1)%particle)

! An ion beam whose bunches hold two charge states of one species is one species record,
! each particle's charge state in chargeState.

bunches(1)%particle%state = alive$
bunches(1)%particle%species = species_id('He+')
bunches(3)%particle%species = species_id('He++')
call hdf5_write_beam ('patches_ion.h5', bunches, .false., error, as_patches = .true.)
write (1, '(a, l1, a)') '"Patches-ion-written" STR "', .not. error, '"'
call hdf5_read_beam ('patches_ion.h5', beam, error)
write (1, '(a, 3i6)') '"Patches-ion-counts" ABS 0 ', size(beam%bunch(1)%particle), &
                                     size(beam%bunch(2)%particle), size(beam%bunch(3)%particle)
write (1, '(a, l1, l1, a)') '"Patches-ion-species" STR "', all(beam%bunch(1)%particle%species == species_id('He+')), &
                                     all(beam%bunch(3)%particle%species == species_id('He++')), '"'

! An id beyond the default-integer ix_user, as a 64-bit writer may store it, is clamped
! with a warning.

call hdf5_open_file ('patches_ion.h5', 'APPEND', f_id, error)
z_id = hdf5_open_group (f_id, 'data/00001/particles/He', error, .true.)
call H5Ldelete_f (z_id, 'id', h5_err)
id8 = [(int(i, int64), i = 1, 20)]
id8(5) = 2_int64**40
call h5ltmake_dataset_f (z_id, 'id', 1, [20_hsize_t], h5kind_to_type(int64, H5_INTEGER_KIND), id8, h5_err)
call H5Gclose_f (z_id, h5_err)
call H5Fclose_f (f_id, h5_err)
call hdf5_read_beam ('patches_ion.h5', beam, error)
write (1, '(a, 3i12)') '"Patches-id64" ABS 0 ', beam%bunch(1)%particle(4)%ix_user, &
                                     beam%bunch(1)%particle(5)%ix_user, beam%bunch(1)%particle(6)%ix_user

bunches(1) = bunch1
bunches(1)%particle%state = alive$
bunches(1)%particle%ix_user = -1
call hdf5_write_beam ('ids.h5', [bunches(1)], .false., error)
call hdf5_open_file ('ids.h5', 'APPEND', f_id, error)
sname = openpmd_species_name(set_species_charge(bunch1%particle(1)%species, 0))
z_id = hdf5_open_group (f_id, 'data/00001/particles/' // trim(sname), error, .true.)
g1 = hdf5_exists(z_id, 'id', error, .false.)
write (1, '(a, l1, a)') '"Unset-id-omitted" STR "', .not. g1, '"'
call H5Gclose_f (z_id, h5_err)
call H5Fclose_f (f_id, h5_err)
call hdf5_read_beam ('ids.h5', beam, error)
write (1, '(a, l1, a)') '"Absent-id-default" STR "', &
               .not. error .and. all(beam%bunch(1)%particle%ix_user == -1), '"'

call hdf5_open_file ('ids.h5', 'APPEND', f_id, error)
z_id = hdf5_open_group (f_id, 'data/00001/particles/' // trim(sname), error, .true.)
call pmd_write_int_to_dataset(z_id, 'id', 'ix_user', unit_1, [(42, i = 1, n_part)], error)
r_id = hdf5_open_group (z_id, 'id', error, .true.)
call H5Adelete_f (r_id, 'unitSI', h5_err)
call H5Gclose_f (r_id, h5_err)
call H5Gclose_f (z_id, h5_err)
call H5Fclose_f (f_id, h5_err)
call hdf5_read_beam ('ids.h5', beam, error)
write (1, '(a, l1, a)') '"Constant-id-without-unitSI" STR "', &
               .not. error .and. all(beam%bunch(1)%particle%ix_user == 42), '"'

call hdf5_open_file ('ids.h5', 'APPEND', f_id, error)
z_id = hdf5_open_group (f_id, 'data/00001/particles/' // trim(sname), error, .true.)
r_id = hdf5_open_group (z_id, 'id', error, .true.)
call hdf5_write_attribute_int(r_id, 'shape', [n_part - 1], error)
call H5Gclose_f (r_id, h5_err)
call H5Gclose_f (z_id, h5_err)
call H5Fclose_f (f_id, h5_err)
call hdf5_read_beam ('ids.h5', beam, error)
write (1, '(a, l1, a)') '"Constant-id-bad-shape" STR "', error, '"'

call hdf5_open_file ('ids.h5', 'APPEND', f_id, error)
z_id = hdf5_open_group (f_id, 'data/00001/particles/' // trim(sname), error, .true.)
call H5Ldelete_f (z_id, 'id', h5_err)
call h5ltmake_dataset_f (z_id, 'id', 1, [9_hsize_t], h5kind_to_type(int64, H5_INTEGER_KIND), id8(1:9), h5_err)
call H5Gclose_f (z_id, h5_err)
call H5Fclose_f (f_id, h5_err)
call hdf5_read_beam ('ids.h5', beam, error)
write (1, '(a, l1, a)') '"Short-id-rejected" STR "', error, '"'

bunches(1)%particle(1)%state = lost$
bunches(1)%particle(1)%species = photon$
call hdf5_write_beam ('filtered.h5', [bunches(1)], .false., error, alive_only = .true., as_patches = .true.)
write (1, '(a, l1, a)') '"Patches-filtered-species-written" STR "', .not. error, '"'
call hdf5_read_beam ('filtered.h5', beam, error)
write (1, '(a, l1, a)') '"Patches-filtered-species-read" STR "', .not. error .and. &
  size(beam%bunch(1)%particle) == n_part - 1 .and. &
  all(beam%bunch(1)%particle%species == bunch1%particle(1)%species), '"'

bunches(1)%particle%state = lost$
call hdf5_write_beam ('filtered.h5', [bunches(1)], .false., error, alive_only = .true., as_patches = .true.)
write (1, '(a, l1, a)') '"Patches-all-filtered-written" STR "', .not. error, '"'
call hdf5_read_beam ('filtered.h5', beam, error)
write (1, '(a, l1, a)') '"Patches-all-filtered-read" STR "', &
               .not. error .and. size(beam%bunch(1)%particle) == 0, '"'

bunches = bunch1
call hdf5_write_beam ('constant_patches.h5', bunches, .false., error, as_patches = .true.)
call hdf5_open_file ('constant_patches.h5', 'APPEND', f_id, error)
z_id = hdf5_open_group (f_id, 'data/00001/particles/' // trim(sname) // '/particlePatches', error, .true.)
call H5Ldelete_f (z_id, 'numParticles', h5_err)
call pmd_write_int_to_dataset(z_id, 'numParticles', 'count', unit_1, [n_part, n_part, n_part], error)
call H5Gclose_f (z_id, h5_err)
call H5Fclose_f (f_id, h5_err)
call hdf5_read_beam ('constant_patches.h5', beam, error)
write (1, '(a, l1, a)') '"Patches-constant-counts" STR "', .not. error .and. size(beam%bunch) == 3 .and. &
               all([(size(beam%bunch(i)%particle) == n_part, i = 1, 3)]), '"'

call hdf5_open_file ('constant_patches.h5', 'APPEND', f_id, error)
z_id = hdf5_open_group (f_id, &
  'data/00001/particles/' // trim(sname) // '/particlePatches/numParticles', error, .true.)
call H5Adelete_f (z_id, 'shape', h5_err)
call H5Gclose_f (z_id, h5_err)
call H5Fclose_f (f_id, h5_err)
call hdf5_read_beam ('constant_patches.h5', beam, error)
write (1, '(a, l1, a)') '"Patches-missing-shape" STR "', error, '"'

!---------------------
! A failed read must close the file. Each case corrupts a fresh file, reads it, and reopens it.

sname = openpmd_species_name(set_species_charge(bunch1%particle(1)%species, 0))

call leak_file_setup ('leak1.h5')
call H5Adelete_f (f_id, 'basePath', h5_err)
call H5Fclose_f (f_id, h5_err)
call leak_check ('leak1.h5', 'Leak-header')

call leak_file_setup ('leak2.h5')
call H5Ldelete_f (f_id, 'data', h5_err)
call H5Fclose_f (f_id, h5_err)
call leak_check ('leak2.h5', 'Leak-base-group')

call leak_file_setup ('leak3.h5')
call H5Adelete_f (f_id, 'basePath', h5_err)
call hdf5_write_attribute_string (f_id, 'basePath', '/missing/', error)
call H5Fclose_f (f_id, h5_err)
call leak_check ('leak3.h5', 'Leak-single-bunch-group')

call leak_file_setup ('leak4.h5')
call H5Adelete_f (f_id, 'particlesPath', h5_err)
call hdf5_write_attribute_string (f_id, 'particlesPath', 'missing/', error)
call H5Fclose_f (f_id, h5_err)
call leak_check ('leak4.h5', 'Leak-particles-group')

call leak_file_setup ('leak5.h5')
z_id = hdf5_open_group (f_id, 'data/00001/particles/' // trim(sname), error, .true.)
call H5Adelete_f (z_id, 'numParticles', h5_err)
call H5Gclose_f (z_id, h5_err)
call H5Fclose_f (f_id, h5_err)
call leak_check ('leak5.h5', 'Leak-no-species')

call leak_file_setup ('leak6.h5')
z_id = hdf5_open_group (f_id, 'data/00001/particles/' // trim(sname), error, .true.)
call H5Adelete_f (z_id, 'totalCharge', h5_err)
call hdf5_write_attribute_string (z_id, 'totalCharge', 'bad', error)
call H5Gclose_f (z_id, h5_err)
call H5Fclose_f (f_id, h5_err)
call leak_check ('leak6.h5', 'Leak-species-attribute')

!---------------------
! Grid_Field Read/Write

ele%key = sbend$
ele%value(rho$) = 4
ele%value(l$) = 3

allocate (gf0(2))
allocate (gf0(1)%ptr, gf0(2)%ptr)
allocate (gf0(1)%ptr%pt(4,2,3), gf0(2)%ptr%pt(2,3,1))

call grid_field_set_values (gf0(1), 0, xyz$)
call grid_field_set_values (gf0(2), 1, rotationally_symmetric_rz$)

call hdf5_write_grid_field ('grid_field.h5', ele, gf0, error)
call hdf5_read_grid_field  ('grid_field.h5', ele, gf1, error)

write (1, '(a, l1, a)') '"Grid1" STR "', grid_field_is_equal (gf0(1), gf1(1)), '"'
write (1, '(a, l1, a)') '"Grid2" STR "', grid_field_is_equal (gf0(2), gf1(2)), '"'

!---------------------
! Old Grid_Field Read

call hdf5_read_grid_field  ('grid_field_old.h5', ele, gf1, error)

write (1, '(a, l1, a)') '"Old-Grid1" STR "', grid_field_is_equal (gf0(1), gf1(1)), '"'
write (1, '(a, l1, a)') '"Old-Grid2" STR "', grid_field_is_equal (gf0(2), gf1(2)), '"'

close (1)

!-----------------------------------------
contains

subroutine leak_file_setup (file)
character(*) file
call hdf5_write_beam (file, [bunch1], .false., error)
call hdf5_open_file (file, 'APPEND', f_id, error)
end subroutine leak_file_setup

!-----------------------------------------

subroutine leak_check (file, label)
character(*) file, label
integer(HID_T) fid
logical read_err, open_err
call hdf5_read_beam (file, beam, read_err)
call hdf5_open_file (file, 'APPEND', fid, open_err)
if (.not. open_err) call H5Fclose_f (fid, h5_err)
write (1, '(3a, 2l1, a)') '"', label, '" STR "', read_err, .not. open_err, '"'
end subroutine leak_check

!-----------------------------------------

function coord_is_equal(p1, p2, ignore_beta) result (equal)

type (coord_struct) p1, p2
logical equal
logical, optional :: ignore_beta

equal = is_eq_rv(p1%vec, p2%vec, 'bunch vec')
equal = equal .and. is_eq_rv(p1%spin, p2%spin, 'bunch spin')
equal = equal .and. is_eq_rv(p1%field, p2%field, 'bunch field')
equal = equal .and. is_eq_rv(p1%phase, p2%phase, 'bunch phase')
equal = equal .and. is_eq_r(p1%s, p2%s, 'bunch s')
equal = equal .and. is_eq_r(real(p1%t, rp), real(p2%t, rp), 'bunch t')
equal = equal .and. is_eq_r(p1%charge, p2%charge, 'bunch charge')
equal = equal .and. is_eq_r(p1%dt_ref, p2%dt_ref, 'bunch dt_ref')
equal = equal .and. is_eq_r(p1%p0c, p2%p0c, 'bunch p0c')
if (.not. logic_option(.false., ignore_beta)) then
  equal = equal .and. is_eq_r(p1%beta, p2%beta, 'bunch beta')
endif

equal = equal .and. is_eq_i(p1%ix_ele, p2%ix_ele, 'bunch ix_ele')
equal = equal .and. is_eq_i(p1%ix_branch, p2%ix_branch, 'bunch ix_branch')
equal = equal .and. is_eq_i(p1%state, p2%state, 'bunch state')
equal = equal .and. is_eq_i(p1%direction, p2%direction, 'bunch direction')
equal = equal .and. is_eq_i(p1%species, p2%species, 'bunch species')
equal = equal .and. is_eq_i(p1%location, p2%location, 'bunch location')

end function coord_is_equal

!-----------------------------------------
! contains

function is_eq_r (r1, r2, who) result (equal)
real(rp), intent(in) :: r1, r2
character(*) who
logical equal

equal = (abs(r1-r2) <= 1d-15 * (abs(r1) + abs(r2)))
if (.not. equal) then
  print *, trim(who) // ' Not equal!'
endif

end function is_eq_r

!-----------------------------------------
! contains

function is_eq_rv (r1, r2, who) result (equal)
real(rp), intent(in) :: r1(:), r2(:)
character(*) who
logical equal

equal = all(abs(r1-r2) <= 1d-15 * (abs(r1) + abs(r2)))
if (.not. equal) then
  print *, trim(who) // ' Not equal!'
endif

end function is_eq_rv

!-----------------------------------------
! contains

function is_eq_i (i1, i2, who) result (equal)
integer, intent(in) :: i1, i2
character(*) who
logical equal

equal = (i1 == i2)
if (.not. equal) then
  print *, trim(who) // ' Not equal!'
endif

end function is_eq_i

!-----------------------------------------
! contains

function is_eq_iv (i1, i2, who) result (equal)
integer, intent(in) :: i1(:), i2(:)
character(*) who
logical equal

equal = all(i1 == i2)
if (.not. equal) then
  print *, trim(who) // ' Not equal!'
endif

end function is_eq_iv

!-----------------------------------------
! contains

function is_eq_l  (l1, l2, who) result (equal)
logical, intent(in) :: l1, l2
character(*) who
logical equal

equal = (l1 .eqv. l2)
if (.not. equal) then
  print *, trim(who) // ' Not equal!'
endif

end function is_eq_l

!-----------------------------------------------------------------
! contains

subroutine grid_field_set_values (gf, ixg, geometry)

type (grid_field_struct), target :: gf
type (grid_field_pt1_struct), pointer :: pt(:,:,:)
integer i, j, k, ixg, geometry, base

!

pt => gf%ptr%pt

do i = 1, size(pt, 1)
do j = 1, size(pt, 2)
do k = 1, size(pt, 3)
  base = 10000*i + 1000*j + 100*k + 10*ixg 
  pt(i,j,k)%B = cmplx(base + [1, 2, 3], 0.0_rp)
  pt(i,j,k)%E = cmplx(0.0_rp, base + [4, 5, 6])
enddo
enddo
enddo

gf%geometry = geometry
gf%harmonic = ixg
gf%phi0_fieldmap = 0
gf%field_scale = ixg + 4
gf%field_type = mixed$
gf%master_parameter = ixg
gf%ele_anchor_pt = ixg + 1
gf%interpolation_order = 2*ixg + 1
gf%dr = [ixg+1, ixg+2, ixg+3]
gf%r0 = [ixg+4, ixg+5, ixg+6]
gf%curved_ref_frame = (ixg == 0)

end subroutine grid_field_set_values

!-----------------------------------------------------------------
! contains

function grid_field_is_equal (gf0, gf1) result (is_equal)

type (grid_field_struct) gf0, gf1
type (grid_field_pt1_struct), pointer :: pt0(:,:,:), pt1(:,:,:)
integer i, j, k
logical is_equal

!

is_equal = .true.

pt0 => gf0%ptr%pt
pt1 => gf1%ptr%pt

is_equal = is_equal .and. is_eq_iv(shape(pt0), shape(pt1), 'Grid_field Shape')
if (.not. is_equal) return

do i = 1, 3
  if (any(int(pt0%B(i)) /= int(pt1%B(i)))) then
    print '(a, i0, a)', 'B', i, ': Not Same' 
    is_equal = .false.
  endif
  if (any(int(pt0%E(i)) /= int(pt1%E(i)))) then
    print '(a, i0, a)', 'E', i, ': Not Same' 
    is_equal = .false.
  endif
enddo

is_equal = is_equal .and. is_eq_i(gf0%geometry, gf1%geometry, 'Grid_field geometry')
is_equal = is_equal .and. is_eq_i(gf0%harmonic, gf1%harmonic, 'Grid_field harmonic')
is_equal = is_equal .and. is_eq_r(gf0%phi0_fieldmap, gf1%phi0_fieldmap, 'Grid_field phi0_fieldmap')
is_equal = is_equal .and. is_eq_r(gf0%field_scale, gf1%field_scale, 'Grid_field field_scale')
is_equal = is_equal .and. is_eq_i(gf0%field_type, gf1%field_type, 'Grid_field field_type')
is_equal = is_equal .and. is_eq_i(gf0%master_parameter, gf1%master_parameter, 'Grid_field master_parameter')
is_equal = is_equal .and. is_eq_i(gf0%ele_anchor_pt, gf1%ele_anchor_pt, 'Grid_field ele_anchor_pt')
is_equal = is_equal .and. is_eq_i(gf0%interpolation_order, gf1%interpolation_order, 'Grid_field interpolation_order')
is_equal = is_equal .and. is_eq_rv(gf0%dr, gf1%dr, 'Grid_field dr')
is_equal = is_equal .and. is_eq_rv(gf0%r0, gf1%r0, 'Grid_field r0')
is_equal = is_equal .and. is_eq_l(gf0%curved_ref_frame, gf1%curved_ref_frame, 'Grid_field curved_ref_frame')

end function grid_field_is_equal 

end program
