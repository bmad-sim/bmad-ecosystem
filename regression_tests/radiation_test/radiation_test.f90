program radiation_test

use bmad
use ptc_map_with_radiation_mod
use ptc_layout_mod
use rad_6d_mod
use radiation_mod

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele
type (coord_struct) orb_start, orb_end, orb0, orb
type (coord_struct), allocatable :: orbit(:)
type (normal_modes_struct) mode
type (rad_int_all_ele_struct), target :: ri_cache, ri_no_cache
type (rad_int1_struct), pointer :: rie1(:), rie2(:), rie3(:), rie4(:)
type (ptc_rad_map_struct) rad_map

real(rp) vec6(6), emit(3), sigma_mat(6,6), m(6,6), vec0(6)

integer i, j, n, ib, ix_cache, n1, n2, ie
logical err
character(100) file_name

!

open (1, file = 'output.now', recl = 200)

!

file_name = 'sigma.bmad'
call bmad_parser(file_name, lat)
do i = 1, lat%n_ele_max
  if (associated(lat%ele(i)%rad_map)) lat%ele(i)%rad_map%stale = .true.
enddo

ptc_com%vertical_kick = 0

call emit_6d(lat%ele(0), .true., mode, sigma_mat)

write (1, '(a, 3es14.6)') '"emit_6d" REL 1e-6', mode%a%emittance, mode%b%emittance, mode%z%emittance
do i = 1, 6
  write (1, '(a, i0, a, 6es14.6)') '"sig_mat', i, '" REL 1e-6', sigma_mat(i,:)
enddo

ele => lat%ele(1)
call radiation_map_setup(ele, err)
call write_rad_map (ele%rad_map%rm0, trim(ele%name) // '-rm0')
call write_rad_map (ele%rad_map%rm1, trim(ele%name) // '-rm1')

ele => lat%ele(2)
call radiation_map_setup(ele, err)
call write_rad_map (ele%rad_map%rm0, trim(ele%name) // '-rm0')
call write_rad_map (ele%rad_map%rm1, trim(ele%name) // '-rm1')

!

bmad_com%radiation_damping_on = .false.
call bmad_parser('radiation_test.bmad', lat)

branch => lat%branch(1)
call reallocate_coord(orbit, branch%n_ele_max)
vec6 = 0
do i = 0, branch%n_ele_max
  call init_coord(orbit(i), vec6, branch%ele(i), upstream_end$)
enddo

call twiss_propagate_all(lat, 1)
call twiss_propagate_all(lat, 2)

ix_cache = 0
call radiation_integrals (lat, orbit, mode, ix_cache, 1, ri_cache)
call radiation_integrals (lat, orbit, mode, -2, 1, ri_no_cache)

ix_cache = 0
call radiation_integrals (lat, orbit, mode, ix_cache, 2, ri_cache)
call radiation_integrals (lat, orbit, mode, -2, 2, ri_no_cache)

rie1 => ri_cache%branch(1)%ele
rie2 => ri_no_cache%branch(1)%ele
rie3 => ri_cache%branch(2)%ele
rie4 => ri_no_cache%branch(2)%ele

do i = 1, branch%n_ele_max
  write (1, '(a, 4es14.6)') '"i0-' // int_str(i) // '" REL 1E-6', rie1(i)%i0, rie2(i)%i0, rie3(i)%i0, rie4(i)%i0
  write (1, '(a, 4es14.6)') '"i1-' // int_str(i) // '" REL 1E-6', rie1(i)%i1, rie2(i)%i1, rie3(i)%i1, rie4(i)%i1
  write (1, '(a, 4es14.6)') '"i2-' // int_str(i) // '" REL 1E-6', rie1(i)%i2, rie2(i)%i2, rie3(i)%i2, rie4(i)%i2
  write (1, '(a, 4es14.6)') '"i3-' // int_str(i) // '" REL 1E-6', rie1(i)%i3, rie2(i)%i3, rie3(i)%i3, rie4(i)%i3
  write (1, '(a, 4es14.6)') '"i4a-' // int_str(i) // '" REL 1E-6', rie1(i)%i4a, rie2(i)%i4a, rie3(i)%i4a, rie4(i)%i4a
  write (1, '(a, 4es14.6)') '"i4b-' // int_str(i) // '" REL 1E-6', rie1(i)%i4b, rie2(i)%i4b, rie3(i)%i4b, rie4(i)%i4b
  write (1, '(a, 4es14.6)') '"i5a-' // int_str(i) // '" REL 1E-6', rie1(i)%i5a, rie2(i)%i5a, rie3(i)%i5a, rie4(i)%i5a
  write (1, '(a, 4es14.6)') '"i5b-' // int_str(i) // '" REL 1E-6', rie1(i)%i5b, rie2(i)%i5b, rie3(i)%i5b, rie4(i)%i5b
  write (1, *)
enddo

!

ele => lat%ele(1)

call init_coord (orb_start, lat%particle_start, ele, upstream_end$)

call track1 (orb_start, ele, lat%param, orb_end)
write (1, '(a, 7es16.8)') '"No-Rad"   ABS 1e-16', orb_end%vec

bmad_com%radiation_damping_on = .true.
call track1 (orb_start, ele, lat%param, orb_end)
write (1, '(a, 7es16.8)') '"Rad-Damp" ABS 1e-16', orb_end%vec

bmad_com%radiation_zero_average = .true.
call track1 (orb_start, ele, lat%param, orb_end)
write (1, '(a, 7es16.8)') '"Rad-Zero" ABS 1e-16', orb_end%vec

close(1)

!---------------------------
contains

subroutine write_rad_map(rm, str)

type (rad_map_struct) :: rm
integer i
character(*) str

!

write (1, *)
write (1, '(2a, 6es14.6)') quote(str // '-ref_orb'), ' ABS 1e-12', rm%ref_orb

do i = 1, 6
  write (1, '(2a, 6es14.6)') quote(str // '-damp_dmat:' // int_str(i)), ' ABS 1e-12', rm%damp_dmat(i,:)
enddo

write (1, '(2a, 6es14.6)') quote(str // '-xfer_damp_vec'), ' ABS 1e-12', rm%xfer_damp_vec

do i = 1, 6
  write (1, '(2a, 6es14.6)') quote(str // '-xfer_damp_mat:' // int_str(i)), ' ABS 1e-12', rm%xfer_damp_mat(i,:)
enddo

do i = 1, 6
  write (1, '(2a, 6es14.6)') quote(str // '-stoc_mat:' // int_str(i)), ' ABS 1e-12', rm%stoc_mat(i,:)
enddo

end subroutine write_rad_map

end program
