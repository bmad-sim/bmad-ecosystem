program geometry_test

use bmad

implicit none

type (lat_struct), target :: lat
type (coord_struct) orbit
type (floor_position_struct) local, floor, pos_ele, floor2, f0
type (ele_struct), pointer :: ele, ele1
type (em_field_struct) f1, f2
type (branch_struct), pointer :: branch

real(rp) w_mat(3,3), s

integer i, status, nargs
logical print_extra

character(100) lat_file

!

lat_file = 'geometry_test.bmad'
print_extra = .false.
nargs = cesr_iargc()

if (nargs > 0) then
  call cesr_getarg(1, lat_file)
  print *, 'Using ', trim(lat_file)
  print_extra = .true.
endif

call bmad_parser (lat_file, lat)

open (1, file = 'output.now')

!

local = floor_position_struct(r0_vec$, w_unit$, 0.0_rp, 0.0_rp, 0.0_rp)
branch => lat%branch(0)

do i = 1, lat%n_ele_max
  local%r = lat%particle_start%vec(1:5:2)

  if (branch%ele(i)%orientation == 1) then
    f0 = branch%ele(i-1)%floor
  else
    f0 = branch%ele(i)%floor
  endif

  floor = coords_local_curvilinear_to_floor(local, branch%ele(i), .true., calculate_angles = .true.)
  write (1, '(a, i0, a, 3f13.8, 2x, 3f13.8)') '"curvi-to-floor-up-T ', i, '" ABS 0 ', floor%r-f0%r, &
                                                       floor%theta-f0%theta, floor%phi-f0%phi, floor%psi-f0%psi
  floor = coords_local_curvilinear_to_floor(local, branch%ele(i), .false., calculate_angles = .true.)
  write (1, '(a, i0, a, 3f13.8, 2x, 3f13.8)') '"curvi-to-floor-up-F ', i, '" ABS 0 ', floor%r-f0%r, &
                                                       floor%theta-f0%theta, floor%phi-f0%phi, floor%psi-f0%psi

  if (branch%ele(i)%orientation == 1) then
    f0 = branch%ele(i)%floor
  else
    f0 = branch%ele(i-1)%floor
  endif

  local%r = lat%particle_start%vec(1:5:2)
  local%r(3) = local%r(3) + lat%ele(i)%value(l$)

  floor = coords_local_curvilinear_to_floor(local, branch%ele(i), .true., calculate_angles = .true.)
  write (1, '(a, i0, a, 3f13.8, 2x, 3f13.8)') '"curvi-to-floor-dn-T ', i, '" ABS 0 ', floor%r-f0%r, &
                                                       floor%theta-f0%theta, floor%phi-f0%phi, floor%psi-f0%psi
enddo

!

do i = 1, lat%n_ele_max
  floor = branch%ele(i-1)%floor
  floor2 = coords_relative_to_floor(floor, lat%particle_start%vec(1:5:2))
  floor%r = floor2%r
  local = coords_floor_to_curvilinear(floor, branch%ele(lat%n_ele_max), ele1, status)
  write (1, '(a, i0, a, i4, 3f13.8, 2x, 3f13.8, i6)') '"Floor-to-curvi-up ', i, '" ABS 0 ', ele1%ix_ele, local%r, &
                                                       local%theta, local%phi, local%psi, status
  floor = branch%ele(i)%floor
  floor2 = coords_relative_to_floor(floor, lat%particle_start%vec(1:5:2))
  floor%r = floor2%r
  local = coords_floor_to_curvilinear(floor, branch%ele(0), ele1, status)
  write (1, '(a, i0, a, i4, 3f13.8, 2x, 3f13.8, i6)') '"Floor-to-curvi-dn ', i, '" ABS 0 ', ele1%ix_ele, local%r, &
                                                       local%theta, local%phi, local%psi, status
enddo

!

do i = 1, lat%n_ele_max
  floor = branch%ele(i-1)%floor
  floor2 = coords_relative_to_floor(floor, lat%particle_start%vec(1:5:2))
  floor%r = floor2%r
  local = coords_floor_to_local_curvilinear(floor, branch%ele(i), status)
  write (1, '(a, i0, a, 3f13.8, 2x, 3f13.8, i6)') '"Floor-to-loc-up ', i, '" ABS 0 ', local%r, &
                                                       local%theta, local%phi, local%psi, status
  floor = branch%ele(i)%floor
  floor2 = coords_relative_to_floor(floor, lat%particle_start%vec(1:5:2))
  floor%r = floor2%r
  local = coords_floor_to_local_curvilinear(floor, branch%ele(i), status)
  write (1, '(a, i0, a, 3f13.8, 2x, 3f13.8, i6)') '"Floor-to-loc-dn ', i, '" ABS 0 ', local%r, &
                                                       local%theta, local%phi, local%psi, status
enddo

!

close (1)

end program
