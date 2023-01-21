program geometry_test

use bmad

implicit none

type (lat_struct), target :: lat
type (coord_struct) orbit
type (floor_position_struct) local0, local, floor, pos_ele, floor0, floor2, f0, body
type (ele_struct), pointer :: ele, ele1
type (em_field_struct) f1, f2
type (branch_struct), pointer :: branch

real(rp) w_mat(3,3), s, w1_mat(3,3), w2_mat(3,3), w3_mat(3,3), dw(3,3), dw1(3,3), t(3,3)
real(rp) floor_ps(6)

integer i, status, nargs
logical print_extra
character(2) :: o_name(-1:1) = ['--', '??', '++']

character(100) lat_file

!

lat_file = 'geometry_test.bmad'
print_extra = .false.
nargs = command_argument_count()

if (nargs > 0) then
  call get_command_argument(1, lat_file)
  print *, 'Using ', trim(lat_file)
  print_extra = .true.
endif

call bmad_parser (lat_file, lat)

open (1, file = 'output.now')

!

local0 = floor_position_struct(lat%particle_start%vec(1:5:2), mat3_unit$, lat%particle_start%vec(2), &
                                                     lat%particle_start%vec(4), lat%particle_start%vec(6))
call floor_angles_to_w_mat (local0%theta, local0%phi, local0%psi, local0%w)

branch => lat%branch(0)

!

do i = 1, lat%n_ele_max
  ele => branch%ele(i)
  body = coords_local_curvilinear_to_body (local0, ele, w1_mat, .true.)
  local = coords_body_to_local(body, ele, w2_mat, .true.)
  dw = w2_mat - w1_mat
  dw1 = local%w - local%w
  write (1, '(4a, 3f13.8, 2x, 3f13.8, 2x, 2f13.8)') '"loc-body-loc', o_name(ele%orientation), trim(ele%name), &
                                '" ABS 1E-12', local%r-local0%r, amod(local%theta-local0%theta), amod(local%phi-local0%phi), &
                                amod(local%psi-local0%psi), maxval(abs(dw)), maxval(abs(dw1))
enddo
print *

do i = 1, lat%n_ele_max
  ele => branch%ele(i)
  floor = coords_local_curvilinear_to_floor (local0, ele, .false., w1_mat, calculate_angles = .true.)
  local = coords_floor_to_local_curvilinear (floor, ele, status, w2_mat)
  dw = w2_mat - w1_mat
  dw1 = local%w - local%w
  write (1, '(4a, 3f13.8, 2x, 3f13.8, 2x, 2f13.8)') '"loc-global-loc', o_name(ele%orientation), trim(ele%name), &
                                '" ABS 1E-12', local%r-local0%r, amod(local%theta-local0%theta), amod(local%phi-local0%phi), &
                                amod(local%psi-local0%psi), maxval(abs(dw)), maxval(abs(dw1))
enddo
print *


do i = 1, lat%n_ele_max
  ele => branch%ele(i)
  floor = coords_local_curvilinear_to_floor (local0, ele, .true., w1_mat, calculate_angles = .true.)
  local = coords_floor_to_local_curvilinear (floor, ele, status, w2_mat)
  local = coords_local_curvilinear_to_body (local, ele, w3_mat, .true.)
  dw = matmul(w2_mat, w3_mat) - w1_mat
  dw1 = local%w - local%w
  write (1, '(4a, 3f13.8, 2x, 3f13.8, 2x, 2f13.8)') '"body-global-body', o_name(ele%orientation), trim(ele%name), &
                                '" ABS 1E-12', local%r-local0%r, amod(local%theta-local0%theta), amod(local%phi-local0%phi), &
                                amod(local%psi-local0%psi), maxval(abs(dw)), maxval(abs(dw1))
enddo
print *

if (print_extra) stop

!

orbit%s = 2.001
orbit%vec = [0.01_rp, 0.01_rp, 0.02_rp, 0.2_rp, 0.03_rp, 2.0_rp]
floor_ps = orbit_to_floor_phase_space (orbit, lat%ele(2))
write (1, '(a, 3f15.11, 4x, 3f15.11)') '"floor_phase_space1" ABS 1E-12', floor_ps(1:5:2), floor_ps(2:6:2)

!

do i = 1, lat%n_ele_max
  ele => branch%ele(i)

  if (ele%orientation == 1) then
    f0 = branch%ele(i-1)%floor
  else
    f0 = ele%floor
  endif

  floor = coords_local_curvilinear_to_floor(local0, ele, .true., calculate_angles = .true.)
  write (1, '(4a, 3f13.8, 2x, 3f13.8)') '"curvi-to-floor-up-T', o_name(ele%orientation), trim(ele%name), '" ABS 0 ', &
                                        floor%r-f0%r, amod(floor%theta-f0%theta), amod(floor%phi-f0%phi), amod(floor%psi-f0%psi)
  floor = coords_local_curvilinear_to_floor(local0, ele, .false., calculate_angles = .true.)
  write (1, '(4a, 3f13.8, 2x, 3f13.8)') '"curvi-to-floor-up-F', o_name(ele%orientation), trim(ele%name), '" ABS 0 ', &
                                        floor%r-f0%r, amod(floor%theta-f0%theta), amod(floor%phi-f0%phi), amod(floor%psi-f0%psi)

  if (ele%orientation == 1) then
    f0 = ele%floor
  else
    f0 = branch%ele(i-1)%floor
  endif

  local = local0
  local%r(3) = local%r(3) + ele%value(l$)

  floor = coords_local_curvilinear_to_floor(local, ele, .true., calculate_angles = .true.)
  write (1, '(4a, 3f13.8, 2x, 3f13.8)') '"curvi-to-floor-dn-T', o_name(ele%orientation), trim(ele%name), '" ABS 0 ', &
                                              floor%r-f0%r, amod(floor%theta-f0%theta), amod(floor%phi-f0%phi), amod(floor%psi-f0%psi)
enddo

!

branch => lat%branch(1)

do i = 1, branch%n_ele_max
  ele => branch%ele(i)

  floor = branch%ele(i-1)%floor
  floor2 = coords_relative_to_floor(floor, lat%particle_start%vec(1:5:2))
  floor%r = floor2%r
  local = coords_floor_to_curvilinear(floor, branch%ele(branch%n_ele_max), ele1, status)
  write (1, '(4a, i4, 3f13.8, 2x, 3f13.8, i6)') '"Floor-to-curvi-up', o_name(ele%orientation), trim(ele%name), '" ABS 0 ', &
                                                       ele1%ix_ele, local%r, local%theta, local%phi, local%psi, status
  floor = ele%floor
  floor2 = coords_relative_to_floor(floor, lat%particle_start%vec(1:5:2))
  floor%r = floor2%r
  local = coords_floor_to_curvilinear(floor, branch%ele(0), ele1, status)
  write (1, '(4a, i4, 3f13.8, 2x, 3f13.8, i6)') '"Floor-to-curvi-dn', o_name(ele%orientation), trim(ele%name), '" ABS 0 ', &
                                                       ele1%ix_ele, local%r, local%theta, local%phi, local%psi, status
enddo

!

do i = 1, branch%n_ele_max
  ele => branch%ele(i)

  floor = branch%ele(i-1)%floor
  floor2 = coords_relative_to_floor(floor, lat%particle_start%vec(1:5:2))
  floor%r = floor2%r
  local = coords_floor_to_local_curvilinear(floor, ele, status)
  write (1, '(4a, 3f13.8, 2x, 3f13.8, i6)') '"Floor-to-loc-up', o_name(ele%orientation), trim(ele%name), '" ABS 0 ', &
                                                  local%r, local%theta, local%phi, local%psi, status
  floor = ele%floor
  floor2 = coords_relative_to_floor(floor, lat%particle_start%vec(1:5:2))
  floor%r = floor2%r
  local = coords_floor_to_local_curvilinear(floor, ele, status)
  write (1, '(4a, 3f13.8, 2x, 3f13.8, i6)') '"Floor-to-loc-dn', o_name(ele%orientation), trim(ele%name), '" ABS 0 ', &
                                                  local%r, local%theta, local%phi, local%psi, status
enddo

!

close (1)

!---------------------------
contains

function amod(angle) result (mod_angle)
real(rp) angle, mod_angle
mod_angle = modulo2(angle, pi)
end function amod

end program
