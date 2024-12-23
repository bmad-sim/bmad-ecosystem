!+
! Regression tests for simulation utility routines
!-

program math_test

use bmad
use random_mod
use naff_mod
use cubic_interpolation_mod

implicit none

type (spline_struct) a_spline(6)
type (coord_struct) orbit
type (field_at_3D_box_struct) box3_field
type (field_at_2D_box_struct) box2_field
type (bicubic_coef_struct) bi_coef
type (tricubic_coef_struct) tri_coef

real(rp) array(4), dE, freq(3), m(6,6)
real(rp) sig1, sig2, sig3, quat(0:3), omega(3), axis2(3), angle2
real(rp) phi1, phi2, phi3, y1, dy1, y2, dy2, dx, dy, dz
real(rp) vec3(3), vec3a(3), vec3b(3), vec3c(3), axis(3), angle, w_mat(3,3), unit_mat(3,3)
real(rp) field2(10:12, 20:22), field3(10:12, 20:22, 30:32), ff, df_dx, df_dy, df_dz, ff0, ff1
real(rp) del, dff_dx, dff_dy, dff_dz, x, y, z, value
real(rp) xx(4), yy(4), zz(4), coef(3), coef_in(3)

integer i, j, k, ie, n, which, where, n_freq, mult, power, width, digits, species, s_mat(6,6)
integer, allocatable :: arr(:)

complex(rp) cdata(32)
complex(rp) amp(3), eval(6), evec(6,6)
complex(rp) amp1, amp2, amp3

logical match, ok, err

character(40) str, sub1, sub2, sub3
character(100) name
character(2) code
character(16) :: extrap(0:3) = [character(16):: 'ZERO', 'LINEAR', 'CONSTANT', 'SYMMETRIC']

!

open (1, file = 'output.now')

!

coef_in = [1.0_rp, 2.0_rp, 3.0_rp]
xx = [0.0_rp, 0.5_rp, 1.0_rp, 1.5_rp]
yy = [-1.0_rp, 0.0_rp, 1.0_rp, -0.5_rp]
zz = coef_in(1) * xx + coef_in(2) * yy + coef_in(3)
call linear_fit_2d(xx, yy, zz, coef)

write (1, '(a, 3es12.4)') '"linear_fit_2d" ABS 1E-12', coef-coef_in

! pointer_to_locations

call pointer_to_locations ('7 n4:n5,4  5:,  :4', arr, n, 3, 10, &
                              [character(3):: 'n3', 'n4', 'n5', 'n6', 'n7', 'n8', 'n9', 'n10'])
write (1, '(a, 99i3)') '"Locs1" ABS 0    ', n, (arr(i), i = 1, n)

! Eigen modes

m(1,:) = [-0.98545332, -0.17845915, -0.02950643, -0.00003723, -0.00128464,  0.00169891]
m(2,:) = [ 0.19413650, -0.97957165, -0.00719935,  0.00089198, -0.00589572, -0.03814434]
m(3,:) = [-0.00085633,  0.00006026, -0.87618978, -0.00841354,  0.00001679,  0.00264910]
m(4,:) = [ 0.01161089,  0.02967948, 26.62729191, -0.88564965,  0.00002272, -0.04043493]
m(5,:) = [-0.08861524, -0.00613454,  0.03485868, -0.00268734,  0.94781760, -8.79515336]
m(6,:) = [ 0.00530822, -0.00029320,  0.00010783, -0.00001597,  0.01159978,  0.94735510]

s_mat = 0
s_mat(1,2) =  1;  s_mat(3,4) =  1;  s_mat(5,6) =  1
s_mat(2,1) = -1;  s_mat(4,3) = -1;  s_mat(6,5) = -1

call mat_eigen (m, eval, evec, err)
do i = 1, 6
  write (1, '(a, 6(2x, 2f12.8))') '"Eigen-val-' // int_str(i) // '" ABS 1e-4', eval(i)
  write (1, '(a, 6(2x, 2f12.8))') '"Eigen-vec-' // int_str(i) // '" ABS 1e-4', evec(i,:)
enddo

write (1, '(a, 3(2x, 2f12.8))') '"vSv-odd"  ABS 1e-4', (sum(conjg(evec(i,:)) * matmul(s_mat, evec(i,:))), i = 1,5,2)
write (1, '(a, 3(2x, 2f12.8))') '"vSv-even" ABS 1e-4', (sum(conjg(evec(i,:)) * matmul(s_mat, evec(i,:))), i = 2,6,2)
do i = 1, 5, 2
  write (1, '(a, 3(2x, 2f12.8))') '"v-pairs-re' // int_str(i) // '"  ABS 1e-4', real(evec(i,:)-evec(i+1,:), rp)
  write (1, '(a, 3(2x, 2f12.8))') '"v-pairs-im' // int_str(i) // '"  ABS 1e-4', aimag(evec(i,:)+evec(i+1,:))
enddo

! random

call ran_seed_put (123)
call ran_gauss(array, sigma_cut = 0.5_rp)
write (1, '(a, 4f14.8)') '"Ran-Gauss-pseudo" ABS 1e-10 ', array
call ran_engine (set = 'quasi')
do i = 1, size(array)
  call ran_gauss(array(i), sigma_cut = 0.5_rp)
enddo
write (1, '(a, 4f14.8)') '"Ran-Gauss-quasi" ABS 1e-10 ', array

! anom mag moment

species = species_id('#3He--')
value = anomalous_moment_of (species)
write (1, '(a, es16.8)')  '"Anom-Moment-He3" REL 1E-6', value

! Tricubic interpolation

do i = 10, 12
do j = 20, 22
 field2(i,j) = i + j
do k = 30, 32
 field3(i,j,k) = i + j + k
enddo
enddo
enddo

del = 0.0001
x = 0.1
y = 0.2
z = 0.3

do i = 10, 12
do j = 20, 22
  ie = modulo(i + j, 4)
  call bicubic_compute_field_at_2D_box(field2, lbound(field2), i, j, extrap(ie), box2_field, err)
  call bicubic_interpolation_coefs(box2_field, bi_coef)
  ff = bicubic_eval(x, y, bi_coef, df_dx, df_dy)
  ff0 = bicubic_eval(x-del, y, bi_coef)
  ff1 = bicubic_eval(x+del, y, bi_coef)
  dff_dx = (ff1 - ff0) / (2 * del)
  ff0 = bicubic_eval(x, y-del, bi_coef)
  ff1 = bicubic_eval(x, y+del, bi_coef)
  dff_dy = (ff1 - ff0) / (2 * del)
  write (1, '(a, 2i3, a, 5f15.8, 4x, l1)') '"BiCubic', i, j, '" ABS 1E-10', ff, df_dx, df_dy
enddo
enddo

print *

do i = 9, 11
do j = 19, 21
do k = 29, 31
  ie = modulo(i + j, 4)
  call tricubic_compute_field_at_3D_box(field3, lbound(field3), i, j, k, extrap(ie), box3_field, err)
  call tricubic_interpolation_coefs(box3_field, tri_coef)
  ff = tricubic_eval(x, y, z, tri_coef, df_dx, df_dy, df_dz)
  ff0 = tricubic_eval(x-del, y, z, tri_coef)
  ff1 = tricubic_eval(x+del, y, z, tri_coef)
  dff_dx = (ff1 - ff0) / (2 * del)
  ff0 = tricubic_eval(x, y-del, z, tri_coef)
  ff1 = tricubic_eval(x, y+del, z, tri_coef)
  dff_dy = (ff1 - ff0) / (2 * del)
  ff0 = tricubic_eval(x, y, z-del, tri_coef)
  ff1 = tricubic_eval(x, y, z+del, tri_coef)
  dff_dz = (ff1 - ff0) / (2 * del)
  write (1, '(a, 2i3, a, 7f13.6, 4x, l1)') '"TriCubic', i, j, '" ABS 1E-10', ff, df_dx, df_dy, df_dz
enddo
enddo
enddo

! Akima spline test

do i = 1, 6
  a_spline(i)%x0 = i + 0.1 * i * i
  a_spline(i)%y0 = i * i
enddo

call spline_akima(a_spline, ok)

call spline_evaluate (a_spline, 1.3_rp, ok, y1, dy1)
call spline_akima_interpolate (a_spline%x0, a_spline%y0, 1.3_rp, ok, y2, dy2)
write (1, '(a, 4f14.8)') '"Spline1" ABS 1E-10 ', y1, dy1, y2 - y1, dy2 - dy1

call spline_evaluate (a_spline, 2.7_rp, ok, y1, dy1)
call spline_akima_interpolate (a_spline%x0, a_spline%y0, 2.7_rp, ok, y2, dy2)
write (1, '(a, 4f14.8)') '"Spline2" ABS 1E-10 ', y1, dy1, y2 - y1, dy2 - dy1

call spline_evaluate (a_spline, 5.0_rp, ok, y1, dy1)
call spline_akima_interpolate (a_spline%x0, a_spline%y0, 5.0_rp, ok, y2, dy2)
write (1, '(a, 4f14.8)') '"Spline3" ABS 1E-10 ', y1, dy1, y2 - y1, dy2 - dy1

call spline_evaluate (a_spline, 8.0_rp, ok, y1, dy1)
call spline_akima_interpolate (a_spline%x0, a_spline%y0, 8.0_rp, ok, y2, dy2)
write (1, '(a, 4f14.8)') '"Spline4" ABS 1E-10 ', y1, dy1, y2 - y1, dy2 - dy1

call spline_evaluate (a_spline, 9.6_rp, ok, y1, dy1)
call spline_akima_interpolate (a_spline%x0, a_spline%y0, 9.6_rp, ok, y2, dy2)
write (1, '(a, 4f14.8)') '"Spline5" ABS 1E-10 ', y1, dy1, y2 - y1, dy2 - dy1

! Parse fortran format tests

call parse_fortran_format('34pf12.3', mult, power, code, width, digits)
write (1, '(a, 2(i0, a), a, 2(i0, a))') '"34pf12.3" STR "', mult, '(', power, 'p', trim(code),  width, '.', digits, ')"'

call parse_fortran_format('x ', mult, power, code, width, digits)
write (1, '(a, 2(i0, a), a, 2(i0, a))') '"x" STR "', mult, '(', power, 'p', trim(code),  width, '.', digits, ')"'

call parse_fortran_format('(7(4pi17)) ', mult, power, code, width, digits)
write (1, '(a, 2(i0, a), a, 2(i0, a))') '"(7(4pi17))" STR "', mult, '(', power, 'p', trim(code),  width, '.', digits, ')"'

call parse_fortran_format('i17 i', mult, power, code, width, digits)
write (1, '(a, 2(i0, a), a, 2(i0, a))') '"BAD-FMT" STR "', mult, '(', power, 'p', trim(code),  width, '.', digits, ')"'

! rotation tests

call mat_make_unit (unit_mat)

w_mat = unit_mat
call rotate_mat(w_mat, x_axis$, 0.37_rp)
call rotate_mat(w_mat, x_axis$, -0.37_rp, .true.)
w_mat = w_mat - unit_mat
write (1, '(a, es11.3)') '"rot X-X" ABS 1E-16 ', maxval(abs(w_mat))

w_mat = unit_mat
call rotate_mat(w_mat, y_axis$, 0.37_rp)
call rotate_mat(w_mat, y_axis$, -0.37_rp, .true.)
w_mat = w_mat - unit_mat
write (1, '(a, es11.3)') '"rot Y-Y" ABS 1E-16 ', maxval(abs(w_mat))

w_mat = unit_mat
call rotate_mat(w_mat, z_axis$, 0.37_rp)
call rotate_mat(w_mat, z_axis$, -0.37_rp, .true.)
w_mat = w_mat - unit_mat
write (1, '(a, es11.3)') '"rot Z-Z" ABS 1E-16 ', maxval(abs(w_mat))

!

axis = [3, 4, 5] / sqrt(50.0_rp)
vec3 = [-2, 3, -4]
angle = 0.67

quat = axis_angle_to_quat(axis, angle)
omega = quat_to_omega(quat)
quat = omega_to_quat(omega)
call quat_to_axis_angle(quat, axis2, angle2)
write (1, '(a, 4es11.3)') '"axis-angle  " ABS 1E-14  ', axis2, angle2
write (1, '(a, 4es11.3)') '"daxis-dangle" ABS 1E-14  ', axis2-axis, angle2-angle

vec3a = rotate_vec_given_axis_angle (vec3, axis, angle)

call axis_angle_to_w_mat (axis, angle, w_mat)
vec3b = matmul(w_mat, vec3)

vec3c = quat_rotate(quat, vec3)

write (1, '(a, 3f11.6)') '"rot vecA" ABS 1E-14  ', vec3a
write (1, '(a, 3es10.2)') '"drot vecB" ABS 1E-14  ', vec3b - vec3a
write (1, '(a, 3es10.2)') '"drot vecC" ABS 1E-14  ', vec3c - vec3a

!

axis = [1.0_rp, 2.0_rp, 3.0_rp]
axis = axis / norm2(axis)
angle  = 0.1_rp
call axis_angle_to_w_mat(axis, angle, w_mat)
quat = w_mat_to_quat(w_mat)
write (1, '(a, 3es10.2)') '"dRot0" ABS 1E-14 ', matmul(w_mat, vec3) - quat_rotate(quat, vec3)

axis = [1.0_rp, 0.01_rp, 0.0_rp]
axis = axis / norm2(axis)
angle  = 3.1_rp
call axis_angle_to_w_mat(axis, angle, w_mat)
quat = w_mat_to_quat(w_mat)
write (1, '(a, 3es10.2)') '"dRot1" ABS 1E-14 ', matmul(w_mat, vec3) - quat_rotate(quat, vec3)

axis = [0.01_rp, 1.0_rp, 0.0_rp]
axis = axis / norm2(axis)
angle  = 3.1_rp
call axis_angle_to_w_mat(axis, angle, w_mat)
quat = w_mat_to_quat(w_mat)
write (1, '(a, 3es10.2)') '"dRot2" ABS 1E-14 ', matmul(w_mat, vec3) - quat_rotate(quat, vec3)

axis = [0.0_rp, 0.01_rp, 1.0_rp]
axis = axis / norm2(axis)
angle  = 3.1_rp
call axis_angle_to_w_mat(axis, angle, w_mat)
quat = w_mat_to_quat(w_mat)
write (1, '(a, 3es10.2)') '"dRot3" ABS 1E-14 ', matmul(w_mat, vec3) - quat_rotate(quat, vec3)

! naff test

amp1 = cmplx(1.8000,0.0000)
sig1 = 0.753262
amp2 = cmplx(0.3000,0.3000)
sig2 = 0.423594
amp3 = cmplx(0.01230,0.1545)
sig3 = 0.173

do i = 1, size(cdata)
  phi1 = twopi*(sig1*(i-1))
  phi2 = twopi*(sig2*(i-1))
  phi3 = twopi*(sig3*(i-1))
  cdata(i) = amp1*exp(cmplx(0.0d0,-phi1)) + amp2*exp(cmplx(0.0d0,-phi2)) + amp3*exp(cmplx(0.0d0,-phi3))
enddo

call naff (cdata, freq, amp)
write (1, '(a, 3es16.8)') '"naff-freq1" REL 2E-6   ', freq(1), real(amp(1)), aimag(amp(1))
write (1, '(a, 3es16.8)') '"naff-freq2" REL 2E-6   ', freq(2), real(amp(2)), aimag(amp(2))
write (1, '(a, 3es16.8)') '"naff-freq3" REL 3E-6   ', freq(3), real(amp(3)), aimag(amp(3))

! Random test

call ran_engine ('quasi')

do i = 1, 10
  call ran_uniform_vector (array)
enddo

write (1, '(a, 4es20.10)') '"QuasiRan" ABS  0', array

! 

write (1, *)

orbit%p0c = 1e6
orbit%vec = 0
orbit%vec(6) = 0.1
orbit%species = positron$

call convert_pc_to (orbit%p0c * (1 + orbit%vec(6)), positron$, beta = orbit%beta)
call apply_energy_kick (1d2, orbit, [0.0_rp, 0.0_rp])
write (1, '(a, 2es20.12)') '"apply_energy_kick:0" REL 1E-12  ', orbit%beta, orbit%vec(6)

call convert_pc_to (orbit%p0c * (1 + orbit%vec(6)), positron$, beta = orbit%beta)
call apply_energy_kick (1d6, orbit, [0.0_rp, 0.0_rp])
write (1, '(a, 2es20.12)') '"apply_energy_kick:1" REL 1E-12  ', orbit%beta, orbit%vec(6)

! Note: run.py script that runs this program defines "Test_EV" environment variable to be "ZZZ"

call fullfilename("$Test_EV/aa${Test_EV}zz$Test_EV", name, ok)
write (1, '(2a)') '"filename:1" STR ', quote(name)

!

close(1)

end program
